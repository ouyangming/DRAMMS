/**
 * @file  ImageOperation.cxx
 * @brief This program performs basic image operations.
 *
 * Note that the input to this program can be any type of image supported by
 * DRAMMS, including, in particular, deformation fields.
 *
 * Copyright (c) 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <iostream>
#include <vector>
#include <memory>
#include <limits>
#include "common/imageio.h"
#include "common/utilities.h"
#include "common/general.h"   // split()
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dramms/basis.h>


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;



// ===========================================================================
// operations
// ===========================================================================

// ---------------------------------------------------------------------------
void PrintImageHeader(const Image* image)
{
    Image::Header hdr = image->hdr;
    if (image->imgfmt == Image::FORMAT_DRAMMS) hdr = PermuteXY(image->hdr);
    else                                       hdr = image->hdr;

    cout << "Note: Fields dim0, dim4 and dim5 may not reflect the actual values" << endl;
    cout << "      of the image header as they are adjusted by the image reader" << endl;
    cout << "      to not be equal to zero." << endl;
    cout << endl;

    cout.precision(10);
    for (int i = 0; i < 7; i++) {
        cout << "dim" << i << "           " << hdr.dim[i] << endl;
    }
    cout << "datatype       " << hdr.datatype << endl;
    cout << "bitpix         " << hdr.bitpix << endl;
    for (int i = 0; i < 7; i++) {
        cout << "pixdim" << i << "        " << hdr.pixdim[i] << endl;
    }
    cout << "scl_slope      " << hdr.scl_slope << endl;
    cout << "scl_inter      " << hdr.scl_inter << endl;
    cout << "intent_code    " << hdr.intent_code << endl;
    cout << "intent_p1      " << hdr.intent_p1 << endl;
    cout << "intent_p2      " << hdr.intent_p2 << endl;
    cout << "intent_p3      " << hdr.intent_p3 << endl;
    cout << "qform_code     " << hdr.qform_code << endl;
    Image::Transform qto_xyz = GetQFormTransform(hdr, Image::FORMAT_ITK);
    for (int r = 0; r < 4; r++) {
        cout << "qto_xyz:" << (r + 1) << "  ";
        for (int c = 0; c < 4; c++) {
            if (signbit(qto_xyz.m[r][c]) == 0) cout << " ";
            cout << "   " << fixed << qto_xyz.m[r][c];
        }
        cout << endl;
    }
    cout << "sform_code     " << hdr.sform_code << endl;
    Image::Transform sto_xyz = GetSFormTransform(hdr, Image::FORMAT_ITK);
    cout << "sto_xyz:1" << "  ";
    for (int c = 0; c < 4; c++) {
        if (signbit(hdr.srow_x[c]) == 0) cout << " ";
        cout << "   " << fixed << hdr.srow_x[c];
    }
    cout << endl;
    cout << "sto_xyz:2" << "  ";
    for (int c = 0; c < 4; c++) {
        if (signbit(hdr.srow_y[c]) == 0) cout << " ";
        cout << "   " << fixed << hdr.srow_y[c];
    }
    cout << endl;
    cout << "sto_xyz:3" << "  ";
    for (int c = 0; c < 4; c++) {
        if (signbit(hdr.srow_z[c]) == 0) cout << " ";
        cout << "   " << fixed << hdr.srow_z[c];
    }
    cout << endl;
    int x_orient, y_orient, z_orient;
    GetImageOrientation(hdr, x_orient, y_orient, z_orient, Image::FORMAT_ITK);
    cout << "qform_orient   " << OrientationToString(x_orient, y_orient, z_orient) << endl;
    cout << "lr_order       " << LROrderToString(GetLeftRightOrder(image->hdr, Image::FORMAT_ITK)) << endl;
    cout << "descrip        " << image->hdr.descrip << endl;
    cout << "magic          " << image->hdr.magic << endl;
}

// ===========================================================================
// command-line parsing
// ===========================================================================

// ---------------------------------------------------------------------------
class RecordOperationVisitor : public TCLAP::Visitor
{
public:
    RecordOperationVisitor() : _op('\0'), _ops(NULL) {}
    RecordOperationVisitor(char op, vector<char>* ops) : _op(op), _ops(ops) {}

    void visit()
    {
        assert(_ops != NULL && _op != '\0');
        if (_ops && _op != '\0') _ops->push_back(_op);
    }

protected:
    char          _op;  ///< Identification of operation to perform.
    vector<char>* _ops; ///< Operations to perform.
}; // class RecordOperationVisitor

// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    bool ok = true; // program execution so far successful

    vector<char> ops; // operations to perform

    // -----------------------------------------------------------------------
    // operation visitors
    const string operations = "acdemprsxyzXYZ"; // one unique character for each
                                                // implemented operation

    map<char, RecordOperationVisitor> record_operation_visitor;
    for (string::const_iterator op = operations.begin(); op != operations.end(); ++op) {
        record_operation_visitor[*op] = RecordOperationVisitor(*op, &ops);
    }

    // -----------------------------------------------------------------------
    // argument constraints
    PositiveValueConstraint<int> pos_int_constraint("<int>");

    // -----------------------------------------------------------------------
    // define command-line arguments
    StringArg type("t", "type",
            "Output datatype. Valid datatype identifiers are 'byte' or 'uchar',"
            " 'short', 'ushort', 'int', 'uint', and 'float'. If not specified,"
            " the datatype of the input image is used.",
            false, "DT_UNKNOWN", "<datatype>");

    MultiStringArg downsample("r", "downsample",
            "Downsample image using an integer valued downsample ratio.",
            false, "<int>", 1, false,
            &record_operation_visitor['r']);

    MultiStringArg extract("e", "extract",
            "Extract voxels inside foreground region specified by mask."
            " The mask can be any scalar image where voxels with an intensity value"
            " greater than zero are considered foreground.",
            false, "<file>", 1, false,
            &record_operation_visitor['e']);

    MultiStringArg mask("m", "mask",
            "Eliminate image values outside of foreground region specified by mask."
            " The mask can be any scalar image where voxels with an intensity value"
            " greater than zero are considered foreground.",
            false, "<file>", 1, false,
            &record_operation_visitor['m']);

    MultiStringArg resample("p", "pixdim",
            "Resample image on image grid with the specified voxel size but the"
            " (approximately) same extent as the current image grid, i.e., the"
            " product of number of voxels and voxel size remains (approximately)"
            " the same as before. Note that due to the integral nature of the"
            " number of voxels and the fact that the voxel size is being fixed,"
            " the image extent may be slightly different due to rounding errors."
            " In particular, the number of voxels are rounded up to the next"
            " integer value to ensure that no image information is lost.",
            false, "<sx>[,<sy>[,<sz>]]", 1, false,
            &record_operation_visitor['p']);

    MultiStringArg resize("d", "dim",
            "Resize image such that the resulting image has the specified number of"
            " voxels, padding and cropping it where necessary. Note that this operation"
            " does not change the voxel size. The --pixdim operation can be used"
            " instead to change the voxel size. As --pixdim also changes the number of"
            " voxels, it has to be given before --dim in order to resample the image"
            " on an image grid with the specified size and resolution.",
            false, "<nx>[,<ny>[,<nz>]]", 1, false,
            &record_operation_visitor['d']);

    MultiStringArg smooth("s", "smooth",
            "Smooth image using Gaussian kernel. The arguments to this operation"
            " are the radii and the sigma values of the kernel in each dimension"
            " or 'auto' if the kernel parameters should be determined automatically."
            " If the argument 'avg' is given, on the other side, the kernel parameters"
            " are choosen automatically and the input image is averaged with the"
            " smoothed image in order to preserve edges. If only two values are given"
            " for each parameter in case of a three-dimensional image,"
            " ry := rx, and sy := sx.",
            false, "<rx>[,<ry>[,<rz>]],<sx>[,<sy>[,<sz>]]", 1, false,
            &record_operation_visitor['s']);

    MultiStringArg threshold("c", "clamp",
            "Threshold image values using the specified threshold(s)."
            " If only one threshold is specified, all voxel components are"
            " thresholded by this value, i.e., components whose absolute value"
            " exceed the absolute value of this threshold are set to the threshold"
            " value without changing their sign. This process is also commonly"
            " referred to as clamping. In case of a vector field, multiple"
            " thresholds, one for each vector dimension, can be specified.",
            false, "<float>[,<float>]...", 1, false,
            &record_operation_visitor['c']);

    MultiStringArg extract_x("x", "extract-x",
            "Extract voxels within specified range of x indices.",
            false, "<int>[,<int>]", 1, false,
            &record_operation_visitor['x']);

    MultiStringArg extract_y("y", "extract-y",
            "Extract voxels within specified range of y indices.",
            false, "<int>[,<int>]", 1, false,
            &record_operation_visitor['y']);

    MultiStringArg extract_z("z", "extract-z",
            "Extract voxels within specified range of z indices.",
            false, "<int>[,<int>]", 1, false,
            &record_operation_visitor['z']);

    MultiStringArg zero_x("X", "zero-x",
            "Set voxel values within specified range of x indices to zero.",
            false, "<int>[,<int>]", 1, false,
            &record_operation_visitor['X']);

    MultiStringArg zero_y("Y", "zero-y",
            "Set voxel values within specified range of y indices to zero.",
            false, "<int>[,<int>]", 1, false,
            &record_operation_visitor['Y']);

    MultiStringArg zero_z("Z", "zero-z",
            "Set voxel values within specified range of z indices to zero.",
            false, "<int>[,<int>]", 1, false,
            &record_operation_visitor['Z']);

    PositionalArg input_file("input_file",
            "Input image file. The input image can either be a scalar intensity image"
            " or a vector field such as in particular a deformation field.",
            true, "", "<input_file>");

    PositionalArg output_file("output_file",
            "Output image file.",
            false, "", "<output_file>");

    MultiSwitchArg verbose("v", "verbose", "Increase verbosity of output messages.");

    // -----------------------------------------------------------------------
    // parse command-line
    try {
        vector<string> examples;

        examples.push_back("EXENAME image.nii.gz"
                           "\nPrints image header information.");

        examples.push_back("EXENAME -s 2,2,1.0,1.0 image.nii"
                           "\nSmooths the two-dimensional image using a Gaussian"
                           " kernel with radius 2 in both dimensions and sigma 1.0.");

        examples.push_back("EXENAME -s auto image.nii image_smoothed.nii.gz"
                           "\nSmooths the image using a Gaussian kernel with automatically"
                           " determined parameters.");

        CmdLine cmd(
                // program identification
                "ImageOperation", PROJECT,
                // program description
                "This program performs certain operations on a given image, either an"
                " intensity image or a deformation field.",
                // example usage
                examples,
                // version information
                RELEASE, "2012 University of Pennsylvania");

        // the constructor of the CmdLine class has already added the standard
        // options for help output and version information

        cmd.add(type);
        cmd.add(downsample);
        cmd.add(resample);
        cmd.add(resize);
        cmd.add(extract);
        cmd.add(extract_x);
        cmd.add(extract_y);
        cmd.add(extract_z);
        cmd.add(mask);
        cmd.add(smooth);
        cmd.add(threshold);
        cmd.add(zero_x);
        cmd.add(zero_y);
        cmd.add(zero_z);
        cmd.add(verbose);
        cmd.add(input_file);
        cmd.add(output_file);

        if (argc == 1) {
            cmd.print_help();
            exit(EXIT_FAILURE);
        }

        cmd.parse(argc, argv);

    } catch (CmdLineException& e) {
        // invalid command-line specification
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    string datatype = type.getValue();

    if      (datatype == "byte" || datatype == "uchar") datatype = "DT_UNSIGNED_CHAR";
    else if (datatype == "short" ) datatype = "DT_INT16";
    else if (datatype == "ushort") datatype = "DT_UINT16";
    else if (datatype == "int"   ) datatype = "DT_INT32";
    else if (datatype == "uint"  ) datatype = "DT_UINT32";
    else if (datatype == "float" ) datatype = "DT_FLOAT";
    
    short dt = static_cast<short>(nifti_datatype_from_string(datatype.c_str()));
    if (datatype != "DT_UNKNOWN" && dt == DT_UNKNOWN) {
        cerr << "Invalid/Unsupported datatype: " << type.getValue() << endl;
        exit(EXIT_FAILURE);
    }

    // -----------------------------------------------------------------------
    // read image

    // use ITK format such that user input can be given in the x and y order
    // as expected by normal users, i.e., the NIfTI order
    Image* image = ReadImage(input_file.getValue().c_str(), Image::FORMAT_ITK);
    if (image == NULL) {
        cerr << "Failed to read image from file " << input_file.getValue() << endl;
        exit(EXIT_FAILURE);
    }

    // -----------------------------------------------------------------------
    // perform operations
    try {
        // if no operation was specified, simply print image header information
        if (ops.empty()) {
            PrintImageHeader(image);
        } else {
            
            if (output_file.getValue().empty()) {
                cerr << "No output image file name specified!" << endl;
                exit(EXIT_FAILURE);
            }

            // ---------------------------------------------------------------
            // cast to output datatype
            if (dt != DT_UNKNOWN && image->hdr.datatype != dt) {
                Image* casted_image = CastImage(image, dt, true, false);
                delete image;
                image = casted_image;
                if (casted_image == NULL) {
                    cerr << "Failed to cast input image to requested output datatype!" << endl;
                    ok = false;
                }
            }

            map<char, int> argidx; // index of argument for given operation
            for (vector<char>::const_iterator op = ops.begin(); ok && op != ops.end(); ++op) {
                if (argidx.find(*op) == argidx.end()) argidx[*op] = 0;
                switch (*op) {
                    // -------------------------------------------------------
                    // downsample image by integral ratio
                    case 'r':
                        {
                            string arg = downsample.getValue()[argidx[*op]];
                            int ratio = 0;
                            if (sscanf(arg.c_str(), "%d", &ratio) != 1) {
                                cerr << "Invalid argument for operation -" << *op << ": " << arg << endl;
                                ok = false;
                                break;
                            }
                            Image* downsampled_image = NULL;
                            try {
                                bool _smooth = false;
                                if (op != ops.begin()) {
                                    vector<char>::const_iterator prev_op = op; prev_op--;
                                    if (*prev_op == 's') {
                                        string smooth_arg = smooth.getValue()[argidx['s'] - 1];
                                        _smooth = (smooth_arg == "auto");
                                    }
                                }
                                if (_smooth) {
                                    downsampled_image = SmoothAndDownsampleImage(image, ratio);
                                } else {
                                    downsampled_image = DownsampleImage(image, ratio);
                                }
                                if (downsampled_image == NULL) {
                                    cerr << "Failed to downsample image!" << endl;
                                    ok = false;
                                }
                            } catch (const exception& e) {
                                cerr << "Failed to downsample image! " << e.what() << endl;
                                ok = false;
                            }
                            if (!ok) break;
                            delete image;
                            image = downsampled_image;
                        }
                        break;
                    // -------------------------------------------------------
                    // resample image, i.e., change voxel size
                    case 'p':
                        {
                            Fvector3d pixdim;
                            pixdim.x = pixdim.y = pixdim.z = 0.0f;

                            string arg = resample.getValue()[argidx[*op]];
                            int n = sscanf(arg.c_str(), "%f,%f,%f", &pixdim.x, &pixdim.y, &pixdim.z);
                            if (n < 1) {
                                cerr << "Invalid argument for operation -" << *op << ": " << arg << endl;
                                ok = false;
                                break;
                            } else if (n == 1) {
                                pixdim.y = pixdim.x;
                                pixdim.z = pixdim.x;
                            } else {
                                if (pixdim.y <= 0) pixdim.y = image->hdr.pixdim[2];
                                if (pixdim.z <= 0) pixdim.z = image->hdr.pixdim[3];
                            }

                            Image* resampled_image = ResampleImage(image, pixdim);
                            if (resampled_image == NULL) {
                                cerr << "Failed to change voxel size of image to ";
                                cerr << pixdim.x << "x" << pixdim.y;
                                if (image->region.nz > 1) cerr << "x" << pixdim.z;
                                cerr << " mm^3!" << endl;
                                ok = false;
                                break;
                            }
                            delete image;
                            image = resampled_image;
                        }
                        break;
                    // -------------------------------------------------------
                    // resize image, i.e., change number of voxels
                    case 'd':
                        {
                            Image::Region region;

                            string arg = resize.getValue()[argidx[*op]];
                            int n = sscanf(arg.c_str(), "%d,%d,%d", &region.nx, &region.ny, &region.nz);
                            if (n < 1) {
                                cerr << "Invalid argument for operation -" << *op << ": " << arg << endl;
                                ok = false;
                                break;
                            } else if (n == 1) {
                                region.ny = region.nx;
                                if (image->region.nz > 1) region.nz = region.nx;
                            } else {
                                if (region.ny <= 0) region.ny = image->region.ny;
                                if (region.nz <= 0) region.nz = image->region.nz;
                            }

                            region.ox = static_cast<int>(round((image->region.nx - region.nx) / 2.0));
                            region.oy = static_cast<int>(round((image->region.ny - region.ny) / 2.0));
                            region.oz = static_cast<int>(round((image->region.nz - region.nz) / 2.0));

                            Image* resized_image = ResizeImage(image, region);
                            if (resized_image == NULL) {
                                cerr << "Failed to change size of image to ";
                                cerr << region.nx << "x" << region.ny;
                                if (image->region.nz > 1) cerr << "x" << region.nz;
                                cerr << " voxels!" << endl;
                                ok = false;
                                break;
                            }
                            delete image;
                            image = resized_image;
                        }
                        break;
                    // -------------------------------------------------------
                    // smooth image components
                    case 's':
                        {
                            string arg = smooth.getValue()[argidx[*op]];

                            // In case of the -r option, we use a different smoothing kernel.
                            // Thus skip the smoothing here, but perform it when the
                            // -r operation is performed.
                            vector<char>::const_iterator next_op = op; next_op++;
                            if (next_op == ops.end() || *next_op != 'r' || arg != "auto") {

                                Image* gauss = NULL;
                                if (image->region.nz <= 1) {
                                    Ivector2d radius;
                                    Fvector2d sigma;
                                    if (arg == "auto") {
                                        GetGaussParameters(image->region.nx, radius.x, sigma.x);
                                        GetGaussParameters(image->region.ny, radius.y, sigma.y);
                                    } else if (arg == "avg") {
                                        Fvector2d pixdim;
                                        pixdim.x = image->hdr.pixdim[1];
                                        pixdim.y = image->hdr.pixdim[2];
                                        GetGaussParameters(pixdim, radius, sigma);
                                    } else {
                                        vector<string> args = split(arg, ',');
                                        if (args.size() == 2) {
                                            ok =       (sscanf(args[0].c_str(), "%d", &radius.x) == 1);
                                            ok = ok && (sscanf(args[1].c_str(), "%f", &sigma .x) == 1);
                                            radius.y = radius.x;
                                            sigma .y = radius.y;
                                        } else if (args.size() == 4) {
                                            ok =       (sscanf(args[0].c_str(), "%d", &radius.x) == 1);
                                            ok = ok && (sscanf(args[1].c_str(), "%d", &radius.y) == 1);
                                            ok = ok && (sscanf(args[2].c_str(), "%f", &sigma .x) == 1);
                                            ok = ok && (sscanf(args[3].c_str(), "%f", &sigma .y) == 1);
                                        } else {
                                            cerr << "Invalid argument for operation -" << *op << ": " << arg << endl;
                                            ok = false;
                                            break;
                                        }
                                    }
                                    gauss = Gauss(radius, sigma, image->imgfmt);
                                } else {
                                    Ivector3d radius;
                                    Fvector3d sigma;
                                    if (arg == "auto") {
                                        GetGaussParameters(image->region.nx, radius.x, sigma.x);
                                        GetGaussParameters(image->region.ny, radius.y, sigma.y);
                                        GetGaussParameters(image->region.nz, radius.z, sigma.z);
                                    } else if (arg == "avg") {
                                        Fvector3d pixdim;
                                        pixdim.x = image->hdr.pixdim[1];
                                        pixdim.y = image->hdr.pixdim[2];
                                        pixdim.z = image->hdr.pixdim[3];
                                        GetGaussParameters(pixdim, radius, sigma);
                                    } else {
                                        vector<string> args = split(arg, ',');
                                        if (args.size() == 2) {
                                            ok =       (sscanf(args[0].c_str(), "%d", &radius.x) == 1);
                                            ok = ok && (sscanf(args[1].c_str(), "%f", &sigma .x) == 1);
                                            radius.z = radius.y = radius.x;
                                            sigma .z = sigma .y = radius.y;
                                        } else if (args.size() == 4) {
                                            ok =       (sscanf(args[0].c_str(), "%d", &radius.x) == 1);
                                            ok = ok && (sscanf(args[1].c_str(), "%d", &radius.z) == 1);
                                            ok = ok && (sscanf(args[2].c_str(), "%f", &sigma .x) == 1);
                                            ok = ok && (sscanf(args[3].c_str(), "%f", &sigma .z) == 1);
                                            radius.y = radius.x;
                                            sigma .y = sigma .x;
                                        } else if (args.size() == 6) {
                                            ok =       (sscanf(args[0].c_str(), "%d", &radius.x) == 1);
                                            ok = ok && (sscanf(args[1].c_str(), "%d", &radius.y) == 1);
                                            ok = ok && (sscanf(args[2].c_str(), "%d", &radius.z) == 1);
                                            ok = ok && (sscanf(args[3].c_str(), "%f", &sigma .x) == 1);
                                            ok = ok && (sscanf(args[4].c_str(), "%f", &sigma .y) == 1);
                                            ok = ok && (sscanf(args[5].c_str(), "%f", &sigma .z) == 1);
                                        } else {
                                            cerr << "Invalid argument for operation -" << *op << ": " << arg << endl;
                                            ok = false;
                                            break;
                                        }
                                    }
                                    gauss = Gauss(radius, sigma, image->imgfmt);
                                }

                                if (gauss == NULL) {
                                    cerr << "Failed to create smoothing kernel!" << endl;
                                    ok = false;
                                    break;
                                }

                                Image* smoothed_image = new Image(image->region.nx,
                                                                  image->region.ny,
                                                                  image->region.nz,
                                                                  image->hdr.datatype,
                                                                  image->GetNumberOfComponents(),
                                                                  image->imgfmt);
                                if (smoothed_image == NULL) {
                                    cerr << "Failed to allocate memory for smoothed image!" << endl;
                                    ok = false;
                                    delete gauss;
                                    break;
                                }
                                smoothed_image->CopyRegion     (image);
                                smoothed_image->CopyTransform  (image);
                                smoothed_image->CopyDataScaling(image);
                                smoothed_image->CopyMetaData   (image);

                                ConvolveImage(image, gauss, smoothed_image);

                                delete gauss;
                                if (arg == "avg") {
                                    for (int k = 0; k < image->region.nz; k++) {
                                        for (int j = 0; j < image->region.ny; j++) {
                                            for (int i = 0; i < image->region.nx; i++) {
                                                for (int n = 0; n < image->GetNumberOfComponents(); n++) {
                                                    image->set(i, j, k, n, 0.5 * (image->get(i, j, k, n) + smoothed_image->get(i, j, k, n)));
                                                }
                                            }
                                        }
                                    }
                                    delete smoothed_image;
                                } else {
                                    delete image;
                                    image = smoothed_image;
                                }
                            }
                        }
                        break;
                    // -------------------------------------------------------
                    // extract voxels inside ROI
                    case 'e':
                    case 'x':
                    case 'y':
                    case 'z':
                        {
                            Image::Region region = image->region;

                            if (*op == 'e') {
                                string mask_file  = extract.getValue()[argidx[*op]];
                                Image* mask_image = ReadImage(mask_file.c_str(), image->imgfmt);
                                if (mask_image == NULL) {
                                    cerr << "Failed to read mask from file " << mask_file << "!" << endl;
                                    ok = false;
                                    break;
                                }
                                region = GetForegroundRegion(mask_image);
                                delete mask_image;
                            } else {
                                string arg;
                                if      (*op == 'x') arg = extract_x.getValue()[argidx[*op]];
                                else if (*op == 'y') arg = extract_y.getValue()[argidx[*op]];
                                else if (*op == 'z') arg = extract_z.getValue()[argidx[*op]];
                                int li = 0, ui = 0;
                                int n = sscanf(arg.c_str(), "%d,%d", &li, &ui);
                                if (n < 1) {
                                    cerr << "Invalid argument for -" << *op << ": " << arg << endl;
                                    ok = false;
                                    break;
                                } else if (n == 1) {
                                    ui = li;
                                }
                                if (*op == 'x') {
                                    region.ox = li;
                                    region.nx = li <= ui ? (ui - li + 1) : (li - ui + 1);
                                } else if (*op == 'y') {
                                    region.oy = li;
                                    region.ny = li <= ui ? (ui - li + 1) : (li - ui + 1);
                                } else if (*op == 'z') {
                                    region.oz = li;
                                    region.nz = li <= ui ? (ui - li + 1) : (li - ui + 1);
                                }
                            }

                            if (verbose.getValue() > 0) {
                                cout << "Extracting voxels within " << region << endl;
                            }

                            Image* tmp = ResizeImage(image, region);
                            if (tmp == NULL) {
                                cerr << "Failed to allocate memory for temporary image!" << endl;
                                ok = false;
                                break;
                            }
                            delete image;
                            image = tmp;
                        }
                        break;
                    // -------------------------------------------------------
                    // set values outside ROI to zero
                    case 'm':
                        {
                            string mask_file  = mask.getValue()[argidx[*op]];
                            Image* mask_image = ReadImage(mask_file.c_str(), image->imgfmt);
                            if (mask_image == NULL) {
                                cerr << "Failed to read mask from file " << mask_file << "!" << endl;
                                ok = false;
                                break;
                            }
                            Image::Region region = GetForegroundRegion(mask_image);
                            delete mask_image;

                            if (verbose.getValue() > 0) {
                                cout << "Setting values outside " << region << " to zero." << endl;
                            }

                            vector<float> zero(image->GetNumberOfComponents(), 0.0f);
                            for (int k = 0; k < region.oz; k++) {
                                for (int j = 0; j < image->region.ny; j++) {
                                    for (int i = 0; i < image->region.nx; i++) {
                                        image->set(i, j, k, zero);
                                    }
                                }
                            }
                            for (int k = region.oz + region.nz; k < image->region.nz; k++) {
                                for (int j = 0; j < image->region.ny; j++) {
                                    for (int i = 0; i < image->region.nx; i++) {
                                        image->set(i, j, k, zero);
                                    }
                                }
                            }
                            for (int j = 0; j < region.oy; j++) {
                                for (int k = 0; k < image->region.nz; k++) {
                                    for (int i = 0; i < image->region.nx; i++) {
                                        image->set(i, j, k, zero);
                                    }
                                }
                            }
                            for (int j = region.oy + region.ny; j < image->region.ny; j++) {
                                for (int k = 0; k < image->region.nz; k++) {
                                    for (int i = 0; i < image->region.nx; i++) {
                                        image->set(i, j, k, zero);
                                    }
                                }
                            }
                            for (int i = 0; i < region.ox; i++) {
                                for (int j = 0; j < image->region.ny; j++) {
                                    for (int k = 0; k < image->region.nz; k++) {
                                        image->set(i, j, k, zero);
                                    }
                                }
                            }
                            for (int i = region.ox + region.nx; i < image->region.nx; i++) {
                                for (int j = 0; j < image->region.ny; j++) {
                                    for (int k = 0; k < image->region.nx; k++) {
                                        image->set(i, j, k, zero);
                                    }
                                }
                            }
                        }
                        break;
                    // -------------------------------------------------------
                    // set values inside ROI to zero
                    case 'X':
                    case 'Y':
                    case 'Z':
                        {
                            Image::Region region = image->region;

                            string arg;
                            if      (*op == 'X') arg = zero_x.getValue()[argidx[*op]];
                            else if (*op == 'Y') arg = zero_y.getValue()[argidx[*op]];
                            else if (*op == 'Z') arg = zero_z.getValue()[argidx[*op]];
                            int li = 0, ui = 0;
                            int n = sscanf(arg.c_str(), "%d,%d", &li, &ui);
                            if (n < 1) {
                                cerr << "Invalid argument for -" << *op << ": " << arg << endl;
                                ok = false;
                                break;
                            } else if (n == 1) {
                                ui = li;
                            }
                            if (*op == 'X') {
                                region.ox = li;
                                region.nx = li <= ui ? (ui - li + 1) : (li - ui + 1);
                            } else if (*op == 'Y') {
                                region.oy = li;
                                region.ny = li <= ui ? (ui - li + 1) : (li - ui + 1);
                            } else if (*op == 'Z') {
                                region.oz = li;
                                region.nz = li <= ui ? (ui - li + 1) : (li - ui + 1);
                            }

                            if (verbose.getValue() > 0) {
                                cout << "Setting values within " << region << " to zero." << endl;
                            }

                            vector<float> zero(image->GetNumberOfComponents(), 0.0f);
                            for (int k = region.oz; k < region.oz + region.nz; k++) {
                                for (int j = region.oy; j < region.oy + region.ny; j++) {
                                    for (int i = region.ox; i < region.ox + region.nx; i++) {
                                        if (0 <= i && i < image->region.nx &&
                                                0 <= j && j < image->region.ny &&
                                                0 <= k && k < image->region.nz) {
                                            image->set(i, j, k, zero);
                                        }
                                    }
                                }
                            }
                        }
                        break;
                    // -------------------------------------------------------
                    // threshold image components
                    case 'c':
                        {
                            vector<float> thresholds;
                            string arg = threshold.getValue()[argidx[*op]];
                            string::size_type pos = 0;
                            for (;;) {
                                thresholds.resize(thresholds.size() + 1);
                                string::size_type n = arg.find(',', pos);
                                if (n < string::npos) n = n - pos;
                                if (sscanf(arg.substr(pos, n).c_str(), "%f", &thresholds[thresholds.size() - 1]) != 1) {
                                    thresholds.clear(); // indicates parse error
                                    break;
                                }
                                if (n == string::npos) break;
                                pos += n + 1;
                            }
                            if (thresholds.size() == 0) {
                                cerr << "Invalid argument for operation -" << *op << ": " << arg << endl;
                                ok = false;
                                break;
                            } else if (thresholds.size() == 1) {
                                while (thresholds.size() < image->GetNumberOfComponents()) {
                                    thresholds.push_back(thresholds[0]);
                                }
                            } else {
                                while (thresholds.size() < image->GetNumberOfComponents()) {
                                    thresholds.push_back(numeric_limits<float>::max());
                                }
                            }
                            vector<float> v;
                            for (int k = 0; k < image->region.nz; k++) {
                                for (int j = 0; j < image->region.ny; j++) {
                                    for (int i = 0; i < image->region.nx; i++) {
                                        image->get(i, j, k, v);
                                        for (int n = 0; n < image->GetNumberOfComponents(); n++) {
                                            // if (fabs(v[n]) > thresholds[n]) {
                                                // if (v[n] < 0) v[n] = - thresholds[n];
                                                // else          v[n] = + thresholds[n];
											if (fabs(v[n]) < thresholds[n]) {
                                                v[n] = 0;
                                            }
                                        }
                                        image->set(i, j, k, v);
                                    }
                                }
                            }
                        }
                        break;
                    // -------------------------------------------------------
                    // invalid operation
                    default:
                        ASSERT(false, "Invalid operation code: " << *op);
                        // ignore error if build in Release configuration
                        break;
                };
                // next time this operation is encountered, use the next argument
                argidx[*op]++;
            }
        }
    } catch (exception& e) {
        cerr << "Caught exception: " << e.what() << endl;
        ok = false;
    }

    // -----------------------------------------------------------------------
    // write output image
    if (ok && !output_file.getValue().empty()) {
        ok = WriteImage(output_file.getValue().c_str(), image);
        if (!ok) {
            cerr << "Failed to write resulting image to file " << output_file.getValue() << "!" << endl;
        }
    }

    // -----------------------------------------------------------------------
    // clean up
    if (image) delete image;

    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
