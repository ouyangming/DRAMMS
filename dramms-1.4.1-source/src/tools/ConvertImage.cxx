/**
 * @file  ConvertImage.cxx
 * @brief Convert image from one datatype to another.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <stdio.h>
#include <nifti1_io.h>
#include <unistd.h>
#include <stdlib.h>
#include <common/imageio.h>
#include <common/utilities.h>

#include <dramms/basis.h>


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;


// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // -----------------------------------------------------------------------
    // define command-line arguments
    StringArg type("t", "type",
            "Output datatype. Valid datatype identifiers are 'byte' or 'uchar',"
            " 'short', 'ushort', 'int', 'uint', and 'float'. If not specified,"
            " the datatype of the input image is used.",
            false, "DT_UNKNOWN", "<datatype>");

    StringArg orient("o", "orient",
            "Change image header such that it specifies the given orientation."
            " Note that this operation does not change the orientation of the image data!"
            " To change the orientation of both image data and header, use"
            " --orient-data instead. The desired orientation has to be specified as"
            " two or three letter orientation code, respectively, such as RAI or LAS,"
            " for example. Note that each letter denotes the location of the *first*"
            " voxel on the respective axis. This is in compliance with the ITK and other"
            " tools for neurological imaging, whereas DICOM would specify the location of"
            " the *last* voxel instead. Hence, to orient the image to the default DICOM"
            " orientation LPS, specify RAI as argument for this operation.",
            false, "", "<ox><oy>[<oz>]");

    StringArg orient_data("", "orient-data",
            "Orient image such that the image data is stored in the given orientation."
            " The desired orientation has to be specified as two or three letter"
            " orientation code, respectively, such as RAI or LAS, for example."
            " Changing the orientation of the data also updates the orientation"
            " matrices of the header such that the images are still shown in the"
            " same orientation. To correct the orientation of an image, use (also) the"
            " -o (--orient) option with the correct orientation as argument."
            " Note that each letter denotes the location of the *first* voxel on the"
            " respective axis. This is in compliance with the ITK and other tools for"
            " neurological imaging, whereas DICOM would specify the location of the"
            " *last* voxel instead. Hence, to orient the image data to the default"
            " DICOM orientation LPS, specify RAI as argument for this operation.",
            false, "", "<ox><oy>[<oz>]");

    StringArg transform_reference_file("", "copy-transform",
            "Copy transformation (qform and/or sform) from given reference image.",
            false, "", "<file>");

    FloatArg min("m", "min",
            "Minimum output intensity. If the given value is outside the range of"
            " intensities that can be represented by the output datatype of the image,"
            " it is set to the closest limit of this range. See help of --scale for"
            " details on how the maximum output intensity affects the rescaling of"
            " intensities.",
            false, 0, "<float>");

    FloatArg max("M", "max",
            "Maximum output intensity. If the given value is outside the range of"
            " intensities that can be represented by the output datatype of the image,"
            " it is set to the closest limit of this range. See help of --scale for"
            " details on how the maximum output intensity affects the rescaling of"
            " intensities.",
            false, 0, "<float>");

    SwitchArg scale("S", "scale",
            "Request the application of the scaling function of the NIfTI-1 header"
            " including a rescaling of intensities to the range [min, max] specified"
            " by the --min and --max arguments.  If the specified minimum output intensity"
            " is greater or equal to the desired maximum output intensity, the default"
            " output range is used, i.e, the intensity range of the input image is preserved"
            " if the datatype does not change. Otherwise, if the desired output datatype"
            " differs from the input datatype, the intensities are rescaled to the maximum"
            " range of the output datatype. Moreover, if the given intensity range is signed,"
            " i.e., min < 0, the intensities are rescaled such that their sign remains"
            " unchanged. Otherwise, a shift (intercept) is applied in order to avoid"
            " the truncation of negative intensities. Note that the header of the output"
            " image by default will encode the applied scaling function such that the"
            " intensities displayed in a viewer remain unchanged by the scaling except for"
            " a possible typical discretization error. See also the help of the --reset"
            " option.");

    SwitchArg smooth("s", "smooth",
            "Smooth image before scaling the intensities.");

    SwitchArg eff("e", "effective",
            "Enable the rescaling of the intensities within the \"effective\" range"
            " of intensities. Intensities outside the effective range are truncated."
            " The effective range is determined by the histogram of intensities.");

    SwitchArg reset_scaling("", "reset-scaling",
            "Request that the scl_slope and scl_inter entries of the NIfTI-1 header"
            " of the output image be reset to 0. Otherwise, these entries will encode"
            " the inverse scaling function that has been applied such that intensities"
            " after the conversion scaled by this function correspond to the intensities"
            " of the input image scaled by the scaling function of this image.");

    PositionalArg input_file("input_file",
            "Input image file. The input image can either be a scalar intensity image"
            " or a vector field such as in particular a deformation field.",
            true, "", "<input_file>");

    PositionalArg output_file("output_file",
            "Output image file.",
            true, "", "<output_file>");

    MultiSwitchArg verbose("v", "verbose", "Increase verbosity of output messages.");

    SwitchArg analyze("a", "analyze",
            "Force output image to be in ANALYZE 7.5 format instead of NIfTI-1."
            " This option is only required if either no file name extension is"
            " specified or if the given extension is .hdr or .img.");

    SwitchArg nifti("n", "nifti",
            "Force output image to be in NIfTI-1 format instead of ANALYZE 7.5."
            " This option is only required if the given extension is .hdr or .img."
            " and the input file is in ANALYZE format. Otherwise, the output"
            " format defaults to NIfTI-1 if no Meta image extension is given"
            " (.mhd, .raw).");

    // -----------------------------------------------------------------------
    // parse command-line
    try {
        CmdLine cmd(
                // program identification
                "ConvertImage", PROJECT,
                // program description
                "This program converts a NIfTI image into another NIfTI image of the desired"
                " orientation, intensity range, and/or output datatype."
                "\n"
                " The scl_slope and scl_inter values of the NIfTI-1 header of the output image are"
                " recalculated to reflect any scaling or rescaling that has been applied to the intensity"
                " values on disk. Hence, by applying this scaling function, one can recover the values"
                " of the input image which correspond to the intensities which result from scaling the"
                " input values according to the scaling function of the NIfTI-1 header of the input image.",
                // example usage
                "",
                // version information
                RELEASE, "2012 University of Pennsylvania");

        // the constructor of the CmdLine class has already added the standard
        // options for help output and version information

        cmd.add(type);
        cmd.add(orient);
        cmd.add(orient_data);
        cmd.add(transform_reference_file);
        cmd.add(min);
        cmd.add(max);
        cmd.add(scale);
        cmd.add(smooth);
        cmd.add(eff);
        cmd.add(reset_scaling);
        cmd.add(analyze);
        cmd.add(nifti);
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

    if      (verbose.getValue() > 1) g_debug = DEBUG_ALL;
    else if (verbose.getValue() > 0) g_debug = DEBUG_WARN;
    else                             g_debug = DEBUG_ERROR;

    if (analyze.getValue() && nifti.getValue()) {
        cerr << "Request the output image to be either in ANALYZE 7.5 or NIfTI-1 format, but not both!" << endl;
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

    int x_orient = 0, y_orient = 0, z_orient = 0;
    if (!orient.getValue().empty() &&
            !StringToOrientation(orient.getValue().c_str(), x_orient, y_orient, z_orient)) {
        cerr << "Invalid image orientation code: " << orient.getValue() << endl;
        exit(EXIT_FAILURE);
    }
    int x_orient_data = 0, y_orient_data = 0, z_orient_data = 0;
    if (!orient_data.getValue().empty() &&
            !StringToOrientation(orient_data.getValue().c_str(), x_orient_data, y_orient_data, z_orient_data)) {
        cerr << "Invalid image orientation code: " << orient_data.getValue() << endl;
        exit(EXIT_FAILURE);
    }

    // -----------------------------------------------------------------------
    // read image - in ITK format such that user input (such as orientation)
    //              are interpreted as expected and not with x and y flipped
    Image* image = ReadImage(input_file.getValue().c_str(), Image::FORMAT_ITK);
    if (image == NULL) {
        cerr << "Failed to read image from file " << input_file.getValue() << '!' << endl;
        exit(EXIT_FAILURE);
    }

    // -----------------------------------------------------------------------
    // cast/scale image intensities
    Image* tmp = NULL;
    if (scale.getValue() || eff.getValue()) {
        if (verbose.getValue() > 0) {
            cout << "Casting and rescaling image data..." << endl;
        }
        short datatype = dt == DT_UNKNOWN ? image->hdr.datatype : dt;
        float lmin, lmax;
        GetDatatypeRange(datatype, lmin, lmax);
        float omin, omax;
        if      (min.getValue() < lmin) omin = lmin;
        else if (min.getValue() > lmax) omin = lmax;
        else                            omin = min.getValue();
        if      (max.getValue() < lmin) omax = lmin;
        else if (max.getValue() > lmax) omax = lmax;
        else                            omax = max.getValue();
        tmp = CastImage(image, datatype, omin, omax, smooth.getValue(), eff.getValue());
    } else if (dt != DT_UNKNOWN && dt != image->hdr.datatype) {
        if (verbose.getValue() > 0) {
            cout << "Casting image data without rescaling (i.e., truncation may occur!)..." << endl;
        }
        tmp = CastImage(image, dt, false);
    } else {
        tmp = image;
    }
    if (tmp != image) {
        delete image;
        image = tmp;
    }
    if (image == NULL) {
        cerr << "Failed to cast image to output datatype!" << endl;
        exit(EXIT_FAILURE);
    }

    // -----------------------------------------------------------------------
    // copy transform from reference image
    if (!transform_reference_file.getValue().empty()) {
        if (verbose.getValue() > 0) {
            cout << "Copying q-form and s-form transformations from " << transform_reference_file.getValue() << "..." << endl;
        }
        Image tmp;
        if (IsNiftiImage(transform_reference_file.getValue().c_str(), &tmp.hdr)) {
            image->CopyTransform(&tmp);
        } else if (IsMetaImage(transform_reference_file.getValue().c_str(), &tmp.hdr)) {
            image->CopyTransform(&tmp);
        } else {
            cerr << "Failed to read header of image file " << transform_reference_file.getValue() << '!' << endl;
            delete image;
            exit(EXIT_FAILURE);
        }
    }

    // -----------------------------------------------------------------------
    // orient image data
    if (x_orient_data != 0) {
        if (verbose.getValue() > 0) {
            cout << "Changing orientation of image data to " << OrientationToString(x_orient_data, y_orient_data, z_orient_data) << "..." << endl;
        }
        if (!OrientImage(image, x_orient_data, y_orient_data, z_orient_data)) {
            cerr << "Failed to change orientation of image data to " << OrientationToString(x_orient_data, y_orient_data, z_orient_data) << "!" << endl;
            delete image;
            exit(EXIT_FAILURE);
        }
    }

    // -----------------------------------------------------------------------
    // orient image header
    if (x_orient != 0) {
        if (verbose.getValue() > 0) {
            cout << "Changing orientation of image to " << orient.getValue() << "..." << endl;
        }
        image->SetOrientation(x_orient, y_orient, z_orient);
    }

    // -----------------------------------------------------------------------
    // reset scaling function
    if (reset_scaling.getValue()) {
        if (verbose.getValue() > 0) {
            cout << "Setting scl_slope and scl_inter of header to zero..." << endl;
        }
        image->hdr.scl_slope = 0.0f;
        image->hdr.scl_inter = 0.0f;
    }
	image->hdr.cal_min = 0.0f;
	image->hdr.cal_max = 0.0f;

    // -----------------------------------------------------------------------
    // write image
    if (verbose.getValue() > 0) {
        cout << "Writing converted image to file " << output_file.getValue() << "..." << endl;
    }
    if (!WriteImage(output_file.getValue().c_str(), image,
                    nifti.getValue() ? Image::FILE_FORMAT_NIFTI_2
                                     : (analyze.getValue() ? Image::FILE_FORMAT_ANALYZE
                                                           : Image::FILE_FORMAT_UNKNOWN))) {
        cerr << "Failed to write image to file " << output_file.getValue() << '!' << endl;
        delete image;
        exit(EXIT_FAILURE);
    }

    // clean up
    delete image;

    exit(EXIT_SUCCESS);
}
