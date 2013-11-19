/**
 * @file  ConvertDeformation.cxx
 * @brief Converts deformation field from one representation to another.
 *
 * Copyright (c) 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <iostream>
#include <string.h> // strncmp()
#include <memory>
#include "common/imageio.h"
#include "common/utilities.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dramms/basis.h>


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;


// ===========================================================================
// help
// ===========================================================================

void print_help()
{
    const string exec_name = exename();
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " [input options] [-F (ITK|DRAMMS|FSL)] [output options]" << endl;
    cout << endl;
    cout << "Description:" << endl;
    cout << "  This program can be used to convert a deformation field given in one representation" << endl;
    cout << "  to another. The known and supported representations, i.e., deformation field formats, are:" << endl;
    cout << endl;
    cout << "  - DRAMMS   The vector field components are stored consecutively, where" << endl;
    cout << "             the order of x and y appears to be flipped compared to the other" << endl;
    cout << "             deformation field formats, i.e., yxzyxzyxzyxz... Note that the x and y" << endl;
    cout << "             axis of the stored deformation field image itself are not flipped." << endl;
    cout << "             Each vector component is given in voxel units of the template space." << endl;
    cout << "             The size of each voxel is stored in the pixdim field of the NIfTI-1 header" << endl;
    cout << "             of the deformation field given in the NIfTI-1 image format." << endl;
    cout << "  - ITK      The vector field components are stored consecutively as xyzxyzxyz..." << endl;
    cout << "             Each vector component is given in physical units, i.e., the displacement" << endl;
    cout << "             is applied to the physical voxel coordinates in LPS orientation." << endl;
    cout << "  - FSL      The vector field components are stored in separate scalar images and" << endl;
    cout << "             these scalar images are concatentated to a single data file, i.e.," << endl;
    cout << "             xxx...yyy...zzz... Each vector component is given in physical units, i.e.," << endl;
    cout << "             is applied to the physical voxel coordinates in LPS orientation." << endl;
    cout << endl;
    cout << "  ATTENTION: If the format of the deformation field is DRAMMS, the arguments for" << endl;
    cout << "             x and y of any input option or operation have to be provided exchanged" << endl;
    cout << "             compared to the other formats. For example, if an image in the ITK format" << endl;
    cout << "             has size 120x100x12 voxels, the size is considered to be 100x120x12 in the" << endl;
    cout << "             DRAMMS format (even though the size of the image memory in each dimension" << endl;
    cout << "             is actually still the same)." << endl;
    cout << endl;
    cout << "Input options:" << endl;
    cout << "  -f (DRAMMS|ITK|FSL) Format of input deformation field. (default: DRAMMS)" << endl;
    cout << "  -d <nx>,<ny>[,<nz>] Size of input deformation field. If the input deformation" << endl;
    cout << "  or -d <file>        field is two-dimensional, either do not specify <nz> or" << endl;
    cout << "                      set it to 0 or 1. Alternatively, you can specify another" << endl;
    cout << "                      ANALYZE or NIfTI image file as argument. In this case, the" << endl;
    cout << "                      image size is extracted from the image header of this file." << endl;
    cout << "                      This option is only required if the deformation field is read" << endl;
    cout << "                      from one or more raw image data file(s), in which case the size" << endl;
    cout << "                      has to be specified before either one of the -i, -x, -y, or" << endl;
    cout << "                      -z options." << endl;
    cout << "  -p <dx>,<dy>[,<dz>] Voxel size in mm of input deformation field. If the input" << endl;
    cout << "  or -p <file>        deformation field is two-dimensional, either do not specify" << endl;
    cout << "                      <dz> or set it to any arbitray value. Alternatively, you can" << endl;
    cout << "                      specify another ANALYZE or NIfTI image file as argument. In" << endl;
    cout << "                      this case, the image size is extracted from the image header" << endl;
    cout << "                      of this file. This option is only required if the deformation" << endl;
    cout << "                      field is read from one or more raw image data file(s), in which" << endl;
    cout << "                      case the voxel size has to be specified before either one of the" << endl;
    cout << "                      -i, -x, -y, or -z options." << endl;
    cout << "  -i <file>           Image file of input deformation field." << endl;
    cout << "                      If the deformation field is stored in separate scalar images," << endl;
    cout << "                      use options -x, -y, and -z instead." << endl;
    cout << "  -x <file>           Image file storing x displacements only. If this option is" << endl;
    cout << "                      used to input the deformation field, also option -y is required" << endl;
    cout << "                      and option -i may not be used. If the deformation field is a" << endl;
    cout << "                      three-dimensional image, the displacements in the z dimension" << endl;
    cout << "                      have to be given as well using the -z option." << endl;
    cout << "  -y <file>           Image file storing y displacements only. If this option is" << endl;
    cout << "                      used to input the deformation field, also option -x is required" << endl;
    cout << "                      and option -i may not be used. If the deformation field is a" << endl;
    cout << "                      three-dimensional image, the displacements in the z dimension" << endl;
    cout << "                      have to be given as well using the -z option." << endl;
    cout << "  -z <file>           Image file storing z displacements only. If this option is" << endl;
    cout << "                      used to input the deformation field, also options -x  and -y are" << endl;
    cout << "                      required and option -i may not be used. If the deformation field is a" << endl;
    cout << "                      two-dimensional image, the value of this option must be '0'." << endl;
    cout << endl;
    cout << "Output options:" << endl;
    cout << "  -F (DRAMMS|ITK|FSL) Change image format to the specified format." << endl;
    cout << "  -3                  Change deformation field vectors to three-dimensional displacements" << endl;
    cout << "                      where the third component is constant zero in case of a two-dimensional" << endl;
    cout << "                      input deformation field." << endl;
    cout << "  -m <file>           Write magnitude of displacement vectors to the specified image file." << endl;
    cout << "  -o <file>           Write deformation field to the specified image file." << endl;
    cout << "  -x <file>           Write x displacements to the specified image file." << endl;
    cout << "  -y <file>           Write y displacements to the specified image file." << endl;
    cout << "  -z <file>           Write z displacements to the specified image file." << endl;
    cout << endl;
    cout << "Standard options:" << endl;
    cout << "  -v --verbose        Increase verbosity of output messages." << endl;
    cout << "  -h --help           Print help and exit." << endl;
    cout << "  --version           Print version information and exit." << endl;
    cout << endl << endl;
	cout << "Example:" << endl;
	cout << "  " << exec_name << " -f FSL -i def_fsl_format.hdr -F DRAMMS -o def_dramms_format.nii.gz" << endl;
	cout << "  " << exec_name << " -f DRAMMS -i def_dramms_format.hdr -F ITK -o def_itk_format.nii.gz" << endl;
    print_contact();
}

// ===========================================================================
// helpers
// ===========================================================================

// ---------------------------------------------------------------------------
const char* getarg(int& argc, char**& argv)
{
    const char* opt = *argv;
    argc--; argv++;
    if (*argv == NULL) {
        cerr << "Missing argument(s) for option " << opt << endl;
        cerr << "See help (-h option) for usage information." << endl;
        exit(1);
    }
    return *argv;
}

// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc,char *argv[])
{
    Ivector3d       dim;                               // input size
    Fvector3d       pixdim;                            // input voxel size
    Image::Format   fmt        = Image::FORMAT_DRAMMS; // format of input images
    Image*          deffield   = NULL;                 // current deformation field
    Image*          x_deffield = NULL;                 // input x displacements only
    Image*          y_deffield = NULL;                 // input y displacements only
    Image*          z_deffield = NULL;                 // input z displacements only
    int             verbose    = 0;                    // verbosity

    dim   .x = dim   .y = dim   .z = 0;
    pixdim.x = pixdim.y = pixdim.z = 0.0f;

    // Note: The auto_ptr instances take care of freeing the allocated
    //       memory even upon exit() due to an error during the execution.
    //       This way we do not need to take care of this any more explicitly.
    auto_ptr<Image> deffield_ptr(deffield);
    auto_ptr<Image> x_deffield_ptr(x_deffield);
    auto_ptr<Image> y_deffield_ptr(y_deffield);
    auto_ptr<Image> z_deffield_ptr(z_deffield);

    // -----------------------------------------------------------------------
    // help and version
    argc--; argv++;
    if (argc == 0) {
        print_help();
        exit(1);
    }
    for (int i = 0; i < argc; i++) {
        // -------------------------------------------------------------------
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose++;
        // -------------------------------------------------------------------
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help();
            exit(0);
        // -------------------------------------------------------------------
        } else if (strcmp(argv[i], "--version") == 0) {
            print_version("ConvertDeformation");
            exit(0);
        }
    }

    // -----------------------------------------------------------------------
    // read deformation field
    while (argc > 0) {
        string opt = argv[0];
        // -------------------------------------------------------------------
        if (opt == "-f") {
            const char* fmtstr = getarg(argc, argv);
            fmt = StringToImageFormat(fmtstr);
            if (fmt == Image::FORMAT_UNKNOWN) {
                cerr << "Invalid argument for option -f: " << fmtstr << endl;
                exit(1);
            }
        // -------------------------------------------------------------------
        } else if (opt == "-d") {
            const char* str = getarg(argc, argv);
            sscanf(str, "%d,%d,%d", &dim.x, &dim.y, &dim.z);
            if (dim.x <= 0 || dim.y <= 0 || dim.z < 0) {
                cerr << "Invalid argument for option -d: " << str << endl;
                exit(1);
            }
            if (dim.z == 0) dim.z = 1;
        // -------------------------------------------------------------------
        } else if (opt == "-p") {
            const char* str = getarg(argc, argv);
            sscanf(str, "%f,%f,%f", &pixdim.x, &pixdim.y, &pixdim.z);
            if (pixdim.x <= 0 || pixdim.y <= 0 || pixdim.z < 0) {
                cerr << "Invalid argument for option -p: " << str << endl;
                exit(1);
            }
        // -------------------------------------------------------------------
        } else if (opt == "-x") {
            if (x_deffield != NULL) {
                cerr << "Input image of x displacements already given!" << endl;
                cerr << "Please specify also remaining component images." << endl;
                exit(1);
            }
            const char* filename = getarg(argc, argv);
            if (IsMetaImage(filename) != -1) {
                x_deffield = ReadMetaDeformationField(filename, fmt);
            } else if (IsNiftiImage(filename) != -1) {
                x_deffield = ReadNiftiDeformationField(filename, fmt);
            } else {
                if (dim.x == 0 || dim.y == 0 || dim.z == 0 ||
                        pixdim.x == 0.0f || pixdim.y == 0.0f || pixdim.z == 0.0f) {
                    cerr << "Input file " << filename << " seems to be a raw image data file." << endl;
                    cerr << "Please specify dimension and voxel size using -d and -p options before option -x." << endl;
                    exit(1);
                }
                x_deffield = new Image(dim.x, dim.y, dim.z,
                                       pixdim.x, pixdim.y, pixdim.z,
                                       DT_FLOAT, 1, fmt);
                if (x_deffield == NULL) {
                    cerr << "Failed to allocate memory for x displacements!" << endl;
                    exit(1);
                }
                if (!ReadRawImage(filename, x_deffield)) {
                    delete x_deffield;
                    x_deffield = NULL;
                }
            }
            if (x_deffield == NULL) {
                cerr << "Failed to read x displacements from file " << filename << endl;
                exit(1);
            }
            if (y_deffield && (z_deffield || dim.z == 1)) break;
        // -------------------------------------------------------------------
        } else if (opt == "-y") {
            if (y_deffield) {
                cerr << "Input image of y displacements already given!" << endl;
                cerr << "Please specify also remaining component images." << endl;
                exit(1);
            }
            const char* filename = getarg(argc, argv);
            if (IsMetaImage(filename) != -1) {
                y_deffield = ReadMetaDeformationField(filename, fmt);
            } else if (IsNiftiImage(filename) != -1) {
                y_deffield = ReadNiftiDeformationField(filename, fmt);
            } else {
                if (dim.x == 0 || dim.y == 0 || dim.z == 0 ||
                        pixdim.x == 0.0f || pixdim.y == 0.0f || pixdim.z == 0.0f) {
                    cerr << "Input file " << filename << " seems to be a raw image data file." << endl;
                    cerr << "Please specify dimension and voxel size using -d and -p options before option -x." << endl;
                    exit(1);
                }
                y_deffield = new Image(dim.x, dim.y, dim.z,
                                       pixdim.x, pixdim.y, pixdim.z,
                                        DT_FLOAT, 1, fmt);
                if (y_deffield == NULL) {
                    cerr << "Failed to allocate memory for y displacements!" << endl;
                    exit(1);
                }
                if (!ReadRawImage(filename, y_deffield)) {
                    delete y_deffield;
                    y_deffield = NULL;
                }
            }
            if (y_deffield == NULL) {
                cerr << "Failed to read y displacements from file " << filename << endl;
                exit(1);
            }
            if (x_deffield && (z_deffield || dim.z == 1)) break;
        // -------------------------------------------------------------------
        } else if (opt == "-z") {
            if (z_deffield) {
                cerr << "Input image of z displacements already given!" << endl;
                cerr << "Please specify also remaining component images." << endl;
                exit(1);
            }
            const char* filename = getarg(argc, argv);
            if (strncmp(filename, "0", 2) == 0) break;
            if (IsMetaImage(filename) != -1) {
                z_deffield = ReadMetaDeformationField(filename, fmt);
            } else if (IsNiftiImage(filename) != -1) {
                z_deffield = ReadNiftiDeformationField(filename, fmt);
            } else {
                if (dim.x == 0 || dim.y == 0 || dim.z == 0 ||
                        pixdim.x == 0.0f || pixdim.y == 0.0f || pixdim.z == 0.0f) {
                    cerr << "Input file " << filename << " seems to be a raw image data file." << endl;
                    cerr << "Please specify dimension and voxel size using -d and -p options before option -x." << endl;
                    exit(1);
                }
                z_deffield = new Image(dim.x, dim.y, dim.z,
                                       pixdim.x, pixdim.y, pixdim.z,
                                       DT_FLOAT, 1, fmt);
                if (z_deffield == NULL) {
                    cerr << "Failed to allocate memory for z displacements!" << endl;
                    exit(1);
                }
                if (!ReadRawImage(filename, z_deffield)) {
                    delete z_deffield;
                    z_deffield = NULL;
                }
            }
            if (z_deffield == NULL) {
                cerr << "Failed to read z displacements from file " << filename << endl;
                exit(1);
            }
            if (z_deffield->region.nz == 1) {
                cerr << "Number of voxels in z dimension set to less or equal to 1." << endl;
                cerr << "Displacement vectors in two dimensions do not have a z component!" << endl;
                exit(1);
            }
            if (x_deffield && y_deffield) break;
        // -------------------------------------------------------------------
        } else if (opt == "-i") {
            if (x_deffield || y_deffield || z_deffield) {
                cerr << "Either specify input deformation field using option -i";
                cerr << " or separate input component images using options -x, -y,";
                cerr << " and -z (in case of three dimensions).";
                cerr << endl;
                exit(1);
            }
            const char* filename = getarg(argc, argv);
            if (IsMetaImage(filename) != -1) {
                deffield = ReadMetaDeformationField(filename, fmt);
            } else if (IsNiftiImage(filename) != -1) {
                deffield = ReadNiftiDeformationField(filename, fmt);
            } else {
                if (dim.x == 0 || dim.y == 0 || dim.z == 0 ||
                        pixdim.x == 0.0f || pixdim.y == 0.0f || pixdim.z == 0.0f) {
                    cerr << "Input file " << filename << " seems to be a raw image data file." << endl;
                    cerr << "Please specify dimension and voxel size using -d and -p options before option -i." << endl;
                    exit(1);
                }
                deffield = new Image(dim.x, dim.y, dim.z,
                                     pixdim.x, pixdim.y, pixdim.z,
                                     DT_FLOAT, dim.z > 1 ? 3 : 2, fmt);
                if (deffield == NULL) {
                    cerr << "Failed to allocate memory for input deformation field!" << endl;
                    exit(1);
                }
                if (!ReadRawImage(filename, deffield)) {
                    delete deffield;
                    deffield = NULL;
                }
            }
            if (deffield == NULL) {
                cerr << "Failed to read deformation field from file " << filename << endl;
                exit(1);
            }
            break;
        // -------------------------------------------------------------------
        } else if (opt != "-v") {
            cerr << "Unknown option or missing deformation field input option before " << opt << endl;
            cerr << "See help (-h option) for usage information." << endl;
            exit(1);
        }
        argc--; argv++;
    }

    if (argc == 0 || (deffield == NULL && x_deffield == NULL)) {
        cerr << "No input deformation field specified or missing component images!" << endl;
        cerr << "See help (-h option) for usage information." << endl;
        exit(1);
    }
    argc--; argv++;

    // -----------------------------------------------------------------------
    // assemble deformation field from component images
    if (deffield == NULL) {
        // ensure (at least) that all component images have the same size
        if (y_deffield->region.nx != x_deffield->region.nx ||
                y_deffield->region.ny != x_deffield->region.ny ||
                y_deffield->region.nz != x_deffield->region.nz) {
            fprintf(stderr, "The image dimensions of the component images do not match!\n");
            exit(1);
        }
        if (z_deffield) {
            if (z_deffield->region.nx != x_deffield->region.nx ||
                    z_deffield->region.ny != x_deffield->region.ny ||
                    z_deffield->region.nz != x_deffield->region.nz) {
                fprintf(stderr, "The image dimensions of the component images do not match!\n");
                exit(1);
            }
        }
        // allocate deformation field
        deffield = new Image(x_deffield->region.nx, x_deffield->region.ny, x_deffield->region.nz,
                             DT_FLOAT, (z_deffield ? 3 : 2), x_deffield->imgfmt);
        if (deffield == NULL) {
            fprintf(stderr, "Failed to allocate memory for combined deformation field!\n");
            exit(1);
        }
        deffield->CopyRegion     (x_deffield);
        deffield->CopyTransform  (x_deffield);
        deffield->CopyDataScaling(x_deffield);
        deffield->CopyMetaData   (x_deffield);
        // copy displacement components to deformation field
        for (int k = 0; k < deffield->region.nz; k++) {
            for (int j = 0; j < deffield->region.ny; j++) {
                for (int i = 0; i < deffield->region.nx; i++) {
                    deffield->set(i, j, k, 0, x_deffield->get(i, j, k));
                    deffield->set(i, j, k, 1, y_deffield->get(i, j, k));
                }
            }
        }
        if (z_deffield) {
            for (int k = 0; k < deffield->region.nz; k++) {
                for (int j = 0; j < deffield->region.ny; j++) {
                    for (int i = 0; i < deffield->region.nx; i++) {
                        deffield->set(i, j, k, 2, z_deffield->get(i, j, k));
                    }
                }
            }
            delete z_deffield; z_deffield = NULL;
        }
        delete x_deffield; x_deffield = NULL;
        delete y_deffield; y_deffield = NULL;
    }

    // reset dim and pixdim as used by following options
    dim.x = dim.y = dim.z = 0;
    pixdim.x = pixdim.y = pixdim.z = 0.0f;

    // -----------------------------------------------------------------------
    // perform operations
    while (argc > 0) {
        string opt = argv[0];
        // -------------------------------------------------------------------
        if (opt == "-3") {
            if (deffield->GetNumberOfComponents() == 2) {
                Image* outdef = ConvertTo3DDeformationField(deffield);
                if (outdef == NULL) {
                    cerr << "Failed to allocate memory for three-dimensional deformation field!" << endl;
                    break;
                }
                delete deffield;
                deffield = outdef;
            }
        // -------------------------------------------------------------------
        } else if (opt == "-F") {
            const char* fmtstr = getarg(argc, argv);
            Image::Format fmt = StringToImageFormat(fmtstr); // local variable!
            if (fmt == Image::FORMAT_UNKNOWN) {
                cerr << "Invalid argument for option -F: " << fmtstr << endl;
                break;
            }
            deffield->SetFormat(fmt);
        // -------------------------------------------------------------------
        } else if (opt == "-x") {
            const char* filename = getarg(argc, argv);
            x_deffield = new Image(deffield->region.nx,
                                   deffield->region.ny,
                                   deffield->region.nz,
                                   DT_FLOAT, 1,
                                   deffield->imgfmt);
            if (x_deffield == NULL) {
                cerr << "Failed to allocate memory for image of x displacements!" << endl;
                break;
            }
            x_deffield->CopyTransform(deffield);
            x_deffield->CopyDataScaling(deffield);
            x_deffield->CopyMetaData(deffield);
            x_deffield->hdr.intent_code = NIFTI_INTENT_NONE; // no longer a vector field
            vector<float> v;
            for (int k = 0; k < deffield->region.nz; k++) {
                for (int j = 0; j < deffield->region.ny; j++) {
                    for (int i = 0; i < deffield->region.nx; i++) {
                        x_deffield->set(i, j, k, deffield->get(i, j, k, 0));
                    }
                }
            }
            if (!WriteImage(filename, x_deffield)) {
                cerr << "Failed to write x displacements to file " << filename << "!" << endl;
                break;
            }
            delete x_deffield;
            x_deffield = NULL;
        // -------------------------------------------------------------------
        } else if (opt == "-y") {
            const char* filename = getarg(argc, argv);
            y_deffield = new Image(deffield->region.nx,
                                   deffield->region.ny,
                                   deffield->region.nz,
                                   DT_FLOAT, 1,
                                   deffield->imgfmt);
            if (y_deffield == NULL) {
                cerr << "Failed to allocate memory for image of y displacements!" << endl;
                break;
            }
            y_deffield->CopyTransform(deffield);
            y_deffield->CopyDataScaling(deffield);
            y_deffield->CopyMetaData(deffield);
            y_deffield->hdr.intent_code = NIFTI_INTENT_NONE; // no longer a vector field
            vector<float> v;
            for (int k = 0; k < deffield->region.nz; k++) {
                for (int j = 0; j < deffield->region.ny; j++) {
                    for (int i = 0; i < deffield->region.nx; i++) {
                        y_deffield->set(i, j, k, deffield->get(i, j, k, 1));
                    }
                }
            }
            if (!WriteImage(filename, y_deffield)) {
                cerr << "Failed to write y displacements to file " << filename << "!" << endl;
                break;
            }
            delete y_deffield;
            y_deffield = NULL;
        // -------------------------------------------------------------------
        } else if (opt == "-z") {
            const char* filename = getarg(argc, argv);
            if (deffield->region.nz <= 1) {
                cerr << "Deformation field is two-dimensional! Cannot write z displacements to file " << filename << "." << endl;
                break;
            }
            z_deffield = new Image(deffield->region.nx,
                                   deffield->region.ny,
                                   deffield->region.nz,
                                   DT_FLOAT, 1,
                                   deffield->imgfmt);
            if (z_deffield == NULL) {
                cerr << "Failed to allocate memory for image of z displacements!" << endl;
                break;
            }
            z_deffield->CopyTransform(deffield);
            z_deffield->CopyDataScaling(deffield);
            z_deffield->CopyMetaData(deffield);
            z_deffield->hdr.intent_code = NIFTI_INTENT_NONE; // no longer a vector field
            vector<float> v;
            for (int k = 0; k < deffield->region.nz; k++) {
                for (int j = 0; j < deffield->region.ny; j++) {
                    for (int i = 0; i < deffield->region.nx; i++) {
                        z_deffield->set(i, j, k, deffield->get(i, j, k, 2));
                    }
                }
            }
            if (!WriteImage(filename, z_deffield)) {
                cerr << "Failed to write z displacements to file " << filename << "!" << endl;
                break;
            }
            delete z_deffield;
            z_deffield = NULL;
        // -------------------------------------------------------------------
        } else if (opt == "-o") {
            const char* filename = getarg(argc, argv);
            if (!WriteImage(filename, deffield)) {
                cerr << "Failed to write deformation field to file " << filename << "!" << endl;
                break;
            }
        // -------------------------------------------------------------------
        } else if (opt == "-m") {
            const char* filename = getarg(argc, argv);
            Image* mag = new Image(deffield->region.nx,
                                   deffield->region.ny,
                                   deffield->region.nz,
                                   DT_FLOAT, 1,
                                   deffield->imgfmt);
            if (mag == NULL) {
                cerr << "Failed to allocate memory for magnitude image!" << endl;
                break;
            }
            mag->CopyTransform(deffield);
            mag->CopyMetaData(deffield);
            mag->hdr.intent_code = NIFTI_INTENT_NONE; // no longer a vector field
            float m;
            vector<float> v;
            for (int k = 0; k < deffield->region.nz; k++) {
                for (int j = 0; j < deffield->region.ny; j++) {
                    for (int i = 0; i < deffield->region.nx; i++) {
                        deffield->get(i, j, k, v);
                        m = 0.0f;
                        for (size_t n = 0; n < v.size(); n++) m += v[n] * v[n];
                        m = sqrt(m);
                        mag->set(i, j, k, m);
                    }
                }
            }
            if (!WriteImage(filename, mag)) {
                cerr << "Failed to write magnitude image to file " << filename << "!" << endl;
                delete mag;
                break;
            }
            delete mag;
        // -------------------------------------------------------------------
        } else if (opt == "-f") {
            cerr << "Deformation field has been loaded already!" << endl;
            cerr << "Input option -f no longer valid. Did you mean -F instead?" << endl;
            cerr << "See help (-h option) for usage information." << endl;
        } else if (opt == "-i") {
            cerr << "Deformation field has been loaded already! Option -i cannot be" << endl;
            cerr << "given more than once or after the input deformation field has" << endl;
            cerr << "been specified using the -x, -y, and -z options." << endl;
            cerr << "See help (-h option) for usage information." << endl;
        } else if (opt == "-d" || opt == "-p") {
            cerr << "Option " << opt << " can only be used before option -i or the" << endl;
            cerr << "options -x, -y, and -z, but not after the deformation field" << endl;
            cerr << "has been loaded already. See help (-h option) for usage information." << endl;
            break;
        // -------------------------------------------------------------------
        } else if (opt != "-v") {
            cerr << "Unknown option: " << opt << endl;
            cerr << "See help (-h option) for a list of available options." << endl;
            break;
        }
        // -------------------------------------------------------------------
        argc--; argv++;
    }

    // -----------------------------------------------------------------------
    // success, if all arguments were processed
    exit(argc);
}
