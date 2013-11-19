/**
 * @file  ApplyTransform.cxx
 * @brief Apply affine or deformable transformation to scalar image.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <common/utilities.h> // ApplyTransform()
#include <common/imageio.h>   // ReadImage(), WriteImage()

#include <dramms/basis.h>  // exename(), print_contact()


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;


// ===========================================================================
// help
// ===========================================================================

// ---------------------------------------------------------------------------
void print_help()
{
    string exec_name = exename();
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_image> <transform_file> <output_image>" << endl;
    cout << endl;
    cout << "Description:" << endl;
    cout << "  This program applies a given transformation on a scalar image. The transformation" << endl;
    cout << "  can be a deformation field or an affine transformation matrix. Two interpolation" << endl;
    cout << "  options are supported, -- trilinear (default) and nearest neighborhood." << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
    cout << "  <input_image>        Input image file." << endl;
    cout << "  <transform_file>     Transformation file, i.e., either an affine transformation matrix" << endl;
    cout << "                       in ASCII format or a deformation field in DRAMMS format." << endl;
    cout << "  <output_image>       Output image file." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
    cout << "  -t <template_file>   Template image. Only required if an affine transformation is applied." << endl;
    cout << "                       In case of a deformation field, the header of this template image is" << endl;
    cout << "                       used instead of the deformation field header if specified." << endl;
    cout << "  -n                   Use nearest neighbor interpolation." << endl;
    cout << "                       (default: linear interpolation)" << endl;
    cout << "  -v                   Increase verbosity of output messages." << endl;
    cout << "  -h                   Print help and exit." << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << exec_name << " subj.nii.gz def.nii.gz warpedsubj.nii.gz" << endl;
    cout << endl;
    cout << "  " << exec_name << " subj.nii.gz A2B_affine.mat A2B_affine.nii.gz" << endl;
    cout << endl;
    print_contact();
}

// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main (int argc,char *argv[])
{
    bool ok = true;

    // default options
    const char* template_file = NULL;
    bool        interpolate   = true;
    int         verbose       = 0;

    // show usage if no arguments were given
    if (argc == 1) {
        print_help();
        exit(1);
    }

    // parse arguments
    int c = -1;
    while ((c = getopt(argc, argv, "nvht:")) != -1) {
        switch (c) {
            case 'n':
                interpolate = false; // use nearest neighbor instead
                break;

            case 't':
                template_file = optarg;
                break;

            case 'v':
                verbose++;
                break;

            case 'h':
                print_help();
                exit(0);

            default:
                // error message printed by getopt()
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 3) {
        cerr << "Missing required arguments!" << endl;
        cerr << "See help (-h option) for a list of required arguments." << endl;
        exit(1);
    }

    const char* input_file     = argv[0];
    const char* transform_file = argv[1];
    const char* output_file    = argv[2];

    // -----------------------------------------------------------------------
    // read input image(s)
    Sequence input_images = ReadSequence(input_file);
    if (input_images.empty()) {
        cerr << "Failed to read image from file " << input_file << endl;
        ok = false;
    }
    Image* input_image = input_images[0];

    // -----------------------------------------------------------------------
    // read header of reference image
    Image::Header* template_hdr = NULL;
    if (template_file) {
        template_hdr = nifti_read_header(template_file, NULL, 1);
        if (template_hdr == NULL) {
            cerr << "Failed to read header of template image " << template_file << "!" << endl;
            cerr << "The file does either not exist or is not a NIfTI-1 or ANALYZE 7.5 image file." << endl;
            exit(EXIT_FAILURE);
        }
        if (template_hdr->dim[3] < 1) {
            template_hdr->dim   [3] = 1;
            template_hdr->pixdim[3] = 1.0f;
        }
        if (input_image->imgfmt == Image::FORMAT_DRAMMS) {
            *template_hdr = PermuteXY(*template_hdr);
        }
    }

    // -----------------------------------------------------------------------
    // read transformation
    Image*           D = NULL;
    Image::Transform T;

    if (ok) {
        if (!ReadTransform(transform_file, D, T, input_image->imgfmt)) {
            cerr << "Failed to read transformation from file " << transform_file << endl;
            ok = false;
        } else if (D && template_hdr) {
            if (template_hdr->dim[1] != D->hdr.dim[1] ||
                    template_hdr->dim[2] != D->hdr.dim[2] ||
                    template_hdr->dim[3] != D->hdr.dim[3]) {
                cerr << "Template image size does not match size of deformation field!" << endl;
                ok = false;
            } else {
                Image template_image;
                template_image.hdr = *template_hdr;
                D->CopyTransform(&template_image);
            }
        }
    }

    // -----------------------------------------------------------------------
    // apply transformation
    Sequence output_images;
    for (size_t i = 0; i < input_images.size() && ok; i++) {
        Image *output_image = NULL;
		input_image = input_images[i];
        if (D) {
            output_image = ApplyTransform(input_image, D, interpolate);
        } else if (template_hdr == NULL) {
            cerr << "Missing reference/template image file!" << endl;
            ok = false;
        } else {
            Image reference;
            reference.hdr       = *template_hdr;
            reference.imgfmt    = input_image->imgfmt;
            reference.region.nx = template_hdr->dim[1];
            reference.region.ny = template_hdr->dim[2];
            reference.region.nz = template_hdr->dim[3];
            reference.UpdateTransforms();
            output_image = ApplyTransform(input_image, T, &reference, interpolate);
        }
        if (ok && output_image == NULL) {
            cerr << "Failed to allocate memory for output image!" << endl;
            ok = false;
        }
        output_images.push_back(output_image);
    }

    // -----------------------------------------------------------------------
    // output intensity range before and after deformation
    if (ok && verbose > 0) {
        float min, max;
        Image *output_image;

		for (size_t i = 0; i < input_images.size() && ok; i++) {
			input_image  = input_images[i];
			output_image = output_images[i];
			cout << "In channel " << i+1 << " out of " << input_images.size() << ":" << endl;
			cout << "Image size before warping: [" << input_image->region.ny << ", "
                << input_image->region.nx << ", " << input_image->region.nz << "]" << endl;
			cout << "Image size after warping:  [" << output_image->region.ny << ", "
                << output_image->region.nx << ", " << output_image->region.nz << "]" << endl;
			GetIntensityRange(input_image, min, max);
			cout << "Intensity range before warping: [" << min << ", " << max << "]" << endl;
			GetIntensityRange(output_image, min, max);
			cout << "Intensity range after warping:  [" << min << ", " << max << "]" << endl << endl;
		}
    }

    // -----------------------------------------------------------------------
    // save result
    if (ok) {
        if (WriteSequence(output_file, output_images)) {
            if (verbose > 0) {
                cout << "Wrote transformed image to file " << output_file << endl;
            }
        } else {
            cerr << "Failed to write transformed image to file " << output_file << endl;
            ok = false;
        }
    }

    // -----------------------------------------------------------------------
    // clean up
    input_images.Clear();
    output_images.Clear();
    if (template_hdr) delete template_hdr;
    if (D)            delete D;

    exit(ok ? 0 : 1);
}
