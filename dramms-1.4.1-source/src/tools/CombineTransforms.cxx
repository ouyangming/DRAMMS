/**
 * @file  CombineTransforms.cxx
 * @brief This program combines two transformations.
 *
 * This program requires two input files and writes one output file.
 * Both input files represent a transformation, either a deformable
 * transformation represented by a displacement vector field or an affine
 * transformation represented by a 4x4 matrix.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <string>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "common/imageio.h"
#include "common/utilities.h"

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
    const Image::Format fmt = Image::FORMAT_DRAMMS; // image format used
    bool                ok  = true;                 // program execution so far successful

    // -----------------------------------------------------------------------
    // define command-line arguments
    SwitchArg concatenate("c", "concatenate",
            "Concatenate two input transformations such that T3(x) = (T2 o T1)(x) = T2[T1(x)]."
            " The input transformations T1 and T2 can be either affine, deformable, or a mixture of both."
            "\n"
            "For example:"
            "\n"
            "1. If T1 is a deformation field warping image A to image B, and T2 is a deformation" 
            " from image B to image C, then T3 will be a deformation field warping image A to image C."
            "\n"
            "2. If T1 is an affine 4x4 transformation matrix (where the fourth row is [0 0 0 1])"
            " from the subject image A to the intermediate affine registered image A2B_affine,"
            " and T2 is a deformation from A2B_affine to the template image B,"
            " then T3 will be the unified deformation which warps the subject image A to the"
            " template image B. This composition is in particular applied by dramms."
            "\n"
            " Note that in the latter case, when composing an affine transformation and a"
            " deformation field, a reference image defining the image space before the affine"
            " transformation (i.e. image A in the example) must be specified using the -f (--affine_from)"
            " option. Moreover, if the image space A2B_affine differs from the image space B,"
            " option -t (--affine_to) has to be used to specify a reference image for the"
            " A2B_affine image space.");

    SwitchArg add("a", "add",
            "Add displacement vectors of deformation fields, i.e., T3(x) = T1(x) + T2(x)."
            " Requires two input transformations having identical dimension.");

    SwitchArg subtract("s", "subtract",
            "Subtract affine or deformable transformation from image A to image B from deformation"
            " field which deforms image A to image C, i.e., T3(x) = T1(x) - T2(x)."
            " If a deformable transformation is subtracted, both deformation fields must have"
            " identical dimension.");

    SwitchArg mean("m", "mean",
            "Compute mean displacement vector field, i.e, T3(x) = (T1(x) + T2(x))/2."
            " Requires two input deformation fields of the same size.");

    StringArg affinefrom_file("f", "affine_from",
            "Reference image file that defines the image space the affine transformation is from."
            " Used only if an affine transformation and a deformation field are being composed."
            " See option -c (--concatenate) for more details.",
            false, "", "<file>");

	StringArg affineto_file("t", "affine_to",
            "Reference image file that defines the image space the affine transformation is to."
            " Used only if an affine transformation and a deformation field are being composed."
            " See option -c (--concatenate) for more details.",
            false, "", "<file>");
			
    MultiSwitchArg verbose("v", "verbose",
            "Increase verbosity of output messages.");

    PositionalArg t1_file    ("T1", "First input transformation file.",  true, "", "<T1>");
    PositionalArg t2_file    ("T2", "Second input transformation file.", true, "", "<T2>");
    PositionalArg output_file("T3", "Output transformation file.",       true, "", "<T3>");

    // -----------------------------------------------------------------------
    // parse command-line
    try {

        CmdLine cmd(
                // program identification
                "CombineTransforms", PROJECT,
                // program description
                "This program combines two input transformations. In case of the concatenation of the"
                " transformations (the default operation), the input transformations T1 and T2 can"
                " be either affine, deformable, or a mixture of both. Otherwise, the input transformations"
                " must be deformation fields.",
                // example usage 1
                "EXENAME -c -f A.nii.gz -t B.nii.gz A2B_affine.mat def_A2Baffine_to_B.nii.gz def_A2B.nii.gz\n"
                "Computes the unified deformation between A and B given affine and deformable components.",
				// example usage 2
                "EXENAME -c -f A.nii.gz -t B.nii.gz A2B_affine.mat def_B2C.nii.gz def_A2C.nii.gz\n"
                "Computes the unified deformation between A and C given affine A->B and deformable B->C.",
                // version information
                RELEASE, "2012 University of Pennsylvania");

        // the constructor of the CmdLine class has already added the standard
        // options for help output and version information

        vector<Arg*> ops;
        ops.push_back(&concatenate);
        ops.push_back(&add);
        ops.push_back(&subtract);
        ops.push_back(&mean);

        cmd.xorAdd(ops);
        cmd.add(verbose);
        //cmd.add(reference_file);
        cmd.add(affinefrom_file);
		cmd.add(affineto_file);
        cmd.add(t1_file);
        cmd.add(t2_file);
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

    // -----------------------------------------------------------------------
    // read image headers of two image in affine registration (from which to which)
    Image::Header* affinefrom_space = NULL;
    if (!affinefrom_file.getValue().empty()) {
        affinefrom_space = nifti_read_header(affinefrom_file.getValue().c_str(), NULL, 1);
        if (affinefrom_space == NULL) {
            cerr << "Failed to read header of reference image " << affinefrom_file.getValue() << "!" << endl;
            cerr << "The file does either not exist or is not a NIfTI-1 or ANALYZE 7.5 image file." << endl;
            exit(EXIT_FAILURE);
        }
        if (fmt == Image::FORMAT_DRAMMS) *affinefrom_space = PermuteXY(*affinefrom_space);
    }
	Image::Header* affineto_space = NULL;
    if (!affineto_file.getValue().empty()) {
        affineto_space = nifti_read_header(affineto_file.getValue().c_str(), NULL, 1);
        if (affineto_space == NULL) {
            cerr << "Failed to read header of reference image " << affineto_file.getValue() << "!" << endl;
            cerr << "The file does either not exist or is not a NIfTI-1 or ANALYZE 7.5 image file." << endl;
            exit(EXIT_FAILURE);
        }
        if (fmt == Image::FORMAT_DRAMMS) *affineto_space = PermuteXY(*affineto_space);
    }

    // -----------------------------------------------------------------------
    // perform operation
    Image* t1def  = NULL;
    Image* t2def  = NULL;
    Image* outdef = NULL;

    Image::Transform t1mat;
    Image::Transform t2mat;
    Image::Transform outmat;

    char op = 'c';
    if      (concatenate.getValue()) op = 'c';
    else if (add        .getValue()) op = 'a';
    else if (subtract   .getValue()) op = 's';
    else if (mean       .getValue()) op = 'm';

    try {
        switch (op) {
            // -------------------------------------------------------------------
            // composition of transformations
            case 'c':
                // read transformations
                ok = ReadTransform(t1_file.getValue().c_str(), t1def, t1mat, fmt);
                ok = ReadTransform(t2_file.getValue().c_str(), t2def, t2mat, fmt) && ok;
                // concatenate transformations
                if (ok) {
                    // compute unified deformation field
                    if (t1def || t2def) {
                        // concatenate two deformation fields
                        if (t1def && t2def) {
                            outdef = ConcatenateTransforms(t1def, t2def);
                        // concatenate deformation and affine transformation
                        } else if (t1def && t2mat.m[3][3] != 0) {
                            // TODO Implement the following function
                            //outdef = ConcatenateTransforms(t1def, t2mat, *reference_space); // reference is image C
                            cerr << "Sorry, the concatenation of a deformable and affine transformation is not yet implemented!" << endl;
                            cerr << "This program can yet only concatenate an affine transformation from image A to B and" << endl;
                            cerr << "a deformable transformation from image B to C. Therefore, exchange the two input" << endl;
                            cerr << " arguments <T1> and <T2> if this is what you actually would like to do." << endl;
                            ok = false;
                        // concatenate affine transformation and deformation
                        } else if (t1mat.m[3][3] != 0 && t2def) {
                            if (affinefrom_space == NULL) {
                                cerr << "Require at least input of image A using option -f (--affine_from) for" << endl;
                                cerr << "concatenation of affine transformation from image A to image B and deformation" << endl;
                                cerr << "from image B to image C. If the image space of image B differs from image C," << endl;
                                cerr << "image B has to be provided as reference as well using option -t (--affine_to)." << endl;
                            }
                            if (affineto_space) {
                                outdef = ConcatenateTransforms(*affinefrom_space, *affineto_space, t1mat, t2def);
                            } else {
                                outdef = ConcatenateTransforms(*affinefrom_space, t1mat, t2def);
                            }
                        }
                        // print error if output deformation field could not be allocated
                        if (ok && outdef == NULL) {
                            cerr << "Failed to allocate memory for output deformation field!" << endl;
                            ok = false;
                        }
                    // concatenate two affine transformations
                    } else {
                        outmat = ConcatenateTransforms(t1mat, t2mat);
                    }
                } else {
                    if (t1def == NULL && t1mat.m[3][3] == 0) {
                        cerr << "File " << t1_file.getValue() << " could either not be opened or is neither a valid affine transformation nor deformation field file!" << endl;
                    }
                    if (t2def == NULL && t2mat.m[3][3] == 0) {
                        cerr << "File " << t2_file.getValue() << " could either not be opened or is neither a valid affine transformation nor deformation field file!" << endl;
                    }
                }
                break;

            // -------------------------------------------------------------------
            case 's':
                t1def = ReadDeformationField(t1_file.getValue().c_str(), fmt);
                if (t1def == NULL) {
                    cerr << "File " << t1_file.getValue() << " could either not be opened or is not a valid deformation field file!" << endl;
                    ok = false;
                }
                ok = ReadTransform(t2_file.getValue().c_str(), t2def, t2mat, fmt) && ok;
                if (ok) {
                    if (t2def == NULL) {
                        if ( (affinefrom_space == NULL)||(affineto_space == NULL) ) {
                            cerr << "Images A and B required for subtraction of affine transformation" << endl;
                            cerr << "from image A to image B from deformation from image A to image C." << endl;
                            cerr << "Use -f and -t options to specify the images A and B." << endl;
                        }
                        outdef = SubtractTransforms(*affinefrom_space, *affineto_space, t1def, t2mat);
                    } else {
                        outdef = SubtractTransforms(t1def, t2def);
                    }
                } else {
                    if (t2def == NULL && t2mat.m[3][3] == 0) {
                        cerr << "File " << t2_file.getValue() << " could either not be opened or is neither a valid affine transformation nor deformation field file!" << endl;
                        ok = false;
                    }
                }
                break;

            // -------------------------------------------------------------------
            // voxel-wise operation applied on two deformation fields
            case 'a':
            case 'm':
                t1def = ReadDeformationField(t1_file.getValue().c_str(), fmt);
                t2def = ReadDeformationField(t2_file.getValue().c_str(), fmt);
                if (t1def == NULL) {
                    cerr << "File " << t1_file.getValue() << " could either not be opened or is not a valid deformation field file!" << endl;
                    ok = false;
                }
                if (t2def == NULL) {
                    cerr << "File " << t2_file.getValue() << " could either not be opened or is not a valid deformation field file!" << endl;
                    ok = false;
                }
                if (ok && (t1def->GetNumberOfComponents() != t2def->GetNumberOfComponents() || !t2def->HasSameSizeAs(t1def))) {
                    cerr << "Input deformation fields must have identical dimensions!" << endl;
                    ok = false;
                }
                if (ok) {
                    outdef = new Image(t1def->region.nx, t1def->region.ny, t1def->region.nz,
                                       t1def->hdr.datatype, t1def->GetNumberOfComponents(),
                                       t1def->imgfmt);
                    if (outdef) {
                        outdef->CopyTransform(t1def);
                        outdef->CopyMetaData (t1def);
                    } else {
                        cerr << "Failed to allocate memory for output deformation field!" << endl;
                        ok = false;
                    }
                }
                if (ok) {
                    const int N = t1def->GetNumberOfComponents();
                    vector<float> v1(N), v2(N), vo(N);
                    for (int k = 0; k < t1def->region.nz; k++) {
                        for (int j = 0; j < t1def->region.ny; j++) {
                            for (int i = 0; i < t1def->region.nx; i++) {
                                t1def->get(i, j, k, v1);
                                t2def->get(i, j, k, v2);
                                for (int n = 0; n < N; n++) {
                                    switch (op) {
                                        case 'a': vo[n] = v1[n] + v2[n]; break;
                                        case 'm': vo[n] = (v1[n] + v2[n]) / 2; break;
                                    };
                                }
                                outdef->set(i, j, k, vo);
                            }
                        }
                    }
                }
                break;

            // -------------------------------------------------------------------
            // invalid operation
            default:
                ASSERT(false, "Invalid operation: " << op);
                ok = false;
        };
    } catch (const invalid_argument& e) {
        cerr << e.what() << endl;
        ok = false;
    }

    // -----------------------------------------------------------------------
    // write output transformation
    if (ok) {
        if (outdef) {
            ok = WriteImage(output_file.getValue().c_str(), outdef);
        } else if (outmat.m[3][3] != 0) {
            ok = WriteAffineTransform(output_file.getValue().c_str(), outmat, fmt);
        }
        if (!ok) cerr << "Failed to write resuling transformation to file " << output_file.getValue() << endl;
    }

    // -----------------------------------------------------------------------
    // clean up
    if (affinefrom_space) delete affinefrom_space;
    if (affineto_space) delete affineto_space;
    if (t1def) delete t1def;
    if (t2def) delete t2def;
    if (outdef) delete outdef;

    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
