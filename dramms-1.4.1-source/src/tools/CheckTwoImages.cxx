/**
 * @file  CheckTwoImages.cxx
 * @brief Checks if two input images share the same image size and voxel size.
 *
 * Copyright (c) 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */
 
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <common/imageio.h> 
#include <common/general.h>

#include <dramms/basis.h> // exename(), print_contact()


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;


// ATTENTION: The exit code of this program is the number of dimensions of
//            the input images! In case of an error or if the images are
//            not defined on the same image grid, the value 0 is returned!
//            On error, the value 255 is returned.


// ===========================================================================
// help
// ===========================================================================

// ---------------------------------------------------------------------------
void print_help()
{
    string exec_name = exename();
    cout << "-------------------------------------------------" << endl;
    cout << "This program checks two intensity images or two features images. For intensity images (by default), the program checks whether two images share the same image and voxel size. In this case, the program returns 0 if two images are not in the same space, 2 if they are in the same 2D space, and 3 if they are in the same 3D space. For feature images (-f option), the program checks whether two feature images have roughly the same number of foreground voxels. It returns 100 if two images do have roughly the same number of foreground voxels, and 200 otherwise" << endl;
    cout << " -h  help" <<endl;
	cout << " -f  check feature images" << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << exec_name << " A.img  B.img" << endl;
    cout << endl;
    print_contact();
}


// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc,char *argv[])
{
	bool checkFeatureImage = false; // by default, check intensity images
 	
    // show usage if no arguments were given
    if (argc == 1) {
        print_help();
        exit(1);
    }

    // parse arguments
    int c = -1;
    while((c = getopt(argc, argv, "hf")) != -1) {
        switch(c) {
            case 'h':
                print_help();
                exit(0);

			case 'f':
				checkFeatureImage=true;
				break;
				
            default:
                break;
        }
    }
  
    argc -= optind;
    argv += optind;

    if (argc < 2) {
        print_help();
        cerr << "Not all required arguments specified!" << endl;
        exit(255);
    }
    if (argc > 2) {
        print_help();
        cerr << "Too many arguments specified!" << endl;
        exit(255);
    }

	
    nifti_1_header *nhdr1 = nifti_read_header(argv[0], NULL, 1);
    if (nhdr1 == NULL) {
        cerr << "Failed to read header of image " << argv[0] << endl;
        exit(255);
    }

    nifti_1_header *nhdr2 = nifti_read_header(argv[1], NULL, 1);
    if (nhdr2 == NULL) {
        cerr << "Failed to read header of image " << argv[1] << endl;
        free(nhdr1);
        exit(255);
    }

	
    if (nhdr1->dim[3] < 1) {
        nhdr1->dim   [3] = 1;
        nhdr1->pixdim[3] = 1.0f;
    }
    if (nhdr2->dim[3] < 1) {
        nhdr2->dim   [3] = 1;
        nhdr2->pixdim[3] = 1.0f;
    }

	if (checkFeatureImage==false)
	{
		bool ShareImageSize = false;
		if ((nhdr1->dim[1] == nhdr2->dim[1]) &&
			(nhdr1->dim[2] == nhdr2->dim[2]) &&
			(nhdr1->dim[3] == nhdr2->dim[3])) ShareImageSize = true;

		bool ShareVoxelSize = false;
		if ((nhdr1->pixdim[1] == nhdr2->pixdim[1]) && 
			(nhdr1->pixdim[2] == nhdr2->pixdim[2]) &&
			(nhdr1->pixdim[3] == nhdr2->pixdim[3])) ShareVoxelSize = true;

		if ( (ShareImageSize==true)&&(ShareVoxelSize==true) ) {
			//cout << "Two images share the same image and voxel sizes." << endl;
			cout << "Image size: (" << nhdr1->dim   [1] << ", " << nhdr1->dim   [2] << ", " << nhdr1->dim   [3] << ")" << endl;
			cout << "Voxel size: (" << nhdr1->pixdim[1] << ", " << nhdr1->pixdim[2] << ", " << nhdr1->pixdim[3] << ")" << endl;
			if (nhdr1->dim[3] == 1) {
				free(nhdr1);
				free(nhdr2);
				exit(2);
			} else if (nhdr1->dim[3] > 1) {
            free(nhdr1);
            free(nhdr2);
            exit(3);
			}
        }
		else
		{
			cout << "Two images vary by image size and/or voxel size." << endl;
			cout << "Image size: (" << nhdr1->dim   [1] << ", " << nhdr1->dim   [2] << ", " << nhdr1->dim   [3] << ") for the first image, (" << nhdr2->dim   [1] << ", " << nhdr2->dim   [2] << ", " << nhdr2->dim   [3] << ") for the second image. " << endl;
			cout << "Voxel size: (" << nhdr1->pixdim   [1] << ", " << nhdr1->pixdim   [2] << ", " << nhdr1->pixdim   [3] << ") for the first image, (" << nhdr2->pixdim   [1] << ", " << nhdr2->pixdim   [2] << ", " << nhdr2->pixdim   [3] << ") for the second image. " << endl;
		}
		free(nhdr1);
		free(nhdr2);
		exit(0);
    }
	else // check feature images
	{
		Image* imageA = ReadImage(argv[0]);
		Image* imageB = ReadImage(argv[1]);
		unsigned char ***imgA, ***imgB;
		imgA = imageA->img.uc;
		imgB = imageB->img.uc;
		int dimx = imageA->region.nx; 
		int dimy = imageA->region.ny; 
		int dimz = imageA->region.nz; 
		
		int i,j,k;
		int minjA, maxjA, minjB, maxjB;
		int countA=0;
		int countB=0;
		int countAfilled=0;
		int countBfilled=0;
		for (k=0;k<dimz;k++)
		  for (i=0;i<dimx;i++)
			{
			minjA=dimy-1;
			minjB=dimy-1;
			maxjA=0;
			maxjB=0;
		    for (j=0;j<dimy;j++)
				{
				if (imgA[k][i][j]>0) {
					countA++;
					minjA=MIN(minjA,j); 
					maxjA=MAX(maxjA,j);
					}
				if (imgB[k][i][j]>0) {
					countB++;
					minjB=MIN(minjB,j); 
					maxjB=MAX(maxjB,j);
					}
				}				
			countAfilled += MAX(maxjA-minjA, 0);
			countBfilled += MAX(maxjB-minjB, 0);
			}
		 // printf("countA=%d, countB=%d, ratio=%f\ncountAfilled=%d, countBfilled=%d, ratio=%f\n", countA, countB, (float)countA/(float)countB, countAfilled, countBfilled, (float)countAfilled/(float)countBfilled);
		
		if ( ( (countA>1.25*countB) || (countA<0.8*countB) )&&( (countAfilled>1.15*countBfilled) || (countAfilled<0.87*countBfilled) ) ){
            delete imageA;
			delete imageB;
            exit(200);
			}
		else
			{
			delete imageA;
			delete imageB;
            exit(100);
			}		
	}
}
