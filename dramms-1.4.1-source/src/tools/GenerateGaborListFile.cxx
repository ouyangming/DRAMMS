/**
 * @file  GenerateGaborListFile.cxx
 * @brief generate a (text) list of Gabor feature images.
 *
 * Copyright (c) 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <string>
#include <iostream>       // cout, cerr, endl
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <dramms/basis.h> // exename(), print_contact()


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
	cout << "-------------------------------------------------" << endl;
    cout << "This program generate a text list of Gabor feature images. " << endl;
	cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_FeaturePrefix> <output_ListFileName>" << endl;
    cout << endl;
	cout << "Optional arguments:" << endl;
	cout << "\t -s <int>     : the number of scales in Gabor features (default: 3)" << endl;
	cout << "\t -o <int>     : the number of orientations in Gabor features (default: 3)" <<endl;
	cout << "\t -H           : horizental features only (default: off, so it will generate both horizental and vertical features.)" << endl;
	cout << "\t -R           : real part of Gabor features only  (default: off)" << endl;
	cout << endl;
    cout << "Example: " << endl;
	cout << "\t " << exec_name <<" A_level3 Gabor_A_level3.lst -s3 -o3 -H " << endl;
	cout << endl;
    print_contact();
}





int main(int argc,char *argv[])
{
  bool horizentalOnly = false;
  bool realPartOnly = false;
  int numScales = 3;
  int numOrientations = 3; // default
  int numFeatures;
  int scaleIndex, orientationIndex;
  FILE *fp;
  
  
  int c=-1;
  while((c=getopt(argc,argv,"s:o:HR")) != -1)
    {
      switch(c)
		{
			case 's':
				sscanf(optarg,"%d",&numScales);
				break;
			
			case 'o':
				sscanf(optarg,"%d", &numOrientations);
				break;
				
			case 'H':
				horizentalOnly = true;
				break;
				
			case 'R':
				realPartOnly = true;
				break;
				
			default:
			    break;
		}
    }
	
  argc -= optind;
  argv += optind;
  if (argc != 2) {
	print_help();
    exit(1);
  }

  const char* GaborPrefix=argv[0];
  const char* GaborFileName=argv[1];
  
  
  fp = fopen(GaborFileName, "w");
  if ( (horizentalOnly==false) && (realPartOnly==false) )
    {
	numFeatures = 4*numScales*numOrientations;
	fprintf(fp, "%d\n", numFeatures);
	for (scaleIndex=0; scaleIndex<numScales; scaleIndex++)
	  for (orientationIndex=0; orientationIndex<numOrientations; orientationIndex++)
	    {
		printf("%s_3dHori_F_imag.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		printf("%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
	    fprintf(fp, "%s_3dHori_F_imag.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		fprintf(fp, "%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		}
		
	for (scaleIndex=0; scaleIndex<numScales; scaleIndex++)
	  for (orientationIndex=0; orientationIndex<numOrientations; orientationIndex++)
	    {
		printf("%s_3dVert_F_imag.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		printf("%s_3dVert_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
	    fprintf(fp, "%s_3dVert_F_imag.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		fprintf(fp, "%s_3dVert_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		}
	fclose(fp);
	}
  else if ( (horizentalOnly == false) && (realPartOnly==true) )
    {
	numFeatures = 2*numScales*numOrientations;
	fprintf(fp, "%d\n", numFeatures);
	for (scaleIndex=0; scaleIndex<numScales; scaleIndex++)
	  for (orientationIndex=0; orientationIndex<numOrientations; orientationIndex++)
	    {
	    printf("%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		printf("%s_3dVert_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
	    fprintf(fp, "%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		fprintf(fp, "%s_3dVert_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		}
	fclose(fp);
	}  
  else if ( (horizentalOnly == true) && (realPartOnly==false) )
    {
	numFeatures = 2*numScales*numOrientations;
	fprintf(fp, "%d\n", numFeatures);
	for (scaleIndex=0; scaleIndex<numScales; scaleIndex++)
	  for (orientationIndex=0; orientationIndex<numOrientations; orientationIndex++)
	    {
	    printf("%s_3dHori_F_imag.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		printf("%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
	    fprintf(fp, "%s_3dHori_F_imag.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		fprintf(fp, "%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		}
	fclose(fp);
	}
  else if ( (horizentalOnly == true) & (realPartOnly==true) )
    {
	numFeatures = numScales*numOrientations;
	fprintf(fp, "%d\n", numFeatures);
	for (scaleIndex=0; scaleIndex<numScales; scaleIndex++)
	  for (orientationIndex=0; orientationIndex<numOrientations; orientationIndex++)
	    {
	    printf("%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
	    fprintf(fp, "%s_3dHori_F_real.%d_%d.nii.gz\n", GaborPrefix, scaleIndex, orientationIndex);
		}
	fclose(fp);
	}

  printf("\nGabor list has been written into file %s\n\n", GaborFileName);	
}

