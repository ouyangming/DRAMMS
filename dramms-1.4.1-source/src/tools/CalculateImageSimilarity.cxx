/**
 * @file  CalculateJacobianMap.cxx
 * @brief Calculate jacobian determinante map.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <string>
#include <iostream> // cout, cerr, endl
#include <math.h>   // M_PI, exp(), pow(), sqrt(),...
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <common/mvcd.h>
#include <common/imageio.h>
#include <common/cres.h>
#include <common/matrix.h>
#include <common/general.h>


#include <dramms/basis.h> // exename(), print_contact()


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;


#define		CC		1
#define 	MI		2
#define		NMI		3
#define		SSD		4


// ===========================================================================
// help
// ===========================================================================

// ---------------------------------------------------------------------------
void print_help()
{
    string exec_name = exename();
	cout << "-------------------------------------------------" << endl;
    cout << "This program calculates similarity metrics (CC/MI/NMI/SSD) of two input images. " << endl;
	cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input image 1> <input image 2>" << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
	cout << "  <input image 1>    file of the first input image." << endl;
    cout << "  <input image 2>    file of the second input image." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
	cout << "  -C   : output correlation coefficient (CC)" <<endl;
	cout << "  -M   : output mutual information (MI)" << endl;
	cout << "  -N   : output normalized mutual information (NMI)" << endl;
	cout << "  -D   : output sum of squared difference (SSD)" << endl;
	cout << endl;
	cout << "Example: " << endl;
	cout << "  " << exec_name <<" A.img B.nii.gz -M" << endl;
	cout << endl;
    print_contact();
}


// -------------------------------------------------
float computeCorrelationCoefficient(Image *source, Image *target)
{
	double sumIntensityA = 0.0;
	double sumIntensityB = 0.0;
	double sumIntensityAxB = 0.0;
	double sumIntensityAsquared = 0.0;
	double sumIntensityBsquared = 0.0;
	double denominator, numerator;
	double cc;
   
	int i,j,k;
	
	Fvector3d resampleRatio;
	resampleRatio.x = static_cast<double>(source->region.nx) / static_cast<double>(target->region.nx) ;
	resampleRatio.y = static_cast<double>(source->region.ny) / static_cast<double>(target->region.ny) ;
	resampleRatio.z = static_cast<double>(source->region.nz) / static_cast<double>(target->region.nz) ;
	
	   // calculate denominator and numerator
	   double intensity;
	   for (k=0; k<target->region.nz; k++)
		 for (i=0; i<target->region.nx; i++)
		   for (j=0; j<target->region.ny; j++)
			 {
			 // resample image1 into target space
			 intensity = source->value(static_cast<float>(i) * resampleRatio.x,
                                       static_cast<float>(j) * resampleRatio.y,
                                       static_cast<float>(k) * resampleRatio.z);
			
			 sumIntensityA        += intensity;
			 sumIntensityB        += static_cast<double>(target->img.uc[k][i][j]);
			 sumIntensityAxB      += static_cast<double>(intensity*target->img.uc[k][i][j]);
			 sumIntensityAsquared += static_cast<double>(intensity*intensity);
			 sumIntensityBsquared += static_cast<double>(target->img.uc[k][i][j]*target->img.uc[k][i][j]);
			 }

      int totalNumVoxels = target->region.nx * target->region.ny * target->region.nz;		
	  numerator          = static_cast<double>(totalNumVoxels)*sumIntensityAxB - sumIntensityA*sumIntensityB;
	  denominator        = sqrt( static_cast<double>(totalNumVoxels)*sumIntensityAsquared - sumIntensityA*sumIntensityA ) * sqrt( static_cast<double>(totalNumVoxels)*sumIntensityBsquared - sumIntensityB*sumIntensityB );

	  if (denominator!=0.0)
		cc = numerator/denominator;
	  
	  printf("\nCC = %f\n\n", cc);
	  
	  return cc;
}



// -------------------------------------------------
float computeMutualInformation(Image *source, Image *target)
{
	  int numBins=64; //default
	  float interval = 256.0/static_cast<float>(numBins);
	  float hist1[numBins], hist2[numBins], jointHist[numBins][numBins];
	  int i,j,k;
	  
	  
	  Fvector3d resampleRatio;
	  resampleRatio.x = static_cast<float>(source->region.nx) / static_cast<float>(target->region.nx) ;
	  resampleRatio.y = static_cast<float>(source->region.ny) / static_cast<float>(target->region.ny) ;
	  resampleRatio.z = static_cast<float>(source->region.nz) / static_cast<float>(target->region.nz) ;
	
	   // initialize histograms
	   for (i=0; i<numBins; i++)
		 {
		 hist1[i] = 0.0;
		 hist2[i] = 0.0;
		 for (j=0; j<numBins; j++)
		   jointHist[i][j] = 0.0;
		 }
		 
		 
	   // calculate histogram
	   float intensity;
	   for (k=0; k<target->region.nz; k++)
		 for (i=0; i<target->region.nx; i++)
		   for (j=0; j<target->region.ny; j++)
			 {
			 // resample image1 into image2 space
			 intensity = source->value(static_cast<float>(i) * resampleRatio.x,
                                       static_cast<float>(j) * resampleRatio.y,
                                       static_cast<float>(k) * resampleRatio.z);

			 hist1[static_cast<int>(intensity/interval)] += 1.0;
			 hist2[static_cast<int>(static_cast<float>(target->img.uc[k][i][j])/interval)] += 1.0;
			 jointHist[static_cast<int>(intensity/interval)][static_cast<int>(static_cast<float>(target->img.uc[k][i][j])/interval)] += 1.0;
			 }
			 
	  // calculate entropies
	  int totalNumVoxels = target->region.nx * target->region.ny * target->region.nz;
	  float entropy1 = 0.0;
	  float entropy2 = 0.0;
	  float jointEntropy = 0.0;
	  for (i=0; i<numBins; i++)
		{
		hist1[i] /= static_cast<float>(totalNumVoxels);
		hist2[i] /= static_cast<float>(totalNumVoxels);
		
		if (hist1[i]>0.0)
			entropy1 -= hist1[i]*log(hist1[i]);
		if (hist2[i]>0.0)
			entropy2 -= hist2[i]*log(hist2[i]);
			
			
		for (j=0; j<numBins; j++)
		  {
		  jointHist[i][j] /= static_cast<float>(totalNumVoxels);
		  if (jointHist[i][j]>0.0)
			jointEntropy -= jointHist[i][j]*log(jointHist[i][j]);
		  }
		}
		
	  
	  // calculate MI and normalized MI
	  float mi = entropy1 + entropy2 - jointEntropy;
	  
	  printf("\nentropy1 = %f, entropy2 = %f, jointEntropy = %f\n", entropy1, entropy2, jointEntropy);
	  printf("\nMI_unnormalized = %f\n", mi);
	  
	  return mi;
}




// -------------------------------------------------
float computeNormalizedMutualInformation(Image *source, Image *target)
{
	  int numBins=64; //default
	  float interval = 256.0/static_cast<float>(numBins);
	  float hist1[numBins], hist2[numBins], jointHist[numBins][numBins];
	  int i,j,k;
	  
	  
	  Fvector3d resampleRatio;
	  resampleRatio.x = static_cast<float>(source->region.nx) / static_cast<float>(target->region.nx) ;
	  resampleRatio.y = static_cast<float>(source->region.ny) / static_cast<float>(target->region.ny) ;
	  resampleRatio.z = static_cast<float>(source->region.nz) / static_cast<float>(target->region.nz) ;
	
	   // initialize histograms
	   for (i=0; i<numBins; i++)
		 {
		 hist1[i] = 0.0;
		 hist2[i] = 0.0;
		 for (j=0; j<numBins; j++)
		   jointHist[i][j] = 0.0;
		 }
		 
		 
	   // calculate histogram
	   float intensity;
	   for (k=0; k<target->region.nz; k++)
		 for (i=0; i<target->region.nx; i++)
		   for (j=0; j<target->region.ny; j++)
			 {
			 // resample image1 into image2 space
			 intensity = source->value(static_cast<float>(i) * resampleRatio.x,
                                       static_cast<float>(j) * resampleRatio.y,
                                       static_cast<float>(k) * resampleRatio.z);

			 hist1[static_cast<int>(intensity/interval)] += 1.0;
			 hist2[static_cast<int>(static_cast<float>(target->img.uc[k][i][j])/interval)] += 1.0;
			 jointHist[static_cast<int>(intensity/interval)][static_cast<int>(static_cast<float>(target->img.uc[k][i][j])/interval)] += 1.0;
			 }
			 
	  // calculate entropies
	  int totalNumVoxels = target->region.nx * target->region.ny * target->region.nz;
	  float entropy1 = 0.0;
	  float entropy2 = 0.0;
	  float jointEntropy = 0.0;
	  for (i=0; i<numBins; i++)
		{
		hist1[i] /= static_cast<float>(totalNumVoxels);
		hist2[i] /= static_cast<float>(totalNumVoxels);
		
		if (hist1[i]>0.0)
			entropy1 -= hist1[i]*log(hist1[i]);
		if (hist2[i]>0.0)
			entropy2 -= hist2[i]*log(hist2[i]);
			
			
		for (j=0; j<numBins; j++)
		  {
		  jointHist[i][j] /= static_cast<float>(totalNumVoxels);
		  if (jointHist[i][j]>0.0)
			jointEntropy -= jointHist[i][j]*log(jointHist[i][j]);
		  }
		}
		
	  
	  float normmi = (entropy1+entropy2)/jointEntropy;
	  printf("\nentropy1 = %f, entropy2 = %f, jointEntropy = %f\n", entropy1, entropy2, jointEntropy);
	  printf("\nMI_normalized = %f\n\n", normmi); 
	  
	  return normmi;
}




// -------------------------------------------------
float computeSumSquaredDifference(Image *source, Image *target)
{
	  int i,j,k;
	  
	  Fvector3d resampleRatio;
	  resampleRatio.x = static_cast<float>(source->region.nx) / static_cast<float>(target->region.nx) ;
	  resampleRatio.y = static_cast<float>(source->region.ny) / static_cast<float>(target->region.ny) ;
	  resampleRatio.z = static_cast<float>(source->region.nz) / static_cast<float>(target->region.nz) ;
	
	   float intensity;
	   float ssd=0;
	   float diff;
	   for (k=0; k<target->region.nz; k++)
		 for (i=0; i<target->region.nx; i++)
		   for (j=0; j<target->region.ny; j++)
			 {
			 // resample image1 into image2 space
			 intensity = source->value(static_cast<float>(i) * resampleRatio.x,
                                       static_cast<float>(j) * resampleRatio.y,
                                       static_cast<float>(k) * resampleRatio.z);

			 diff = static_cast<float>(intensity - target->img.uc[k][i][j]);
			 ssd += diff*diff;
			 }
			 
	  printf("\nSSD = %f\n\n", ssd);
	  return ssd;
}


// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc,char *argv[])
{
    if (argc == 1) {
        print_help();
        exit(1);
    }

  int metric=CC; // default
  
  
  int c=-1;
  while((c=getopt(argc,argv,"CMNDv")) != -1)
    {
      switch(c)
	{
	case 'C':
	    metric=CC;
		break;
	
	case 'M':
		metric=MI;
		break;
		
	case 'N':
		metric=NMI;
		break;
		
	case 'D':
		metric=SSD;
		break;
		
    case 'v':
        // ignore
        break;

    case 'h':
        print_help();
        exit(0);

	default:
        // error message printed by getopt() already
        exit(EXIT_FAILURE);
		break;
	}
    }

  
  
  argc -= optind;
  argv += optind;

    if (argc < 2) {
        cerr << "Not all required arguments specified!" << endl;
        cerr << "See help (-h option) for a list of required arguments." << endl;
        exit(1);
    }
    if (argc > 2) {
        cerr << "Too many arguments specified!" << endl;
        cerr << "See help (-h option) for usage information." << endl;
        exit(1);
    }

  Image* imageA = ReadNiftiImage(argv[0], DT_UNSIGNED_CHAR);
  Image* imageB = ReadNiftiImage(argv[1], DT_UNSIGNED_CHAR);
  
  switch (metric)
	{
	case CC:
		computeCorrelationCoefficient(imageA, imageB);
		break;
		
	case MI:
		computeMutualInformation(imageA, imageB);
		break;
		
	case NMI:
		computeNormalizedMutualInformation(imageA, imageB);
		break;
		
	case SSD:
		computeSumSquaredDifference(imageA, imageB);
		break;
	}
	
	
  delete imageA;
  delete imageB;
}



