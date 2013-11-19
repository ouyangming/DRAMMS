/**
 * @file  CombineFeatureLists.cxx
 * @brief combine two feature lists into one.
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
    cout << "This program combines two input feature lists into a single list. " << endl;
	cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_ListFileName1> <input_ListFileName2> <output_CombinedFileName>" << endl;
    cout << endl;
    cout << "Example: " << endl;
	cout << "  " << exec_name <<" GaborFeature1.lst GaborFeature2.lst GaborFeature_all.lst " << endl;
	cout << endl;
    print_contact();
}



int main(int argc,char *argv[])
{
  int numFeatures1, numFeatures2, numFeaturesCombined;
  int idx;
  
  char featureName[1000];
  FILE *fp1, *fp2, *fpCombined;
  
  argc -= optind;
  argv += optind;
  if (argc != 3) {
	print_help();
    exit(1);
  }

  const char* fileName1 = argv[0];
  const char* fileName2 = argv[1];
  const char* fileNameCombined = argv[2];
  
  
  if (NULL==(fp1=fopen(fileName1,"rb"))){
       printf("File %s doesn't exist!\n",fileName1);
	   exit(1);
	   }
  if (NULL==(fp2=fopen(fileName2,"rb"))){
       printf("File %s doesn't exist!\n",fileName2);
	   exit(1);
	   }
  fscanf(fp1,"%d",&numFeatures1);  
  fscanf(fp2,"%d",&numFeatures2);
  numFeaturesCombined = numFeatures1+numFeatures2;
  fpCombined = fopen(fileNameCombined, "w");
  fprintf(fpCombined, "%d\n", numFeaturesCombined);
  
  
  for (idx=1;idx<=numFeatures1;idx++)
    {
	fscanf(fp1,"%s",featureName);
	fprintf(fpCombined, "%s\n", featureName);
	}
  fclose(fp1);
  for (idx=1;idx<=numFeatures2;idx++)
	{
	fscanf(fp2,"%s",featureName);
	fprintf(fpCombined, "%s\n", featureName);
	}
  fclose(fp2);
  fclose(fpCombined);
  
  printf("\n%d features in file %s\n", numFeatures1, fileName1);	
  printf("%d features in file %s\n", numFeatures2, fileName2);	
  printf("=>\n");
  printf("%d featuers combined in file %s\n\n", numFeaturesCombined, fileNameCombined);
}



