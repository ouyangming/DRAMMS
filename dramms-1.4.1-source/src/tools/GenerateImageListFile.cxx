/**
 * @file  GenerateImageListFile.cxx
 * @brief generate a text file that lists intensity image as feature image to be used in deformation.
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
    cout << "This program generate a text file that lists intensity image as feature image to be used in deformation. " << endl;
	cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_imageName> <output_ListFileName>" << endl;
    cout << endl;
	cout << "Example: " << endl;
	cout << "  " << exec_name <<" A.img ImageFeatures.lst " << endl;
	cout << endl;
    print_contact();
}




int main(int argc,char *argv[])
{
  FILE *fp;
  
  argc -= optind;
  argv += optind;
  if (argc != 2) {
	print_help();
    exit(1);
  }

  const char* ImageName=argv[0];
  const char* ListFileName=argv[1];
  
  
  fp = fopen(ListFileName, "w");
  fprintf(fp, "1\n");
  fprintf(fp, "%s\n", ImageName);
  fclose(fp);
	
  printf("\nImage list has been written into file %s\n\n", ListFileName);	
}

