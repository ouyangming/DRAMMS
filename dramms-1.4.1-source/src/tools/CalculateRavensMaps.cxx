/**
 * @file  CalculateRavensMaps.cxx
 * @brief Calculate RAVENS map from deformation field.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <iostream>
#include <string>
#include <set>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <common/cres.h>
#include <common/general.h>
#include <common/imageio.h> 

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
    cout << "This program takes two inputs 1) label/segmented image in subject space and 2) deformation field (generated when registering template image to subject image), and outputs RAVENS maps (in signed short datatype) in template space (one RAVENS map per tissue type, or per ROI). Up to 5 RAVENS maps for up to 5 tissue types can be calculated. " << endl;
    cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <labelImage> <templateImage> <DefField> <RAVENS_prefix>" << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
    cout << "  <labelImage>      : label image in subject space. Note that the label image must be in uint8 (aka byte, or uchar) datatype; other input images can be in any datatype" << endl;
    cout << "  <templateImage>   : template image, where RAVENS map will reside." << endl; 
    cout << "  <DefField>        : deformation field (generated when registering template to subject image)." << endl;
    cout << "  <RAVENS_prefix>   : prefix for all RAVENS maps." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
    cout << "  -m <int>,<int>,<int>,<int>,<int>   : labels of up to 5 ROIs where RAVENS maps will be calculated (will output one RAVENS map for each ROI)." << endl;
    cout << "  -f <int>          : scale factor (default: 1000)" << endl;
    cout << "  -h                : help; usage of this program." << endl;
    cout << "Example:" << endl;
    cout << "  " << exec_name << " labelImg.img DF.def RAVENSprefix -m100,200,255" << endl;
    cout << endl;
    print_contact();
}

// ===========================================================================
// auxiliary functions
// ===========================================================================

// ---------------------------------------------------------------------------
void maskoutImage(unsigned char ***Img, Ivector3d imageSize, int *foregroundIntensity, int num)
{
  int i,j,k,n;
  unsigned char ***labelImg;
  labelImg = UCalloc3d(imageSize.x, imageSize.y,imageSize.z);
  
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
        {
        labelImg[k][i][j]=0;
        for (n=0;n<num;n++)
        {
        if (Img[k][i][j]==foregroundIntensity[n])
           labelImg[k][i][j]=n+1;
        }
        }
  
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
        Img[k][i][j]=labelImg[k][i][j];
        
  UCfree3d(labelImg, imageSize.z,imageSize.x);
}

// ---------------------------------------------------------------------------
float InterpolatedIntensity(float ii, float jj, float kk, unsigned char ***Img, int x_size, int y_size, int z_size)
{
  float CurrentV ;
  float b,c,d, b1,c1,d1;
  int   ni,nj,nk, niP1,njP1,nkP1;
  float GreyValue ;


  ni = (int)ii ;
  nj = (int)jj ;
  nk = (int)kk ;
  
  niP1 = ni+1 ;
  njP1 = nj+1 ;
  nkP1 = nk+1 ;
  
  GreyValue = 0;
  if(ni>=0 && ni<x_size-1  &&  nj>=0 && nj<y_size-1  &&  nk>=0 && nk<z_size-1 )
    {
      b = ii-ni ;        b1 = 1.-b ;
      c = jj-nj ;        c1 = 1.-c ;
      d = kk-nk ;        d1 = 1.-d ;

      CurrentV = ( d1*(Img[nk][ni][nj]*(b1*c1) + Img[nk][niP1][nj]*(b*c1) + Img[nk][ni][njP1]*(b1*c) + Img[nk][niP1][njP1]*(b*c)) + d*(Img[nkP1][ni][nj]*(b1*c1) + Img[nkP1][niP1][nj]*(b*c1) + Img[nkP1][ni][njP1]*(b1*c) + Img[nkP1][niP1][njP1]*(b*c)) )/( d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c)) ) ; 

     
      GreyValue = CurrentV ;
    }
  
  if ( (ni==x_size-1 && nj>=0 && nj<y_size-1 && nk>=0 && nk<z_size-1) || (ni>=0 && ni<x_size-1 && nj==y_size-1 && nk>=0 && nk<z_size-1)  || (ni>=0 && ni<x_size-1 && nj>=0 && nj<y_size-1 && nk==z_size-1) )
      GreyValue = (float)Img[nk][ni][nj] ;

  return GreyValue ;
}

// ---------------------------------------------------------------------------
void InterpolatedDisplacement(Fvector3d *Displace_subVoxel, float ii, float jj, float kk, Fvector3d ***DeformFld, int x_size, int y_size, int z_size)
{
  float CurrentV ;
  float b,c,d, b1,c1,d1;
  int   ni,nj,nk, niP1,njP1,nkP1, GreyValue ;


  ni = (int)ii ;
  nj = (int)jj ;
  nk = (int)kk ;
  
  niP1 = ni+1 ;
  njP1 = nj+1 ;
  nkP1 = nk+1 ;
  
  if(ni>=0 && ni<x_size-1  &&  nj>=0 && nj<y_size-1  &&  nk>=0 && nk<z_size-1 )
    {
      b = ii-ni ;        b1 = 1.-b ;
      c = jj-nj ;        c1 = 1.-c ;
      d = kk-nk ;        d1 = 1.-d ;

      (*Displace_subVoxel).x = ( d1*(DeformFld[nk][ni][nj].x*(b1*c1) + DeformFld[nk][niP1][nj].x*(b*c1) + DeformFld[nk][ni][njP1].x*(b1*c) + DeformFld[nk][niP1][njP1].x*(b*c)) + d*(DeformFld[nkP1][ni][nj].x*(b1*c1) + DeformFld[nkP1][niP1][nj].x*(b*c1) + DeformFld[nkP1][ni][njP1].x*(b1*c) + DeformFld[nkP1][niP1][njP1].x*(b*c)) )/( d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c)) ) ; 

      (*Displace_subVoxel).y = ( d1*(DeformFld[nk][ni][nj].y*(b1*c1) + DeformFld[nk][niP1][nj].y*(b*c1) + DeformFld[nk][ni][njP1].y*(b1*c) + DeformFld[nk][niP1][njP1].y*(b*c)) + d*(DeformFld[nkP1][ni][nj].y*(b1*c1) + DeformFld[nkP1][niP1][nj].y*(b*c1) + DeformFld[nkP1][ni][njP1].y*(b1*c) + DeformFld[nkP1][niP1][njP1].y*(b*c)) )/( d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c)) ) ; 

      (*Displace_subVoxel).z = ( d1*(DeformFld[nk][ni][nj].z*(b1*c1) + DeformFld[nk][niP1][nj].z*(b*c1) + DeformFld[nk][ni][njP1].z*(b1*c) + DeformFld[nk][niP1][njP1].z*(b*c)) + d*(DeformFld[nkP1][ni][nj].z*(b1*c1) + DeformFld[nkP1][niP1][nj].z*(b*c1) + DeformFld[nkP1][ni][njP1].z*(b1*c) + DeformFld[nkP1][niP1][njP1].z*(b*c)) )/( d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c)) ) ; 
    }
  
  if((ni==x_size-1 && nj>=0 && nj<y_size-1 && nk>=0 && nk<z_size-1) || (ni>=0 && ni<x_size-1 && nj==y_size-1 && nk>=0 && nk<z_size-1) || (ni>=0 && ni<x_size-1 && nj>=0 && nj<y_size-1 && nk==z_size-1))
    {
      (*Displace_subVoxel).x = DeformFld[nk][ni][nj].x ;
      (*Displace_subVoxel).y = DeformFld[nk][ni][nj].y ;
      (*Displace_subVoxel).z = DeformFld[nk][ni][nj].z ;
    }
}

// ---------------------------------------------------------------------------
void DistributingVolume(float weight, float ii, float jj, float kk, float ***Mass, int x_size, int y_size, int z_size)
{
  float CurrentV ;
  float b,c,d, b1,c1,d1, combined;
  int   ni,nj,nk, niP1,njP1,nkP1, GreyValue ;


  ni = (int)ii ;
  nj = (int)jj ;
  nk = (int)kk ;
  
  niP1 = ni+1 ;
  njP1 = nj+1 ;
  nkP1 = nk+1 ;
  
  
  if(ni>=0 && ni<x_size-1  &&  nj>=0 && nj<y_size-1  &&  nk>=0 && nk<z_size-1 )
    {
      b = ii-ni ;        b1 = 1.-b ;
      c = jj-nj ;        c1 = 1.-c ;
      d = kk-nk ;        d1 = 1.-d ;
      
      
      combined = ( d1*((b1*c1)+(b*c1)+(b1*c)+(b*c)) + d*((b1*c1)+(b*c1)+(b1*c)+(b*c)) ) ;

      Mass[nk][ni][nj]      += weight*d1*(b1*c1)/combined ;
      Mass[nk][niP1][nj]    += weight*d1*(b*c1)/combined ;
      Mass[nk][ni][njP1]    += weight*d1*(b1*c)/combined ;
      Mass[nk][niP1][njP1]  += weight*d1*(b*c)/combined ;
      Mass[nkP1][ni][nj]    += weight*d*(b1*c1)/combined ;
      Mass[nkP1][niP1][nj]  += weight*d*(b*c1)/combined ;
      Mass[nkP1][ni][njP1]  += weight*d*(b1*c)/combined ;
      Mass[nkP1][niP1][njP1]+= weight*d*(b*c)/combined ;
    }
  
  if((ni==x_size-1 && nj>=0 && nj<y_size-1 && nk>=0 && nk<z_size-1) || (ni>=0 && ni<x_size-1 && nj==y_size-1 && nk>=0 && nk<z_size-1) || (ni>=0 && ni<x_size-1 && nj>=0 && nj<y_size-1 && nk==z_size-1))
    Mass[nk][ni][nj]      += weight*1 ;
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

  Ivector3d     imageSizeSubj, imageSizeTemp;
  Fvector3d     voxelSizeSubj, voxelSizeTemp;
  float 	voxelSizeRatioS2T;
  float     ***Mass;
  short         ***RAVENSMap;
  Fvector3d ***DeformFld;
  int           i,j,k,ravensIndex;
  int       x,y,z;
  int       scaleFactor;
  float         level, weight, thres ;
  int           sampleNum;
  bool DeformationSmoothOrNot=false;
  bool DRAMMSDeformationOrNot=true;
  bool templateSizeGivenOrNot=false;
  bool maskOrNot = false;
  int foregroundIntensity[5] = {0, 0, 0, 0, 0};
  int numRAVENSMaps=1;
  int smoothWindowSizeXY, smoothWindowSizeZ;
  float smoothKernelXY, smoothKernelZ;
  float volumeOriginal=0.0f;
  float volumeRAVENS=0.0f;
  float warpedX, warpedY, warpedZ;
  Fvector3d  Subj_subVoxel, Displace_subVoxel ;
  float ii,jj,kk;
  float volumePreservePercentage;
  float scaledMass;
  float massOverflow;
  
  
  Image* labelimage  = NULL;
  Image* tempimage   = NULL;
  Image* deformation = NULL;
  Image* ravensimage = NULL;

  
  sampleNum = 2;      // was 1 in earlier version in SBIA
  float interval = 1.0/(sampleNum*2+1) ;
  scaleFactor = 1000;  // to scale original RAVENS values (typically <10) into short values
  
  
  int c=-1;
  while ( (c=getopt(argc,argv,"m:f:h"))!=-1 )
    {
      switch(c)
    {
    case 'm':   // can calculate RAVENS for up to five ROIs/VOIs
        sscanf(optarg, "%d,%d,%d,%d,%d", &foregroundIntensity[0],&foregroundIntensity[1],&foregroundIntensity[2],&foregroundIntensity[3],&foregroundIntensity[4]);
        maskOrNot = true;
        break;
        
    case 'f':
        sscanf(optarg, "%d", &scaleFactor);
        break;
        
    case 'h':
        print_help();
        exit(0);

    default:
        break;
    }
    }

  argc -= optind;
  argv += optind;

    if (argc < 4) {
        cerr << "Not all required arguments specified!" << endl;
        cerr << "See help (-h option) for a list of required arguments." << endl;
        exit(1);
    }
    if (argc > 4) {
        cerr << "Too many arguments specified!" << endl;
        cerr << "See help (-h option) for usage information." << endl;
        exit(1);
    }
 
  const char* labelImgName=argv[0];
  const char* tempImgName=argv[1];
  const char* deformationfile=argv[2];

  string outputRAVENSMapPrefix=argv[3];
  string outputRAVENSMapExt;

  set<string> exts;
  exts.insert(".hdr");
  exts.insert(".img");
  exts.insert(".nii");
  exts.insert(".nii.gz");
  exts.insert(".hdr.gz");
  exts.insert(".img.gz");
  ::basis::os::path::splitext(outputRAVENSMapPrefix, outputRAVENSMapPrefix, outputRAVENSMapExt, &exts);
  if (outputRAVENSMapExt.empty()) outputRAVENSMapExt = ::basis::os::path::splitext(labelImgName, &exts)[1];
  if (outputRAVENSMapExt.empty()) outputRAVENSMapExt = ::basis::os::path::splitext(tempImgName,  &exts)[1];
  if (outputRAVENSMapExt.empty()) outputRAVENSMapExt = ".nii.gz";

  printf("\n");
  printf("input label image         = %s\n", labelImgName);
  printf("template image            = %s\n", tempImgName);
  printf("input deformation         = %s\n", deformationfile);
  printf("outupt RAVENS prefix      = %s\n", outputRAVENSMapPrefix.c_str());
  printf("outupt RAVENS extension   = %s\n", outputRAVENSMapExt.c_str());
 
  //load label images (should be in subject space)
  labelimage = ReadImage(labelImgName);
  if (labelimage == NULL) {
    cerr << "Failed to read reference image file " << labelimage << endl;
    exit(1);
  }
  if (labelimage->hdr.datatype != DT_UNSIGNED_CHAR) {
    cerr << "Input label image must be of data type DT_UNSIGNED_CHAR/DT_UINT8!" << endl;
    exit(1);
  }
  imageSizeSubj.x = labelimage->region.nx;
  imageSizeSubj.y = labelimage->region.ny;
  imageSizeSubj.z = labelimage->region.nz;
  voxelSizeSubj.x = labelimage->hdr.pixdim[1];
  voxelSizeSubj.y = labelimage->hdr.pixdim[2];
  voxelSizeSubj.z = labelimage->hdr.pixdim[3];
 
  // load template image (for image and voxel dimension in template space)
  tempimage = ReadNiftiImage(tempImgName);
  if (tempimage == NULL) {
    cerr << "Failed to read reference image file " << tempimage << endl;
    exit(1);
  }
  imageSizeTemp.x = tempimage->region.nx;
  imageSizeTemp.y = tempimage->region.ny;
  imageSizeTemp.z = tempimage->region.nz;
  voxelSizeTemp.x = tempimage->hdr.pixdim[1];
  voxelSizeTemp.y = tempimage->hdr.pixdim[2];
  voxelSizeTemp.z = tempimage->hdr.pixdim[3];

  // calculate voxel size ratio
  voxelSizeRatioS2T = (voxelSizeSubj.x * voxelSizeSubj.y * voxelSizeSubj.z)/(voxelSizeTemp.x * voxelSizeTemp.y * voxelSizeTemp.z);
  if ( (voxelSizeRatioS2T != voxelSizeRatioS2T)||(voxelSizeRatioS2T==0) )
	voxelSizeRatioS2T=1.0;
	
	
  // load deformation (should be in subject space)
  deformation = ReadDeformationField(deformationfile);
  if (deformation == NULL) {
    cerr << "Failed to read deformation field from file " << deformationfile << endl;
    delete labelimage;
    delete tempimage;
    exit(1);
  }

  // verify that the input deformation is in the same space as the input label image (both in subject space)
  if (deformation->region.nx != imageSizeSubj.x ||
        deformation->region.ny != imageSizeSubj.y ||
        deformation->region.nz != imageSizeSubj.z) {
    printf("\nError: input label image and input deformation should be in the same space (both subject space)!\n");
    printf("Note: input image must be in byte datatype.\n");
    printf("      deformation must be in DRAMMS/DRAMMS format.\n\n");
    delete labelimage;
    delete tempimage;
    delete deformation;
    exit(1);
  }
  printf("image size of RAVENS maps = (%d,%d,%d)\n\n", imageSizeTemp.x, imageSizeTemp.y, imageSizeTemp.z);
 
  DeformFld   = deformation->img.v3;
  Mass        = Falloc3d(imageSizeTemp.x, imageSizeTemp.y, imageSizeTemp.z);
  ravensimage = new Image(imageSizeTemp.x, imageSizeTemp.y, imageSizeTemp.z, DT_SIGNED_SHORT);
  ravensimage->CopyTransform(tempimage);
  ravensimage->CopyMetaData(tempimage);
  RAVENSMap   = ravensimage->img.ss;
 
  if (maskOrNot) {
    numRAVENSMaps=0;
    for (i=0;i<5;i++) {
      if (foregroundIntensity[i] != 0) numRAVENSMaps++;
    }
    maskoutImage(labelimage->img.uc, imageSizeSubj, foregroundIntensity, numRAVENSMaps);
  }

  //------------------------------
  //------------------------------
  printf("Calculate RAVENS values...\n");
  fflush(stdout);
  char filename[1024];

  for (ravensIndex=0;ravensIndex<numRAVENSMaps;ravensIndex++) {
    if (numRAVENSMaps > 1) {
      printf("------------\n");
      printf("RAVENS map for tissue type %d out of %d\n", ravensIndex+1, numRAVENSMaps);
      fflush(stdout);
    }
    volumeOriginal=0;
    volumeRAVENS=0;

    for (k=0;k<imageSizeTemp.z;k++)
      for (i=0;i<imageSizeTemp.x;i++)
        for (j=0;j<imageSizeTemp.y;j++)
          {
          Mass[k][i][j]=0.0;
          RAVENSMap[k][i][j] =0;
          }


    for (k=0;k<imageSizeSubj.z;k++)
      for (i=0;i<imageSizeSubj.x;i++)
        for (j=0;j<imageSizeSubj.y;j++)
          {
          if (labelimage->img.uc[k][i][j]==ravensIndex+1)
            {
            volumeOriginal+=1.0;
            
            for(z=-sampleNum; z<=sampleNum; z++)
              for(x=-sampleNum; x<=sampleNum; x++)
                for(y=-sampleNum; y<=sampleNum; y++)
                  {
                  Subj_subVoxel.x = x*interval + i ;
                  Subj_subVoxel.y = y*interval + j ;
                  Subj_subVoxel.z = z*interval + k ;
            
                  InterpolatedDisplacement(&Displace_subVoxel, Subj_subVoxel.x, Subj_subVoxel.y, Subj_subVoxel.z, DeformFld, imageSizeSubj.x, imageSizeSubj.y, imageSizeSubj.z) ;
                 
                  ii = Subj_subVoxel.x + Displace_subVoxel.x ;
                  jj = Subj_subVoxel.y + Displace_subVoxel.y ;
                  kk = Subj_subVoxel.z + Displace_subVoxel.z ;   // (ii,jj,kk) is in template space now
                  
                  //DistributingVolume(1.0/pow((sampleNum*2+1),3.0), ii, jj, kk, Mass, imageSizeTemp.x, imageSizeTemp.y, imageSizeTemp.z) ;
                  DistributingVolume(voxelSizeRatioS2T/pow((sampleNum*2+1),3.0), ii, jj, kk, Mass, imageSizeTemp.x, imageSizeTemp.y, imageSizeTemp.z) ;
                  }
            } // if
          } // for...for...for...

    // ***** added on Sept 16, 2011  (start)
    // avoid overflow of short value; typically the original ravens value in Mass should not be greater than 32.767,
    // otherwise it might be due to over-aggressive deformations. In practice, we have observed this happening in
    // about 30 subjects out of 606 subjects (at baseline) in accord dataset. So below is a way to avoid negative
    // RAVENS caused by data overflow.
    float numNbhVoxels;
    float upperlim=32767.0/(float)scaleFactor;
    for (k=0;k<imageSizeTemp.z;k++)
      for (i=0;i<imageSizeTemp.x;i++)
        for (j=0;j<imageSizeTemp.y;j++)
          {
          if (Mass[k][i][j]>upperlim)
              {
              massOverflow=Mass[k][i][j]-upperlim;
              Mass[k][i][j]=upperlim;   
                
              // in rear cases when overflow happens, distribute the mass to nearby 3*3*3 box
              numNbhVoxels=0.0;
              for (x=MAX(-i,-1);x<=MIN(1,(imageSizeTemp.x-1-i));x++)
                for (y=MAX(-j,-1);y<=MIN(1,(imageSizeTemp.y-1-j));y++)
                  for (z=MAX(-k,-1);z<=MIN(1,(imageSizeTemp.z-1-k));z++)
                    {
                    if (Mass[k+z][i+x][j+y]<upperlim)
                      numNbhVoxels += 1.0;    // count how many voxels in the neighborhood can be used to distribute mass to 
                    }
              if (numNbhVoxels>0.0)
              {
              massOverflow /= numNbhVoxels;
                
              for (x=MAX(-i,-1);x<=MIN(1,(imageSizeTemp.x-1-i));x++)
                for (y=MAX(-j,-1);y<=MIN(1,(imageSizeTemp.y-1-j));y++)
                  for (z=MAX(-k,-1);z<=MIN(1,(imageSizeTemp.z-1-k));z++)
                    {
                    Mass[k+z][i+x][j+y] += massOverflow;
                      
                    if (Mass[k+z][i+x][j+y]>upperlim)
                      Mass[k+z][i+x][j+y]=upperlim;
                    }
              }
              }
          }
    // ***** added on Sept 16, 2011  (end)
    

    for (k=0;k<imageSizeTemp.z;k++)
      for (i=0;i<imageSizeTemp.x;i++)
        for (j=0;j<imageSizeTemp.y;j++)
          {
          scaledMass = Mass[k][i][j]*(float)scaleFactor+0.5;
          if (scaledMass>32767) scaledMass=32767;
          if (scaledMass<0) scaledMass=0;
          RAVENSMap[k][i][j] = static_cast<short>(scaledMass);
          volumeRAVENS += Mass[k][i][j];
          }
     
    volumeOriginal *= (voxelSizeSubj.x * voxelSizeSubj.y * voxelSizeSubj.z);
    volumeRAVENS   *= (voxelSizeTemp.x * voxelSizeTemp.y * voxelSizeTemp.z);
    volumePreservePercentage = volumeRAVENS/volumeOriginal;
    printf("\n");
    printf("\tVolume of original image (in subject space) = %f\n", volumeOriginal);
    printf("\tVolume of RAVENS map (in template space)    = %f\n", volumeRAVENS);
    printf("\t%2.2f%% of the volume has been preserved in RAVENS calculation.\n", volumePreservePercentage*100);
    printf("\tMinor difference (<1%%) should be due to interpolation errors and shouldn't matter.\n");
    printf("\tMajor difference (>5%%) should be due to overflow of signed short datatype at places where ravens values are large as a result of aggressive deformations.\n");
    printf("\n");

    if (numRAVENSMaps > 1) {
      snprintf(filename, 1024, "%s_%d%s", outputRAVENSMapPrefix.c_str(), foregroundIntensity[ravensIndex], outputRAVENSMapExt.c_str());
      WriteImage(filename, ravensimage);
      printf("RAVENS map for label %d has been saved as \n   %s \nin *SIGNED SHORT* data type.\n\n", foregroundIntensity[ravensIndex], filename);
    } else {
      snprintf(filename, 1024, "%s%s", outputRAVENSMapPrefix.c_str(), outputRAVENSMapExt.c_str());
      WriteImage(filename, ravensimage);
      printf("RAVENS map has been saved as \n   %s \nin *SIGNED SHORT* data type.\n\n", filename);
    }
    fflush(stdout);
  }
 
  // clean up
  if (labelimage)   delete labelimage;
  if (tempimage)    delete tempimage;
  if (deformation)  delete deformation;
  if (ravensimage)  delete ravensimage;  
}

