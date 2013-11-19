/**
 * @file  Deform2D.cxx
 * @brief Deformabley register 2D images given a list of feature images in both image spaces.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#include <memory>
#include <common/matrix.h>  
#include <common/imageio.h>
#include <common/general.h>
#include <common/mvcd.h>
#include <common/cres.h>
#include <common/image.h>
#include <common/utilities.h>


#include <dramms/basis.h> // exename(), print_contact()


// acceptable in .cxx file
using namespace std;
using namespace dramms;


// ===========================================================================
// help
// ===========================================================================

// ---------------------------------------------------------------------------
void print_help()
{
    string exec_name = exename();
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_subj_image> <input_temp_image> <prefix_featureFile> <output_image> <def_field> " << endl;
    cout << endl;
    cout << "Description:" << endl;
    cout << "  This program outputs the registered image and the corresponding deformation" << endl;
    cout << "  field given two 2D images and files listing feature images computed from these." << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
    cout << "  <input_imageA>        Image to be registered, i.e., subject image." << endl;
    cout << "  <input_imageB>        The fixed image to be registered to, i.e., template image." << endl;
    cout << "  <prefix_featurefile>  Prefix of file that contains list of feature images." << endl;
    cout << "  <output_image>        Output registered image in space of image B."  << endl;
    cout << "  <def_field>           Output deformation field which warps image A to the space of image B." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
    cout << "  -r <int>              The number of multi-resolution pyramid levels used for optimization. (default: 3)" << endl;
    cout << "  -b <int>              The number of voxels between neigboring control points. (default: 7)" << endl;
    cout << "  -I <dir>              Request storage of intermediate images and deformations to" << endl;
    cout << "                        files in the specified directory. (default: none)" << endl;
    cout << endl;
    cout << "Standard arguments:" << endl;
    cout << "  -v                    Increase verbosity of output messages." << endl;
    cout << "  -h                    Print help and exit." << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << exec_name << " A.nii.gz B.nii.gz Gabor_ A2B.nii.gz def_A2B.nii.gz -r3 -b7" << endl;
    cout << endl;
    print_contact();
}

// ===========================================================================
// list of sub-functions (details in the end of this file)
// ===========================================================================

void FFD2D(unsigned char **sliceA, unsigned char **sliceB, unsigned char ***featureMapA, unsigned char ***featureMapB, unsigned char **sliceMask, Ivector2d sliceSize, int numFeatures, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defB2A, float **sliceA2B, int levelIndex, int numLevels);
float Bspline(float u, int level);
float calculateIntensityDifferenceSquared(unsigned char **tempSlice, float **movingSlice, Ivector2d sliceSize);
float calculateEnergyOnFeatures(unsigned char ***featureMapB, float ***featureMapA2B, Ivector2d sliceSize, int numFeatures);
void gauss2D(int numPoints1, float std1, int numPoints2, float std2, float theta, float **h);
float gauss1D(float x, float std);
float searchDisplacementAtThisControlPoint(unsigned char **sliceA, unsigned char **sliceB, float **sliceA2BTemp, float **sliceA2BIterative, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, int controlPointIndex1, int controlPointIndex2, int NthIter, int NthLevel, float initEnergy);
float searchDisplacementAtThisControlPointBasedOnFeatures(unsigned char **sliceA, unsigned char **sliceB, unsigned char ***featureMapA, unsigned char ***featureMapB, float ***featureMapA2BTemp, float ***featureMapA2BIterative, Ivector2d sliceSize, int numFeatures, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defField, int controlPointIndex1, int controlPointIndex2, int NthIter, int NthLevel, float initEnergyAB);
void GenerateFloatImageAndDeformationFieldByFFD(unsigned char **sliceA, unsigned char **sliceB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, float **sliceA2BFloat, Fvector2d **defField, int levelIndex, int numLevels);
float GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(unsigned char ***featureMapA, unsigned char ***featureMapB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, int numFeatures, float ***featureMapA2B, Fvector2d **defField, int levelIndex, int numLevels);
float FFDLocalEffect(unsigned char **sliceA, int xcoor, int ycoor, Fvector2d **controlPoints, Fvector2d **defField, Ivector2d sliceSize, int distBetweenControlPoints, int controlPointIndex1, int controlPointIndex2, float increment, char orientation, int tempOrNot, int featureIndex, int numFeatures);
void CopyFloatImage2D(float **imgA, float**imgB, Ivector2d imgSize);
void CopyFeatureMapsFloat(float ***featureMapOrig, float ***featureMapCopy, Ivector2d sliceSize, int numFeatures);
float sign(float x);
void CopyDisplacementAtControlPoints2D(Fvector2d **controlPoints, Fvector2d **controlPointsBackup, int numControlPointsX, int numControlPointsY);
void calculateDisplacementIncrementAtControlPoints(Fvector2d **controlPointsBackup, Fvector2d **controlPoints, Fvector2d **controlPointsIncrement, int numControlPointsX, int numControlPointsY);
void smoothDisplacementAtControlPoints2D(Fvector2d **controlPoints, Fvector2d **controlPointsSmoothed, int numControlPointsX, int numControlPointsY, int distBetweenControlPoints, int levelIndex, int numLevels, float factor);
void UpdateControlPointsWithSmoothIncrement2(Fvector2d **controlPointsBackup, Fvector2d **controlPointsIncrement, Fvector2d **controlPointsIncrementSmoothed, Fvector2d **controlPoints, int numControlPointsX, int numControlPointsY, float alpha, unsigned char **mask, int distBetweenControlPoints, int levelIndex, int numLevels);
void UpdateControlPointsWithSmoothIncrement(Fvector2d **controlPointsBackup, Fvector2d **controlPoints, Fvector2d **controlPointsUpdated, int numControlPointsX, int numControlPointsY, int IterIndex, int maxNumIterInResolution, int levelIndex, int numLevels, int distBetweenControlPoints, unsigned char **mask);
void LinearCombinationOfTwoControlPointsMats(Fvector2d **controlPointsNew, Fvector2d **controlPointsOld, Fvector2d **controlPointsUpdated, float weighting, int numControlPointsX, int numControlPointsY);
void Convolution2DFloat(float** input,Ivector2d inputSize,float** mask,Ivector2d maskSize,float** output);
void smoothDeformationField(Fvector2d **dfField, Fvector2d **dfFieldSmoothed, Ivector2d dfSize, int levelIndex, int numLevels);
void UpsampleDisplacementAtControlPoints(Fvector2d **controlPointsPreviousLevel, Fvector2d **controlPointsThisLevel, int numControlPointsXPreviousLevel, int numControlPointsYPreviousLevel, int numControlPointsXThisLevel, int numControlPointsYThisLevel, int distBetweenControlPoints, unsigned char **maskForControlPointsThisLevel, int levelIndex, int numLevels);
void GenerateMaskForControlPoints(unsigned char **sliceA, unsigned char **sliceB, Ivector2d sliceSize, int distBetweenControlPoints, unsigned char **mask, int levelIndex, int numLevels);
void saveDisplacementAtControlPoints2D(Fvector2d **controlPoints, int numControlPointsXThisLevel, int numControlPointsYThisLevel);
void SmoothUCImage2D(unsigned char **slice, unsigned char **sliceSmoothed, Ivector2d sliceSize);
bool ReadFeaturesFromFeatureList(const char* featureImageListFile,unsigned char ***featureMap, int x_size,int y_size);
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(unsigned char ***featureMapA, unsigned char ***featureMapB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defField, int numFeatures, float ***featureMapA2BUpdated, float ***featureMapA2BOld, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind, int controlPointIndex1, int controlPointIndex2, float increment, char orientation, int NthLevel, int updateFeaturesOrNot);
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(unsigned char ***featureMapA, unsigned char ***featureMapB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defFiled, int numFeatures, float ***featureMapA2BUpdated, float ***featureMapA2BOld, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind, int controlPointIndex1, int controlPointIndex2, float increment, char orientation, int NthLevel);
void CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(float ***featureMapOrig, float ***featureMapCopied, Ivector2d sliceSize, int numFeatures, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind);
void CopyDispAtTheBlockInducedByThisControlPoint(Fvector2d **defFieldOrig, Fvector2d **defFieldCopied, Ivector2d sliceSize, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind);
void GenerateFloatImageFromDeformationField(unsigned char **sliceA, unsigned char **sliceB, float **sliceA2BFloat, Fvector2d **defField, Ivector2d sliceSize);
float bilinearInterpolation(unsigned char **img, float x, float y, int sizeX, int sizeY);

// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc,char *argv[])
{
  int         numImgLevels                     = 3;
  int         distBetweenControlPoints         = 7;
  const char* folderNameForIntermediateResults = NULL;

  if (argc == 1) {
    print_help();
    exit(1);
  }

  int c = -1;
  while((c=getopt(argc,argv,"r:b:I:hv")) != -1)
    {
      switch(c)
        {
            case 'r':
                sscanf(optarg, "%d", &numImgLevels);
                break;
                
            case 'b':
                sscanf(optarg, "%d", &distBetweenControlPoints);
                break;
                
            case 'I':
                folderNameForIntermediateResults = optarg;
                break;

            case 'h':
                print_help();
                exit(0);

            case 'v':
                // ignored
                break;

            default:
                // error message printed by getopt() already
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (argc < 5) {
        cerr << "Missing required arguments!" << endl;
        cerr << "See help (-h option) for a list of required arguments." << endl;
        exit(1);
    }
    if (argc > 5) {
        cerr << "Too many arguments specified!" << endl;
        cerr << "See help (-h option) for usage information." << endl;
        exit(1);
    }

   const char* inputImgNameA              = argv[0];
   const char* inputImgNameB              = argv[1];
   const char* featureListPrefix          = argv[2];
   const char* outputImgNameA2B           = argv[3];
   const char* outputDeformationFieldName = argv[4];

   printf("\ninput image A = %s \ninput image B = %s\n", inputImgNameA, inputImgNameB);
   printf("\nfeature list prefix = %s\n", featureListPrefix);
   printf("\n");
   if (folderNameForIntermediateResults != NULL) {
     printf("save intermediate results, into folder %s\n", folderNameForIntermediateResults);
   }
   printf("\nmulti-resolution: number of levels = %d\n\nnumber of voxels between two adjacent control points = %d\n", numImgLevels, distBetweenControlPoints);

   struct tms startTime;
   times(&startTime);

   // read image A
   Image* imageA = ReadImage(inputImgNameA);
   auto_ptr<Image> imageA_ptr(imageA); // auto-free memory
   if (imageA == NULL) {
     fprintf(stderr, "Failed to read image A from file %s!\n", inputImgNameA);
     exit(1);
   }
   if (imageA->hdr.datatype != DT_UNSIGNED_CHAR) {
     fprintf(stderr, "Image %s has datatype %d, but DT_UNSIGNED_CHAR is required!\n", inputImgNameA, imageA->hdr.datatype);
     exit(1);
   }
   if (imageA->region.nz > 1) {
     fprintf(stderr, "Image %s has more than two dimensions! Use Deform3D to register three dimensional images.\n", inputImgNameA);
     exit(1);
   }
   imageA->SetFormat(Image::FORMAT_DRAMMS);
   // read image B
   Image* imageB = ReadImage(inputImgNameB);
   auto_ptr<Image> imageB_ptr(imageB); // auto-free memory
   if (imageB == NULL) {
     fprintf(stderr, "Failed to read image B from file %s!\n", inputImgNameB);
     exit(1);
   }
   if (imageB->hdr.datatype != DT_UNSIGNED_CHAR) {
     fprintf(stderr, "Image %s has datatype %d, but DT_UNSIGNED_CHAR is required!\n", inputImgNameB, imageB->hdr.datatype);
     exit(1);
   }
   if (imageB->region.nz > 1) {
     fprintf(stderr, "Image %s has more than two dimensions! Use Deform3D to register three dimensional images.\n", inputImgNameB);
     exit(1);
   }
   imageB->SetFormat(Image::FORMAT_DRAMMS);
   // check image size
   Ivector2d imageSize;
   imageSize.x = imageB->region.nx;
   imageSize.y = imageB->region.ny;
   if (imageA->region.nx != imageSize.x || imageA->region.ny != imageSize.y) {
     fprintf(stderr, "Size of images A and B does not match! Preregister the images using an affine registration method.\n");
     exit(1);
   }
   // multi-resolution implementation   
   Fvector2d** controlPointsPreviousLevel     = NULL;
   int         numControlPointsXPreviousLevel = 0;
   int         numControlPointsYPreviousLevel = 0;
   char        filename[1024];

   for (int i = 0; i < numImgLevels; i++) {
     // level index
     int levelIndex = numImgLevels - i;
     // image size at this resolution level
     Ivector2d imageSizeThisLevel;
     imageSizeThisLevel.x = imageSize.x / (int)pow(2, (levelIndex-1));
     imageSizeThisLevel.y = imageSize.y / (int)pow(2, (levelIndex-1));
     printf("\033[31m\n---------------------------------------------------------------\n\033[m");
     printf("\033[31m\nSearch in resolution level %d, image size = (%d, %d)\n\033[m", levelIndex, imageSizeThisLevel.x, imageSizeThisLevel.y);
     printf("\033[31m\n---------------------------------------------------------------\n\033[m");
     printf("\n");
     // downsample input image to the resolution at this level
     Image* imageAThisLevel = DownsampleImage(imageA, (int)pow(2, (levelIndex-1)));
     Image* imageBThisLevel = DownsampleImage(imageB, (int)pow(2, (levelIndex-1)));
     auto_ptr<Image> imageAThisLevel_ptr(imageAThisLevel); // auto-free memory at end of scope
     auto_ptr<Image> imageBThisLevel_ptr(imageBThisLevel); // auto-free memory at end of scope
     // allocate memory for images at this level
     Image imageA2B   (imageSizeThisLevel.x, imageSizeThisLevel.y, 1, DT_FLOAT,         1, Image::FORMAT_DRAMMS);
     Image deformation(imageSizeThisLevel.x, imageSizeThisLevel.y, 1, DT_FLOAT,         2, Image::FORMAT_DRAMMS);
     Image maskImage  (imageSizeThisLevel.x, imageSizeThisLevel.y, 1, DT_UNSIGNED_CHAR, 1, Image::FORMAT_DRAMMS);
     imageA2B   .CopyRegion(imageBThisLevel); // important for cropping as done implicitly by WriteNiftiImage()
     deformation.CopyRegion(imageBThisLevel); // currently not used as images are not padded (see Deform3D)
     maskImage  .CopyRegion(imageBThisLevel);
     imageA2B   .CopyTransform(imageBThisLevel); // defines image to be in space B
     deformation.CopyTransform(imageBThisLevel);
     maskImage  .CopyTransform(imageBThisLevel);
     // get pointers to image data
     unsigned char** imgA   = imageAThisLevel->img.uc[0];
     unsigned char** imgB   = imageBThisLevel->img.uc[0];
     float**         imgA2B = imageA2B.img.fl[0];
     Fvector2d**     defB2A = deformation.img.v2[0];
     unsigned char** mask   = maskImage.img.uc[0];
     // generate mask
     GenerateMaskForControlPoints(imgA, imgB, imageSizeThisLevel, distBetweenControlPoints, mask, levelIndex, numImgLevels);    
     // read feature images....
     printf("\nReading feature images ...\n");
     FILE* fp = NULL;
     // ...of image A
     int numFeaturesA = 0;
     sprintf(filename, "%sA_level%d.lst", featureListPrefix, levelIndex);
     printf("feature list name A = %s\n", filename);
     fp = fopen(filename, "r");
     if (fp == NULL) {
       fprintf(stderr, "Features list file %s could not be opened! Please check if file exists and its permissions.\n", filename);
       exit(1);
     }
     if (fscanf(fp, "%d", &numFeaturesA) != 1) {
       fprintf(stderr, "Failed to read number of features from file %s!\n", filename);
       fclose(fp);
       exit(1);
     }
     printf("numFeatureA = %d\n", numFeaturesA);
     Image featureMapAImage(imageSizeThisLevel.x, imageSizeThisLevel.y, numFeaturesA, DT_UNSIGNED_CHAR, 1, Image::FORMAT_DRAMMS);
     unsigned char*** featureMapA = featureMapAImage.img.uc;
     if (!ReadFeaturesFromFeatureList(filename, featureMapA, imageSizeThisLevel.x, imageSizeThisLevel.y)) {
       fprintf(stderr, "Failed to read all feature images for image A!\n");
       exit(1);
     }
     // ...of image B
     sprintf(filename, "%sB_level%d.lst", featureListPrefix, levelIndex);
     printf("feature list name B = %s\n", filename);
     int numFeaturesB = 0;
     fp = fopen(filename, "r");
     if (fp == NULL) {
       fprintf(stderr, "Features list file %s could not be opened! Please check if file exists and its permissions.\n", filename);
       exit(1);
     }
     if (fscanf(fp, "%d", &numFeaturesB) != 1) {
       fprintf(stderr, "Failed to read number of features from file %s!\n", filename);
       fclose(fp);
       exit(1);
     }
     printf("numFeatureB = %d\n", numFeaturesB);
     Image featureMapBImage(imageSizeThisLevel.x, imageSizeThisLevel.y, numFeaturesB, DT_UNSIGNED_CHAR, 1, Image::FORMAT_DRAMMS);
     unsigned char*** featureMapB = featureMapBImage.img.uc;
     if (!ReadFeaturesFromFeatureList(filename, featureMapB, imageSizeThisLevel.x, imageSizeThisLevel.y)) {
       fprintf(stderr, "Failed to read all feature images for image B!\n");
       exit(1);
     }
     // allocate memory for displacement at each control point
     int numControlPointsXThisLevel = (int)(ceil((float)imageSizeThisLevel.x/(float)distBetweenControlPoints));
     int numControlPointsYThisLevel = (int)(ceil((float)imageSizeThisLevel.y/(float)distBetweenControlPoints));
     printf(" *** number of control points at this level = %d (%d*%d)\n\n", (numControlPointsXThisLevel*numControlPointsYThisLevel), numControlPointsXThisLevel, numControlPointsYThisLevel);
     Fvector2d** controlPointsThisLevel = Fvector2dalloc2d(numControlPointsXThisLevel, numControlPointsYThisLevel);
     if (controlPointsPreviousLevel != NULL) {
       // initialize control points of this level using those of previous level
       UpsampleDisplacementAtControlPoints(controlPointsPreviousLevel, controlPointsThisLevel, numControlPointsXPreviousLevel, numControlPointsYPreviousLevel, numControlPointsXThisLevel, numControlPointsYThisLevel, distBetweenControlPoints, mask, levelIndex, numImgLevels);
       // the factor of 0.3 decreases the smooth kernel by 30%
       smoothDisplacementAtControlPoints2D(controlPointsThisLevel, controlPointsThisLevel, numControlPointsXThisLevel, numControlPointsYThisLevel, distBetweenControlPoints, levelIndex, numImgLevels, 0.3);
       // free previous control points and set pointer for next level to control points of this level
       Fvector2dfree2d(controlPointsPreviousLevel, numControlPointsXPreviousLevel);
       controlPointsPreviousLevel     = controlPointsThisLevel;
       numControlPointsXPreviousLevel = numControlPointsXThisLevel;
       numControlPointsYPreviousLevel = numControlPointsYThisLevel;
     } else {
       // initial control points of first level
       for (int x = 0; x < numControlPointsXThisLevel; x++) {
         for (int y = 0; y < numControlPointsYThisLevel; y++) {
           controlPointsThisLevel[x][y].x = 0;
           controlPointsThisLevel[x][y].y = 0;
         }
       }
     }
     // save initial transformed image and deformation field
     if (folderNameForIntermediateResults) {
       GenerateFloatImageAndDeformationFieldByFFD(imgA, imgB, imageSizeThisLevel, distBetweenControlPoints, controlPointsThisLevel, imgA2B, defB2A, levelIndex, numImgLevels);
       Image* castedA2B = CastImage(&imageA2B, DT_UNSIGNED_CHAR, false);
       sprintf(filename, "%s/A2B_level%d_init.nii.gz", folderNameForIntermediateResults, levelIndex);
       WriteImage(filename, castedA2B);
       delete castedA2B;
       printf("Saving initial deformation field in level %d\n", levelIndex);
       sprintf(filename, "%s/DField_level%d_init.nii.gz", folderNameForIntermediateResults, levelIndex);
       WriteImage(filename, &deformation);
     }
     // search displacements at control points in this level, then use the displacement at control points to generate warped image and deformation field
     FFD2D(imgA, imgB, featureMapA, featureMapB, mask, imageSizeThisLevel, numFeaturesA, distBetweenControlPoints, controlPointsThisLevel, defB2A, imgA2B, levelIndex, numImgLevels);
     // save intermediate results in this level
     if (folderNameForIntermediateResults) {
        printf("saving intermediate warped images and deformation fields into folder %s\n", folderNameForIntermediateResults);
        // image A
        sprintf(filename, "%s/A_level%d.nii.gz", folderNameForIntermediateResults, levelIndex);
        WriteImage(filename, imageAThisLevel);
        // image B
        sprintf(filename, "%s/B_level%d.nii.gz", folderNameForIntermediateResults, levelIndex);
        WriteImage(filename, imageBThisLevel);
        // deformed image
        Image* castedA2B = CastImage(&imageA2B, DT_UNSIGNED_CHAR, false);
        sprintf(filename, "%s/A2B_level%d.nii.gz", folderNameForIntermediateResults, levelIndex);
        WriteImage(filename, castedA2B);
        delete castedA2B;
        // deformation field
        sprintf(filename, "%s/DField_level%d.nii.gz", folderNameForIntermediateResults, levelIndex);
        WriteImage(filename, &deformation);
      }
      // save final result at finest resolution
      if (i == (numImgLevels - 1)) {
        printf("\n\nAt the finest resolution,  save final results: \n");
        Image* castedA2B = CastImage(&imageA2B, DT_UNSIGNED_CHAR, false);
        WriteImage(outputImgNameA2B, castedA2B);
        delete castedA2B;
        printf("... %s ", outputImgNameA2B);
        WriteImage(outputDeformationFieldName, &deformation);
        printf(" ... %s !\n\n", outputDeformationFieldName);
      }
   }
   if (controlPointsPreviousLevel) Fvector2dfree2d(controlPointsPreviousLevel, numControlPointsXPreviousLevel);
   // print total running time
   struct tms endTime;
   times(&endTime);
   double durTime=((double)endTime.tms_utime-(double)startTime.tms_utime);
   printf("\n\nDuration=%.2f seconds\n",durTime/100.0); 
}

// ===========================================================================
// auxiliary functions
// ===========================================================================

// sub-function 7
void FFD2D(unsigned char **sliceA, unsigned char **sliceB, unsigned char ***featureMapA, unsigned char ***featureMapB, unsigned char **sliceMask, Ivector2d sliceSize, int numFeatures, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defB2A, float **sliceA2B, int levelIndex, int numLevels)
{
  int maxNumIterInResolution=5;  // maxIter1
  int i;
  int IterIndex, controlPointIndex1, controlPointIndex2, featureIndex;
  Ivector2d realCoorThisControlPoint;
  int numControlPointsX = (int)(ceil((float)sliceSize.x/distBetweenControlPoints));
  int numControlPointsY = (int)(ceil((float)sliceSize.y/distBetweenControlPoints));
  Fvector2d **controlPointsBackup;
  controlPointsBackup = Fvector2dalloc2d(numControlPointsX, numControlPointsY);
      
  float energyAB[maxNumIterInResolution+1];
  for (i=0;i<=maxNumIterInResolution;i++)
    energyAB[i]=0;
    
  // get initial enerygy
  float ***featureMapA2BTemp, ***featureMapA2BIterative;
  featureMapA2BTemp = (float***)malloc(sizeof(float**)*numFeatures);       
  featureMapA2BIterative = (float***)malloc(sizeof(float**)*numFeatures);
  for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
    {
     featureMapA2BTemp[featureIndex]=Falloc2d(sliceSize.x, sliceSize.y); 
     featureMapA2BIterative[featureIndex]=Falloc2d(sliceSize.x, sliceSize.y); 
    }
  energyAB[0] = GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BTemp, defB2A, levelIndex, numLevels);
  CopyFeatureMapsFloat(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures);
 
  float diffRatio, energyABIterative, energyABTemp;
  energyABIterative = energyAB[0];
  energyABTemp = energyAB[0];
  printf("energy = %f\n", energyAB[0]);
 
  for (IterIndex=0;IterIndex<maxNumIterInResolution;IterIndex++)
    {
      CopyDisplacementAtControlPoints2D(controlPoints, controlPointsBackup, numControlPointsX, numControlPointsY);  // copy: controlPoints -> controlPointsBackup
      
      
      //----------------------------
      // Calculate displacement at all control points
      //----------------------------
      for (controlPointIndex1=0;controlPointIndex1<numControlPointsX;controlPointIndex1++)
        for (controlPointIndex2=0;controlPointIndex2<numControlPointsY;controlPointIndex2++)
          {
          // check if this control point falls into the mask
          realCoorThisControlPoint.x = controlPointIndex1*distBetweenControlPoints;
          realCoorThisControlPoint.y = controlPointIndex2*distBetweenControlPoints;
          
          if (sliceMask[realCoorThisControlPoint.x][realCoorThisControlPoint.y]>0)
            {
            printf("\033[33m\n\nUpdate control point (%d, %d) => the real coordinate (%d, %d)\n\033[m", controlPointIndex1, controlPointIndex2, realCoorThisControlPoint.x, realCoorThisControlPoint.y);  //no display
            energyABIterative = searchDisplacementAtThisControlPointBasedOnFeatures(sliceA, sliceB, featureMapA, featureMapB, featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, distBetweenControlPoints, controlPoints, defB2A, controlPointIndex1, controlPointIndex2, IterIndex, levelIndex, energyABTemp); // "controlPoints" will be updated at this control point
            energyABTemp = energyABIterative;
            } // if
          else
            {
            controlPoints[controlPointIndex1][controlPointIndex2].x=0;
            controlPoints[controlPointIndex1][controlPointIndex2].y=0;
            } // else
          } // search displacement at each control point
      
      //----------------------------------------------------------
      // Check: after one iteration of search displacement at all control points, check the total energy. If it is relatively stable, no need for further iterations in this resolution.
      //----------------------------------------------------------
      UpdateControlPointsWithSmoothIncrement(controlPointsBackup, controlPoints, controlPoints, numControlPointsX, numControlPointsY, IterIndex, maxNumIterInResolution, levelIndex, numLevels, distBetweenControlPoints, sliceMask);
      smoothDisplacementAtControlPoints2D(controlPoints, controlPoints, numControlPointsX, numControlPointsY, distBetweenControlPoints, levelIndex, numLevels, 0.04*pow(1.5,(float)(3-levelIndex)));  // this factor decreases the smooth kernel to its 4% 
      energyAB[IterIndex+1] = GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative, defB2A, levelIndex, numLevels);
      diffRatio = (energyAB[IterIndex+1]-energyAB[IterIndex])/energyAB[IterIndex];
      energyABTemp = energyAB[IterIndex+1];
      
      printf("\033[35m\n** %d-th iteration in this level\n\033[m", IterIndex+1);;
      printf("\033[35m** energy: \033[m");
      for (i=0;i<=maxNumIterInResolution;i++)
        printf("\033[35m %f  \033[m", energyAB[i]);
      printf("\033[35m\n**\n\n\n\033[m");
      
      // when total energy doesn't decrease too much, it becomes stable, then terminate iterations at this resolution/level, this requirement gets stricker as resolutions gets finer
      printf("\nEnergy reduced by %f%%\n\n", -100*diffRatio);
      if ( (fabs(diffRatio)<0.075*exp(-pow((numLevels-levelIndex),2.0)/((numLevels+1)*2.0)))  || (IterIndex==(maxNumIterInResolution-1)))
        {
        GenerateFloatImageFromDeformationField(sliceA, sliceB, sliceA2B, defB2A, sliceSize);
        break;
        }
        
      //when total energy increases too much, go back to the previous and terminate iteration at this resolution/level; This requirement gets stricker as resolution gets finer
      else if ( diffRatio>0.01*exp(-pow((numLevels-levelIndex),2.0)/(float)numLevels) )
        {
        LinearCombinationOfTwoControlPointsMats(controlPoints, controlPointsBackup, controlPoints, 0.85, numControlPointsX, numControlPointsY);
        GenerateFloatImageFromDeformationField(sliceA, sliceB, sliceA2B, defB2A, sliceSize);
        break;
        }
        
    } // for (IterIndex=0;IterIndex<maxNumIterInResolution;IterIndex++)
    
  Fvector2dfree2d(controlPointsBackup, numControlPointsX);
  free(featureMapA2BTemp);
  free(featureMapA2BIterative);
}


// sub-function 8
float Bspline(float u, int level)
{
  float value;
  switch (level)
    {
     case 0:
        value = pow((1-u),3)/6;
        break;
     case 1:
        value = (3*pow(u,3)-6*pow(u,2)+4)/6;
        break;
     case 2:
        value = (-3*pow(u,3)+3*pow(u,2)+3*u+1)/6;
        break;
     case 3:
        value = pow(u,3)/6;
        break;
     default:
        value = 0.0;
        break;
    }
  
  return value;
}


// sub-function 10
float calculateIntensityDifferenceSquared(unsigned char **tempSlice, float **movingSlice, Ivector2d sliceSize)
{
  int x,y;
  float energy=0.0;
  
  for (x=0;x<sliceSize.x;x++)
    for (y=0;y<sliceSize.y;y++)
      {
        energy += pow( ((float)tempSlice[x][y]-movingSlice[x][y]), 2.0);
      }
      
  return energy;
}


//sub-function 10-f
float calculateEnergyOnFeatures(unsigned char ***featureMapB, float ***featureMapA2B, Ivector2d sliceSize, int numFeatures)
{
  int x,y;
  int featureIndex;
  float totalEnergy=0.0;
  float energySquaredAtThisPoint;
  
  for (x=0;x<sliceSize.x;x++)
    for (y=0;y<sliceSize.y;y++)
      {
      energySquaredAtThisPoint = 0;
      for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
        energySquaredAtThisPoint += pow( ((float)featureMapB[featureIndex][x][y] - featureMapA2B[featureIndex][x][y]), 2.0 );
        
      totalEnergy += energySquaredAtThisPoint;
      }
      
  totalEnergy /= numFeatures;   // Ou added on 10/09/2008
  
  return totalEnergy;
}

// sub-function 11
void gauss2D(int numPoints1, float std1, int numPoints2, float std2, float theta, float **h)
{
  float rotation[2][2];
  float u[2];
  int i,j;
  float sumH=0.0;
  
  rotation[0][0]=cos(theta);
  rotation[0][1]=-sin(theta);
  rotation[1][0]=sin(theta);
  rotation[1][1]=cos(theta);

  for (i=0;i<numPoints2;i++)
    for (j=0;j<numPoints1;j++)
      {
        u[0]=rotation[0][0]*(j-(numPoints1-1)/2)+rotation[0][1]*(i-(numPoints2-1)/2);
        u[1]=rotation[1][0]*(j-(numPoints1-1)/2)+rotation[1][1]*(i-(numPoints2-1)/2);
        h[i][j]=gauss1D(u[0],std1)*gauss1D(u[1],std2);
        //printf("u = (%f, %f), h[%d][%d] = %f\n", u[0], u[1],i,j,h[i][j]);
        sumH += h[i][j];
      }
        
  for (i=0;i<numPoints2;i++)
    for (j=0;j<numPoints1;j++)
      {
        h[i][j]/= sumH;
      }
}


// sub-function 12
float gauss1D(float x, float std)
{
  float value;
  value = exp(-pow(x,2.0)/(2.0*pow(std,2.0))) / (std*sqrt(2*G_PI));
  
  return value;
}


//sub-function 13-f
float searchDisplacementAtThisControlPointBasedOnFeatures(unsigned char **sliceA, unsigned char **sliceB, unsigned char ***featureMapA, unsigned char ***featureMapB, float ***featureMapA2BTemp, float ***featureMapA2BIterative, Ivector2d sliceSize, int numFeatures, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defField, int controlPointIndex1, int controlPointIndex2, int NthIter, int NthLevel, float initEnergyAB)
{
  int maxIter1=2;
  //if (NthIter==0)  maxIter1=2;
  int maxIter2=1;
  int iterIndex1, iterIndex2;
  int featureIndex;
  int i,j;
  int xcoor, ycoor;
  //int xcoor, ycoor;
  float rate;
    
  float maxDisplacementAllowed;
  if (distBetweenControlPoints<14)
    maxDisplacementAllowed = distBetweenControlPoints * pow(1.35, 4-NthLevel);
  else 
    maxDisplacementAllowed = distBetweenControlPoints * pow(1.25, 3-NthLevel);
    
  
  float deltaX = 0.3 / exp(-pow((float)NthIter,2.0)/25.0) * (float)distBetweenControlPoints/3.0;
  float deltaY = 0.3 / exp(-pow((float)NthIter,2.0)/25.0) * (float)distBetweenControlPoints/3.0;
  
  float minIncrementX = 0.5*pow(1.25,(3-NthLevel))*exp((float)distBetweenControlPoints/50.0);
  float minIncrementY = 0.5*pow(1.25,(3-NthLevel))*exp((float)distBetweenControlPoints/50.0);
  float maxIncrementX = (float)distBetweenControlPoints*0.5*pow(1.1, 3-NthLevel);
  float maxIncrementY = (float)distBetweenControlPoints*0.5*pow(1.1, 3-NthLevel);
  
  printf("\n(deltax, deltay) = (%f, %f)\n",deltaX, deltaY);  //no display
  printf("minimum increment required in each step = (%f, %f)\n", minIncrementX, minIncrementY); // no display
  printf("maximum increment allowed in each step = (%f, %f)\n", maxIncrementX, maxIncrementY); // no display
  
  Ivector2d realCoorThisControlPoint;
  realCoorThisControlPoint.x = controlPointIndex1*distBetweenControlPoints;
  realCoorThisControlPoint.y = controlPointIndex2*distBetweenControlPoints;
  
  float dxIterative=controlPoints[controlPointIndex1][controlPointIndex2].x;
  float dyIterative=controlPoints[controlPointIndex1][controlPointIndex2].y;
  float diffEnergyX, diffEnergyY, deltaEnergyXSquared, deltaEnergyYSquared;
  float gradientX, gradientY;
  float incrementX, incrementY;
  float actualIncrementX, actualIncrementY;
  float previousIncrementX=0, previousIncrementY=0;
  float eta;
  
  int affectingRadiusAhead = distBetweenControlPoints*2;
  int affectingRadiusBehind = distBetweenControlPoints*2;
  int numControlPointsX = (int)(ceil((float)sliceSize.x/distBetweenControlPoints));
  int numControlPointsY = (int)(ceil((float)sliceSize.y/distBetweenControlPoints));
  
  Fvector2d** defFieldTemp = Fvector2dalloc2d(sliceSize.x, sliceSize.y);
  for (i=-affectingRadiusAhead;i<affectingRadiusBehind;i++)
    for (j=-affectingRadiusAhead;j<affectingRadiusBehind;j++)
      {
      xcoor = realCoorThisControlPoint.x + i;
      ycoor = realCoorThisControlPoint.y + j;
      if ( (xcoor>=0)&(xcoor<sliceSize.x)&(ycoor>=0)&(ycoor<sliceSize.y) )
        {
        defFieldTemp[xcoor][ycoor].x = defField[xcoor][ycoor].x;
        defFieldTemp[xcoor][ycoor].y = defField[xcoor][ycoor].y;
        }
      }
  
  
  //float energyAB = calculateIntensityDifferenceSquared(sliceB, sliceA2BIterative, sliceSize);
  float energyAB = initEnergyAB;
  float previousEnergyAB = energyAB;
  float ENERGY[maxIter1+1];
  
    
    
  for (iterIndex1=0;iterIndex1<maxIter1;iterIndex1++)
    {
    ENERGY[iterIndex1]=energyAB;
    previousEnergyAB = energyAB;
    
    //------------------------------------
    // search in x direction
    //-----------------------------------
    for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
      {
        printf("X: Iter #%d, %d\n", iterIndex1, iterIndex2); // no display
        controlPoints[controlPointIndex1][controlPointIndex2].x += deltaX;
        
        diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, deltaX, 'X', NthLevel);
        // gradient
        gradientX = diffEnergyX/deltaX;
        printf("\tdiffEnergyX/energyAB = %f, gradientX = %f\n", diffEnergyX/energyAB, gradientX);  // no display
        
        // if change is too small, no need to update displacement at this control point, break;
        if ( fabs(diffEnergyX/energyAB) <= 0.000002/ ((float)distBetweenControlPoints/3)* exp(-pow(NthIter,2.0)/5.0) )
          {
           controlPoints[controlPointIndex1][controlPointIndex2].x = dxIterative;
           printf("\tGradientX is too small!\n"); // no display
           printf("\tdx = %f, dy = %f, diffEnergyX = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyX); // no display
           break;
          } // if gradient is too small, break;
        else // else, proceed along the gradient direction
          {
           // get a proper eta
           //eta = 850*exp(-pow(iterIndex1,2.0)/100)*exp(-pow(iterIndex2,2.0)/100)* 1/energyAB * pow(12, (3-NthLevel));
           eta = 500*exp(-pow(iterIndex1,2.0)/100.0)*exp(-pow(iterIndex2,2.0)/100.0)* pow(10, -8)*(energyAB/(sliceSize.x*sliceSize.y*distBetweenControlPoints*distBetweenControlPoints)) * pow(2.0, (3-NthLevel));
           
           if (diffEnergyX/energyAB > 0.5)   
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*0.1*diffEnergyX/energyAB;   // in case of big gradient, do not behave too aggressively
           else if ( (diffEnergyX/energyAB>0.2) & (diffEnergyX/energyAB<=0.5) )
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*0.5*diffEnergyX/energyAB;
           else if ( (diffEnergyX/energyAB >0.1) & (diffEnergyX/energyAB<=0.2) )
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*diffEnergyX/energyAB;
           
           // get a proper increment 
           incrementX = - eta*gradientX;
                   
           if (fabs(incrementX)<minIncrementX)
             incrementX = sign(incrementX)*minIncrementX;
           else if (fabs(incrementX)>maxIncrementX)
             incrementX = sign(incrementX)*maxIncrementX;  // do not move too much in one iteration
           printf("\teta = %f, gradientX = %f\n\ttentative incrementX = %f", eta, gradientX, incrementX); // no display
           
           // contraint 1 on increment : it can not exceed maximum displacement allowed
           dxIterative = dxIterative + incrementX;
           if (fabs(dxIterative)>maxDisplacementAllowed)
             dxIterative = sign(dxIterative)*maxDisplacementAllowed;
             
           // constraint 2 on increment: it must guarantee that dx proceeds in the right direction; otherwise, increment=0.1*increment
           actualIncrementX = dxIterative - controlPoints[controlPointIndex1][controlPointIndex2].x + deltaX;
           printf(", actual incrementX = %f\n", actualIncrementX);
           controlPoints[controlPointIndex1][controlPointIndex2].x = dxIterative;
            
           if ( ((incrementX==-previousIncrementX)&(fabs(incrementX)==maxIncrementX)) || ((incrementX==-previousIncrementX)&(fabs(incrementX)==minIncrementX)) )
            {// abondon this move
            printf("\tAbondon this move\n"); // no display
            controlPoints[controlPointIndex1][controlPointIndex2].x -= incrementX;
            dxIterative -= incrementX;
            printf("\tdx = %f, dy = %f\n\tenergy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, previousEnergyAB); // no display
            break;
            }
            
           // calculate energy change brought about by incrementX
           CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
           diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementX, 'X', NthLevel, NNO);
           energyAB = previousEnergyAB + diffEnergyX;
           rate = diffEnergyX/previousEnergyAB;
           printf("\tenergyAB = %f, \n\tpreviousEnergyAB = %f, \n\tchangeRate = %f\n", energyAB, previousEnergyAB, rate );   // no display
           
           
           
           if ( (fabs(rate)<0.00005*pow(1.5,(NthLevel-3))) || ((fabs(dxIterative)==maxDisplacementAllowed)&&(rate<0.01)) )
            {
            // too little change to be continued or the displacement has reached the maximum allowed
            printf("\tToo little change or the displacement has reached maximum allowed -- no need to move in this direction at this stage!\n"); // no display
            printf("\tdx = %f, dy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y); // no display
            previousEnergyAB = energyAB;
            previousIncrementX = 0;
            CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            break;
            }
           else if ( rate > 0.01 * pow(0.8, (3-NthLevel)) )
            {
            // moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
            dxIterative = dxIterative - incrementX + 0.2*incrementX*0.01* pow(0.8, (3-NthLevel))/rate;
            actualIncrementX = 0.2*incrementX*0.01* pow(0.8, (3-NthLevel))/rate;
            controlPoints[controlPointIndex1][controlPointIndex2].x = dxIterative;
            //GenerateFloatImageByFFD(sliceA, sliceB, sliceSize, distBetweenControlPoints, controlPoints, sliceA2BIterative);
            //GenerateFeaturesByFFD(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative);
            CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementX, 'X', NthLevel, NNO);
            energyAB = previousEnergyAB + diffEnergyX;
            //energyAB = calculateIntensityDifferenceSquared(sliceB, sliceA2BIterative, sliceSize);
            //energyAB = calculateEnergyOnFeatures(featureMapB, featureMapA2BIterative, sliceSize, numFeatures);
            printf("\thas moved too much.. Aho, move back!\n"); // no display
            printf("\tdx = %f, dy = %f, diffEnergyX = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyX, energyAB); // no display
            if (energyAB>1.01*previousEnergyAB)
              {
              dxIterative = dxIterative - 0.2*incrementX*0.01* pow(0.8, (3-NthLevel))/rate;
              actualIncrementX = 0;
              controlPoints[controlPointIndex1][controlPointIndex2].x = dxIterative;
              //GenerateFeaturesByFFD(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative);
              printf("\thas moved too much.. Aho, move back!\n"); // no display
              printf("\tdx = %f, dy = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, previousEnergyAB); // no display
              diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementX, 'X', NthLevel, YYES);
              break;
              }
            else if (energyAB>0.9995*previousEnergyAB && energyAB<=1.01*previousEnergyAB)
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              break;            
              }
            else
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              }
            }
           else if ( (rate > 0.000025 * pow(0.8, (3-NthLevel))) && (rate <= 0.01 * pow(0.8, (3-NthLevel))) )
            {
            // moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
            dxIterative = dxIterative - 0.5*incrementX;
            actualIncrementX = 0.5*incrementX;
            controlPoints[controlPointIndex1][controlPointIndex2].x = dxIterative;
            //GenerateFloatImageByFFD(sliceA, sliceB, sliceSize, distBetweenControlPoints, controlPoints, sliceA2BIterative);
            //GenerateFeaturesByFFD(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative);
            //energyAB = calculateIntensityDifferenceSquared(sliceB, sliceA2BIterative, sliceSize);
            //energyAB = calculateEnergyOnFeatures(featureMapB, featureMapA2BIterative, sliceSize, numFeatures);
            CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementX, 'X', NthLevel, NNO);
            energyAB = previousEnergyAB + diffEnergyX;
            printf("\thas moved too much.. Aho, move back!\n"); // no display
            printf("\tdx = %f, dy = %f, diffEnergyX = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyX, energyAB); // no display
            if (energyAB>1.01*previousEnergyAB)
              {
              dxIterative = dxIterative - 0.5*incrementX;
              actualIncrementX = 0;
              controlPoints[controlPointIndex1][controlPointIndex2].x = dxIterative;
              //GenerateFeaturesByFFD(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative);
              printf("\thas moved too much.. Aho, move back!\n"); // no display
              printf("\tdx = %f, dy = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, previousEnergyAB); // no display
              diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementX, 'X', NthLevel, YYES);
              break;
              }
            else if (energyAB>0.9995*previousEnergyAB && energyAB<=1.01*previousEnergyAB)
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              break;            
              }
            else
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              break;
              }
            }
           else
            {
            previousEnergyAB = energyAB; 
            previousIncrementX = incrementX;
            //CopyFloatImage2D(sliceA2BTemp, sliceA2BIterative, sliceSize);
            //CopyFeatureMapsFloat(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures);
            CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            } // else, increment is proper
           printf("\tdx = %f, dy = %f, diffEnergyX = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyX); // no display
          } // else, proceed along the gradient direction
      } //for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
      
    
    //-----------------------------------
    // search in y direction
    //-----------------------------------
    for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
      {
        printf("Y: Iter #%d, %d\n", iterIndex1, iterIndex2); // no display
        controlPoints[controlPointIndex1][controlPointIndex2].y += deltaY;
        
        /* // commented on July 10
        diffEnergyY=0;
  
        // calculate energy change brought about by detlaY
        for (i=-affectingRadiusAhead;i<affectingRadiusBehind;i++)
          for (j=-affectingRadiusAhead;j<affectingRadiusBehind;j++)
            {
            xcoor = realCoorThisControlPoint.x + j;
            ycoor = realCoorThisControlPoint.y + i;
        
            if (xcoor>=0 && xcoor<sliceSize.x && ycoor>=0 && ycoor<sliceSize.y)
              {
              deltaEnergyYSquared = 0.0;
              for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
                {
                featureMapA2BTemp[featureIndex][xcoor][ycoor] = FFDLocalEffect(featureMapA[featureIndex], xcoor, ycoor, controlPoints, sliceSize, distBetweenControlPoints);
                deltaEnergyYSquared += (featureMapA2BTemp[featureIndex][xcoor][ycoor]+featureMapA2BIterative[featureIndex][xcoor][ycoor]-2*(float)featureMapB[featureIndex][xcoor][ycoor])*(featureMapA2BTemp[featureIndex][xcoor][ycoor]-featureMapA2BIterative[featureIndex][xcoor][ycoor]);
                }
              diffEnergyY += deltaEnergyYSquared; 
              
              }
            } // calculate energy change brought by deltaX
        */ // commented on July 10
        
        diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, deltaY, 'Y', NthLevel);
        
        // gradient
        gradientY = diffEnergyY/deltaY;
        printf("\tdiffEnergyX/energyAB = %f, gradientY = %f\n", diffEnergyY/energyAB, gradientY);  // no display
        
        if ( fabs(diffEnergyY/energyAB) <= 0.000002/ ((float)distBetweenControlPoints/3)* exp(-pow(NthIter,2.0)/5.0) )
          {
           controlPoints[controlPointIndex1][controlPointIndex2].y = dyIterative;
           printf("\tGradientY is too small!\n"); // no display
           printf("\tdx = %f, dy = %f, diffEnergyY = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyY); // no display
           break;
          } // if gradient is too small, break;
        else // else, proceed along the gradient direction
          {
           // get a proper eta
           //eta = 850*exp(-pow(iterIndex1,2.0)/100)*exp(-pow(iterIndex2,2.0)/100)* 1/energyAB * pow(12, (3-NthLevel)) ;
           eta = 500*exp(-pow(iterIndex1,2.0)/100.0)*exp(-pow(iterIndex2,2.0)/100.0)* pow(10, -8)*(energyAB/(sliceSize.x*sliceSize.y*distBetweenControlPoints*distBetweenControlPoints)) * pow(2.0, (3-NthLevel));
           
           if (diffEnergyY/energyAB > 0.5)   
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*0.1*diffEnergyY/energyAB;   // in case of big gradient, do not behave too aggressively
           else if ( (diffEnergyY/energyAB>0.2) & (diffEnergyY/energyAB<=0.5) )
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*0.5*diffEnergyY/energyAB;
           else if ( (diffEnergyY/energyAB > 0.1) & (diffEnergyY/energyAB<=0.2) )
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*diffEnergyY/energyAB;
           
           // get a proper increment 
           incrementY = - eta*gradientY;
           
           if (fabs(incrementY)<minIncrementY)
             incrementY = sign(incrementY)*minIncrementY;
           else if (fabs(incrementY)>maxIncrementY)
             incrementY = sign(incrementY)*maxIncrementY;  // do not move too little or too much in one iteration
           printf("\teta = %f, gradientY = %f\n\ttentative incrementY = %f", eta, gradientY, incrementY); // no display
           
           
           // contraint 1 on increment : it can not exceed 90% of the distance between two adjacent control points
           dyIterative = dyIterative + incrementY;
           if (fabs(dyIterative)>maxDisplacementAllowed)
             dyIterative = sign(dyIterative)*maxDisplacementAllowed;
             
           
           // constraint 2 on increment: it must guarantee that dx proceeds in the right direction; otherwise, increment=0.1*increment
           actualIncrementY = dyIterative - controlPoints[controlPointIndex1][controlPointIndex2].y + deltaY;
           printf(", actual incrementY = %f\n", actualIncrementY);
           controlPoints[controlPointIndex1][controlPointIndex2].y = dyIterative;
           
           if ( ((incrementY==-previousIncrementY)&(fabs(incrementY)==maxIncrementY)) || ((incrementY==-previousIncrementY)&(fabs(incrementY)==minIncrementY)) )
            {// abondon this move
            printf("\tAbondon this move\n"); // no display
            controlPoints[controlPointIndex1][controlPointIndex2].y -= incrementY;
            dyIterative -= incrementY;
            printf("\tdx = %f, dy = %f\n\tenergy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, previousEnergyAB); // no display
            break;
            }
            
           // calculate energy change brought about by incrementY
           CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
           diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementY, 'Y', NthLevel, NNO);
               
           //GenerateFloatImageByFFD(sliceA, sliceB, sliceSize, distBetweenControlPoints, controlPoints, sliceA2BTemp);
           //energyAB = calculateIntensityDifferenceSquared(sliceB, sliceA2BTemp, sliceSize);
           energyAB = previousEnergyAB + diffEnergyY;
           rate = diffEnergyY/previousEnergyAB;
           printf("\tenergyAB = %f, \n\tpreviousEnergyAB = %f, \n\tchangeRate = %f\n", energyAB, previousEnergyAB, rate );   // no display
           
           
           
           if ( (fabs(rate)<0.00005*pow(1.5,(NthLevel-3))) || ((fabs(dyIterative)==maxDisplacementAllowed)&&(rate<0.01)) )
            {
            // too little change to be continued or the displacement has reached the maximum allowed
            printf("\tToo little change or the displacement has reached maximum allowed -- no need to move in this direction at this stage!\n"); // no display
            printf("\tdx = %f, dy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y); // no display
            previousEnergyAB = energyAB;
            previousIncrementY = 0;
            CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            break;
            }
           else if ( rate > 0.01 * pow(0.8, (3-NthLevel)) )
            {
            // moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
            dyIterative = dyIterative - incrementY + 0.2*incrementY*0.01* pow(0.8, (3-NthLevel))/rate;
            actualIncrementY = 0.2*incrementY*0.01* pow(0.8, (3-NthLevel))/rate;
            controlPoints[controlPointIndex1][controlPointIndex2].y = dyIterative;
            //GenerateFloatImageByFFD(sliceA, sliceB, sliceSize, distBetweenControlPoints, controlPoints, sliceA2BIterative);
            //GenerateFeaturesByFFD(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative);
            //energyAB = calculateIntensityDifferenceSquared(sliceB, sliceA2BIterative, sliceSize);
            //energyAB = calculateEnergyOnFeatures(featureMapB, featureMapA2BIterative, sliceSize, numFeatures);
            CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementY, 'Y', NthLevel, NNO);
            energyAB = previousEnergyAB + diffEnergyY;
            printf("\thas moved too much.. Aho, move back!\n"); // no display
            printf("\tdx = %f, dy = %f, diffEnergyY = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyX, energyAB); // no display
            if (energyAB>1.01*previousEnergyAB)
              {
              dyIterative = dyIterative - 0.2*incrementY*0.01* pow(0.8, (3-NthLevel))/rate;
              actualIncrementY = 0;
              controlPoints[controlPointIndex1][controlPointIndex2].y = dyIterative;
              printf("\thas moved too much.. Aho, move back!\n"); // no display
              printf("\tdx = %f, dy = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, previousEnergyAB); // no display
              diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementY, 'Y', NthLevel, YYES);
              break;
              }
            else if (energyAB>0.9995*previousEnergyAB && energyAB<=1.01*previousEnergyAB)
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              break;            
              }
            else
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              }
            }
           else if ( (rate > 0.000025 * pow(0.8, (3-NthLevel))) && (rate <= 0.01 * pow(0.8, (3-NthLevel))) )
            {
            // moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
            dyIterative = dyIterative - 0.5*incrementY;
            actualIncrementY = 0.5*incrementY;
            controlPoints[controlPointIndex1][controlPointIndex2].y = dyIterative;
            //GenerateFloatImageByFFD(sliceA, sliceB, sliceSize, distBetweenControlPoints, controlPoints, sliceA2BIterative);
            //GenerateFeaturesByFFD(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, numFeatures, featureMapA2BIterative);
            //energyAB = calculateIntensityDifferenceSquared(sliceB, sliceA2BIterative, sliceSize);
            //energyAB = calculateEnergyOnFeatures(featureMapB, featureMapA2BIterative, sliceSize, numFeatures);
            CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementY, 'Y', NthLevel, NNO);
            energyAB = previousEnergyAB + diffEnergyY;
            printf("\thas moved too much.. Aho, move back!\n"); // no display
            printf("\tdx = %f, dy = %f, diffEnergyX = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyX, energyAB); // no display
            if (energyAB>1.01*previousEnergyAB)
              {
              dyIterative = dyIterative - 0.5*incrementY;
              actualIncrementY = 0;
              controlPoints[controlPointIndex1][controlPointIndex2].y = dyIterative;
              printf("\thas moved too much.. Aho, move back!\n"); // no display
              printf("\tdx = %f, dy = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, previousEnergyAB); // no display
              diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, sliceSize, distBetweenControlPoints, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind, controlPointIndex1, controlPointIndex2, actualIncrementY, 'Y', NthLevel, YYES);
              break;
              }
            else if (energyAB>0.9995*previousEnergyAB && energyAB<=1.01*previousEnergyAB)
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              break;            
              }
            else
              {
              previousEnergyAB = energyAB;
              CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
              break;
              }
            }
           else
            {
            previousEnergyAB = energyAB; 
            previousIncrementY = incrementY;
            //CopyFloatImage2D(sliceA2BTemp, sliceA2BIterative, sliceSize);
            //CopyFeatureMapsFloat(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures);
            CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, sliceSize, numFeatures, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, sliceSize, realCoorThisControlPoint, affectingRadiusAhead, affectingRadiusBehind);
            } // else, increment is proper
           printf("\tdx = %f, dy = %f, diffEnergyY = %f\n\n", controlPoints[controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex1][controlPointIndex2].y, diffEnergyY); // no display
          } // else, proceed along the gradient direction
      } //for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
    
    //----------------------------------
    // Check if continue searching at this control point
    //----------------------------------
    energyAB = previousEnergyAB;
    ENERGY[iterIndex1+1]=energyAB;
    if ( (fabs((ENERGY[iterIndex1+1]-ENERGY[iterIndex1])/ENERGY[iterIndex1])<0.001*pow(0.8,(3-NthLevel))) | (ENERGY[iterIndex1+1]-ENERGY[iterIndex1]>0) )
        break;
    } //for (iterIndex1=0;iterIndex1<maxIter1;iterIndex1++)

  // release memory
  Fvector2dfree2d(defFieldTemp, sliceSize.x);
  
  return energyAB;

}

// sub-function 14
void GenerateFloatImageAndDeformationFieldByFFD(unsigned char **sliceA, unsigned char **sliceB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, float **sliceA2BFloat, Fvector2d **defField, int levelIndex, int numLevels)
{
  int x,y;
  int indexEffectiveControlPointX, indexEffectiveControlPointY;
  int m,n;
  float dx, dy;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY;
  float u,v;
  int numControlPointsX = (int)(ceil((float)sliceSize.x/distBetweenControlPoints));
  int numControlPointsY = (int)(ceil((float)sliceSize.y/distBetweenControlPoints));
  
  for (x=0;x<sliceSize.x;x++)
    for (y=0;y<sliceSize.y;y++)
      { //for each point
        dx=0;
        dy=0;
        
        // find out effecting control points
        indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPoints)-1;
        indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPoints)-1;
        u = (float)x/(float)distBetweenControlPoints - floor((float)x/(float)distBetweenControlPoints);
        v = (float)y/(float)distBetweenControlPoints - floor((float)y/(float)distBetweenControlPoints);
        
        // calculate displacement at this voxel (x,y) under the influence of all effective control points
        for (m=0;m<=3;m++)
          for (n=0;n<=3;n++)
            {
            indexEffectiveControlPointX = m+indexFirstEffectiveControlPointX;
            indexEffectiveControlPointY = n+indexFirstEffectiveControlPointY;
            
            if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY) )
              {
                dx += Bspline(u,m)*Bspline(v,n)*controlPoints[indexEffectiveControlPointX][indexEffectiveControlPointY].x;
                dy += Bspline(u,m)*Bspline(v,n)*controlPoints[indexEffectiveControlPointX][indexEffectiveControlPointY].y;
              }
            }
        
        //save deformation field
        defField[x][y].x = dx;
        defField[x][y].y = dy;
        //interpolate intensity (float)
        sliceA2BFloat[x][y] = bilinearInterpolation(sliceA, (float)x+dx, (float)y+dy, sliceSize.x, sliceSize.y);        
        if ( (sliceA2BFloat[x][y]==0.0)&&(x+dx<0 || x+dx>=(sliceSize.x-1) || y+dy<0 || y+dy>=(sliceSize.y-1)) )
          sliceA2BFloat[x][y] = MAX(0.0, sliceA[x][y]);
      } // for each point
      
  // smooth (twice) deformation in each direction
  smoothDeformationField(defField, defField, sliceSize, levelIndex, numLevels);
}

// sub-function 14-f
float GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(unsigned char ***featureMapA, unsigned char ***featureMapB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, int numFeatures, float ***featureMapA2BFloat, Fvector2d **defField, int levelIndex, int numLevels)
{
  int x,y;
  int indexEffectiveControlPointX, indexEffectiveControlPointY;
  int m,n;
  int featureIndex;
  float feature;
  float dx, dy;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY;
  float u,v;
  int numControlPointsX = (int)(ceil((float)sliceSize.x/distBetweenControlPoints));
  int numControlPointsY = (int)(ceil((float)sliceSize.y/distBetweenControlPoints));
  
  float totalEnergy=0.0;
  float energySquaredAtThisPoint;
  
  for (x=0;x<sliceSize.x;x++)
    for (y=0;y<sliceSize.y;y++)
	  { //for each point
		dx=0;
		dy=0;
		
		energySquaredAtThisPoint = 0.0;
		
		// find out effecting control points
		indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPoints)-1;
		indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPoints)-1;
		u = (float)x/(float)distBetweenControlPoints - floor((float)x/(float)distBetweenControlPoints);
		v = (float)y/(float)distBetweenControlPoints - floor((float)y/(float)distBetweenControlPoints);
		
		// calculate displacement at this voxel (x,y) under the influence of all effective control points
		for (m=0;m<=3;m++)
		  for (n=0;n<=3;n++)
		    {
			indexEffectiveControlPointX = m+indexFirstEffectiveControlPointX;
			indexEffectiveControlPointY = n+indexFirstEffectiveControlPointY;
			
			if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY) )
			  {
				dx += Bspline(u,m)*Bspline(v,n)*controlPoints[indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += Bspline(u,m)*Bspline(v,n)*controlPoints[indexEffectiveControlPointX][indexEffectiveControlPointY].y;
			  }
			}
		/*
		//interpolate intensity (float)
		sliceA2BFloat[x][y] = bilinearInterpolation(sliceA, (float)x+dx, (float)y+dy, sliceSize);		
		*/
		// save deformation field
		defField[x][y].x = dx;
		defField[x][y].y = dy;
		
		
		//interpolate features (float)
		for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
		  {
		  //feature = bilinearInterpolationFeatures(featureMapA[featureIndex], (float)x+dx, (float)y+dy, sliceSize);
		  feature = bilinearInterpolation(featureMapA[featureIndex], (float)x+dx, (float)y+dy, sliceSize.x, sliceSize.y);
		  if ( (feature==0.0)&&(x+dx<0 || x+dx>=(sliceSize.x-1) || y+dy<0 || y+dy>=(sliceSize.y-1)) )
		    feature = MAX(0.0, featureMapA[featureIndex][x][y]);
			
		  //printf("featureMapA2BFloat[%d][%d][%d] = %f\n", featureIndex, x,y,featureMapA2BFloat[featureIndex][x][y]);
		  featureMapA2BFloat[featureIndex][x][y] = feature;
		  energySquaredAtThisPoint += pow( (featureMapA2BFloat[featureIndex][x][y] - (float)featureMapB[featureIndex][x][y]), 2.0 );
		  }
		totalEnergy += energySquaredAtThisPoint;
	  } // for each point
	  
  // smooth (twice) deformation in each direction
  smoothDeformationField(defField, defField, sliceSize, levelIndex, numLevels);
  
  totalEnergy /= (float)numFeatures;
  return totalEnergy;
}

// sub-function 15
float FFDLocalEffect(unsigned char **sliceA, int xcoor, int ycoor, Fvector2d **controlPoints, Fvector2d **defField, Ivector2d sliceSize, int distBetweenControlPoints, int controlPointIndex1, int controlPointIndex2, float increment, char orientation, int tempOrNot, int featureIndex, int numFeatures)
{
  float intensity;
  int indexEffectiveControlPointX, indexEffectiveControlPointY;
  int m,n;
  float dx, dy;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY;
  float u,v;
  int numControlPointsX = (int)(ceil((float)sliceSize.x/distBetweenControlPoints));
  int numControlPointsY = (int)(ceil((float)sliceSize.y/distBetweenControlPoints));
  
  //calculate displacement at this voxel (xcoor, ycoor) under the influence of all effective control points
  dx=defField[xcoor][ycoor].x;
  dy=defField[xcoor][ycoor].y;
 
  // find out effecting control points
  indexFirstEffectiveControlPointX = (int)floor((float)xcoor/(float)distBetweenControlPoints)-1;
  indexFirstEffectiveControlPointY = (int)floor((float)ycoor/(float)distBetweenControlPoints)-1;
  u = (float)xcoor/(float)distBetweenControlPoints - floor((float)xcoor/(float)distBetweenControlPoints);
  v = (float)ycoor/(float)distBetweenControlPoints - floor((float)ycoor/(float)distBetweenControlPoints);
    
  m = controlPointIndex1-indexFirstEffectiveControlPointX;
  n = controlPointIndex2-indexFirstEffectiveControlPointY;


  switch (orientation)
    {
    case 'X':  // 0 for direction 'X'
      dx += Bspline(u,m)*Bspline(v,n)*increment;
      break;
      
    case 'Y':  // 1 for direction 'Y'
      dy += Bspline(u,m)*Bspline(v,n)*increment;
      break;
    
    default:
      break;
    }
  
  
  if (tempOrNot==0 & featureIndex==(numFeatures-1))  // 1) not temp, need to update defField 2) update deformation at this point only once, not numFeatures times !!!
    {
    defField[xcoor][ycoor].x = dx;
    defField[xcoor][ycoor].y = dy;
    }

  //interpolate intensity (float)
  intensity = bilinearInterpolation(sliceA, (float)xcoor+dx, (float)ycoor+dy, sliceSize.x, sliceSize.y);

  return intensity;
}


// sub-function 16
void CopyFloatImage2D(float **imgA, float**imgB, Ivector2d imgSize)    // copy A to B
{
  int i,j;
  
  for (i=0;i<imgSize.x;i++)
    for (j=0;j<imgSize.y;j++)
      {
      imgB[i][j] = imgA[i][j];
      }
}


//sub-function 16-f
void CopyFeatureMapsFloat(float ***featureMapOrig, float ***featureMapCopy, Ivector2d sliceSize, int numFeatures)
{
  int i,j;
  int featureIndex;
  
  for (i=0;i<sliceSize.x;i++)
    for (j=0;j<sliceSize.y;j++)
      for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
        {
        featureMapCopy[featureIndex][i][j] = featureMapOrig[featureIndex][i][j];
        }
}

// sub-function 17
float sign(float x)
{
  float signX;
  signX = (float)((x>0)-(x<0));
  return signX;
}

// sub-function 18
void CopyDisplacementAtControlPoints2D(Fvector2d **controlPoints, Fvector2d **controlPointsBackup, int numControlPointsX, int numControlPointsY)
{
  int i,j;
  
  for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
	  controlPointsBackup[i][j].x=controlPoints[i][j].x;
	  controlPointsBackup[i][j].y=controlPoints[i][j].y;
	  }
}

// sub-function 19
void calculateDisplacementIncrementAtControlPoints(Fvector2d **controlPointsBackup, Fvector2d **controlPoints, Fvector2d **controlPointsIncrement, int numControlPointsX, int numControlPointsY)
{
  int i,j;
  for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
      {
      controlPointsIncrement[i][j].x = controlPoints[i][j].x - controlPointsBackup[i][j].x;
      controlPointsIncrement[i][j].y = controlPoints[i][j].y - controlPointsBackup[i][j].y;
      }
}


//sub-function 20
void smoothDisplacementAtControlPoints2D(Fvector2d **controlPoints, Fvector2d **controlPointsSmoothed, int numControlPointsX, int numControlPointsY, int distBetweenControlPoints, int levelIndex, int numLevels, float factor)
{
  int i,j;
  float **smoothMat;
  float gaussianKernel;
  int gaussianRadius;
  float **controlPointsX, **controlPointsY;
  float **controlPointsXSmoothed, **controlPointsYSmoothed;
  Ivector2d inputSize, smoothMatSize;
  
  if (numControlPointsX<=15 && numControlPointsX*distBetweenControlPoints<120)
   {
   gaussianRadius = 1;
   gaussianKernel = 0.45 * factor;
   }
  else if (numControlPointsX>15 || numControlPointsX*distBetweenControlPoints>200)
   {
   gaussianRadius = 4;
   gaussianKernel = 5.5 * factor;
   }
  else
   {
    //printf("smoothing matrix size (5, 5)\n");
    gaussianRadius = 3;
    gaussianKernel = 4.2 * factor;  // the smoothing kernel is adaptive to the distance between two adjacent control points  
    
    if (gaussianKernel<0.35)
      gaussianKernel = 0.35;  // if kernel is less than 0.35, the gaussian template will be 1 at the center but 0 everywhere else. This should generally be avoided.
   }
  
  
  inputSize.x = numControlPointsX+2*gaussianRadius;
  inputSize.y = numControlPointsY+2*gaussianRadius;
  smoothMatSize.x = 2*gaussianRadius+1;
  smoothMatSize.y = 2*gaussianRadius+1;
  
  smoothMat = Falloc2d(smoothMatSize.x, smoothMatSize.y);
  gauss2D(smoothMatSize.x, gaussianKernel, smoothMatSize.y, gaussianKernel, 0, smoothMat);
  
  controlPointsX = Falloc2d(inputSize.x, inputSize.y);
  controlPointsY = Falloc2d(inputSize.x, inputSize.y);
  controlPointsXSmoothed = Falloc2d(inputSize.x, inputSize.y);
  controlPointsYSmoothed = Falloc2d(inputSize.x, inputSize.y);
  
  for (i=0;i<inputSize.x;i++)
    for (j=0;j<inputSize.y;j++)
      {
      if ( (i<gaussianRadius) | (j<gaussianRadius) | (i>=(numControlPointsX+gaussianRadius)) | (j>=(numControlPointsY+gaussianRadius)) )
        {
        controlPointsX[i][j] = 0.0;
        controlPointsY[i][j] = 0.0;
        controlPointsXSmoothed[i][j] = 0.0;
        controlPointsYSmoothed[i][j] = 0.0;
        }
      else
        {
        controlPointsX[i][j] = controlPoints[i-gaussianRadius][j-gaussianRadius].x;
        controlPointsY[i][j] = controlPoints[i-gaussianRadius][j-gaussianRadius].y;
        controlPointsXSmoothed[i][j] = 0.0;
        controlPointsYSmoothed[i][j] = 0.0;
        }
      }

  Convolution2DFloat(controlPointsX, inputSize, smoothMat, smoothMatSize, controlPointsXSmoothed);
  Convolution2DFloat(controlPointsY, inputSize, smoothMat, smoothMatSize, controlPointsYSmoothed);

  for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
      {
        controlPointsSmoothed[i][j].x = controlPointsXSmoothed[i+gaussianRadius][j+gaussianRadius];
        controlPointsSmoothed[i][j].y = controlPointsYSmoothed[i+gaussianRadius][j+gaussianRadius];
      }
      
  // release momery
  Ffree2d(controlPointsX, inputSize.x);
  Ffree2d(controlPointsY, inputSize.x);
  Ffree2d(controlPointsXSmoothed, inputSize.x);
  Ffree2d(controlPointsYSmoothed, inputSize.x);
  Ffree2d(smoothMat, smoothMatSize.x);
}



//sub-function 21
void UpdateControlPointsWithSmoothIncrement2(Fvector2d **controlPointsBackup, Fvector2d **controlPointsIncrement, Fvector2d **controlPointsIncrementSmoothed, Fvector2d **controlPoints, int numControlPointsX, int numControlPointsY, float alpha, unsigned char **mask, int distBetweenControlPoints, int levelIndex, int numLevels)
{
  int i,j, ii,jj;
  Ivector2d realCoordinate;
  float maxDisplacementAllowed;
  float offset;
  
  if (distBetweenControlPoints<14)
    maxDisplacementAllowed = distBetweenControlPoints * pow(1.35, numLevels+1-levelIndex);
  else
    maxDisplacementAllowed = distBetweenControlPoints * pow(1.25, numLevels-levelIndex);
    
  for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
      {
      realCoordinate.x = i*distBetweenControlPoints;
      realCoordinate.y = j*distBetweenControlPoints;
      
      controlPoints[i][j].x = controlPointsBackup[i][j].x + ((1-alpha)*controlPointsIncrement[i][j].x + alpha*controlPointsIncrementSmoothed[i][j].x);
      controlPoints[i][j].y = controlPointsBackup[i][j].y + ((1-alpha)*controlPointsIncrement[i][j].y + alpha*controlPointsIncrementSmoothed[i][j].y);
      }
    
  for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
      {
      // last step: check if displacement at a control point has exceeded the maximum allowed
      if (fabs(controlPoints[i][j].x) > maxDisplacementAllowed)
        {
        offset = controlPoints[i][j].x - sign(controlPoints[i][j].x) * maxDisplacementAllowed;
        controlPoints[i][j].x = sign(controlPoints[i][j].x) * maxDisplacementAllowed;
        
        ii = i+(int)sign(controlPoints[i][j].x);
        while (ii>=0 && ii<numControlPointsX)
          {
          //printf("X: %d->%d, offset = %f\n", i, ii, offset);
          controlPoints[ii][j].x += offset*pow(0.8, fabs((float)(ii-i)));
          if (fabs(controlPoints[ii][j].x) > maxDisplacementAllowed)
            {
            offset = controlPoints[ii][j].x - sign(controlPoints[ii][j].x)*maxDisplacementAllowed;
            controlPoints[ii][j].x = sign(controlPoints[ii][j].x) * maxDisplacementAllowed;
            ii += (int)sign(controlPoints[ii][j].x);
            }
          else 
            break;
          }//while
        } // if
      if (fabs(controlPoints[i][j].y) > maxDisplacementAllowed)
        {
        offset = controlPoints[i][j].y - sign(controlPoints[i][j].y)*maxDisplacementAllowed;
        controlPoints[i][j].y = sign(controlPoints[i][j].y) * maxDisplacementAllowed;
                
        jj = j+(int)sign(controlPoints[i][j].y);
        while (jj>=0 && jj<numControlPointsY)
          {
          //printf("Y: %d->%d, offset = %f\n", j, jj, offset);
          controlPoints[i][jj].y += offset*pow(0.8, fabs((float)(jj-j)));
          if (fabs(controlPoints[i][jj].y) > maxDisplacementAllowed)
            {
            offset = controlPoints[i][jj].y - sign(controlPoints[i][jj].y)*maxDisplacementAllowed;
            controlPoints[i][jj].y = sign(controlPoints[i][jj].y) * maxDisplacementAllowed;
            jj += (int)sign(controlPoints[i][jj].y);
            }
          else 
            break;
          }//while
        } //if
      }
}


// sub-function 22
void UpdateControlPointsWithSmoothIncrement(Fvector2d **controlPointsBackup, Fvector2d **controlPoints, Fvector2d **controlPointsUpdated, int numControlPointsX, int numControlPointsY, int IterIndex, int maxNumIterInResolution, int levelIndex, int numLevels, int distBetweenControlPoints, unsigned char **mask)
{
  // idea: controlPointsUpdated = controlPointsBackup + (1-alpha)*(controlPoints-controlPointsBackup) + alpha*smooth(controlPoints-controlPointsBackup)
  
  Fvector2d **controlPointsIncrement, **controlPointsIncrementSmoothed;
  controlPointsIncrement = Fvector2dalloc2d(numControlPointsX, numControlPointsY);
  controlPointsIncrementSmoothed = Fvector2dalloc2d(numControlPointsX, numControlPointsY);
  float alpha;   // smooth weighting
  
  calculateDisplacementIncrementAtControlPoints(controlPointsBackup, controlPoints, controlPointsIncrement, numControlPointsX, numControlPointsY);  // calculate increment: controlPoints - controlPointsBackup
  smoothDisplacementAtControlPoints2D(controlPointsIncrement, controlPointsIncrementSmoothed, numControlPointsX, numControlPointsY, distBetweenControlPoints, levelIndex, numLevels, 1.0); // smooth: controlPointsIncrement -> controlPointsIncrementSmoothed
     
  alpha=0.5*exp(-pow((float)IterIndex,2.0)/(2*pow((float)maxNumIterInResolution,2.0))); // Parameter "alpha" controls how much to smooth: the bigger alpha is, the more smoothness there will be. At coarse stage, smooth more; as displacement gets more and more accurate at each control point, smooth less
  UpdateControlPointsWithSmoothIncrement2(controlPointsBackup, controlPointsIncrement, controlPointsIncrementSmoothed, controlPointsUpdated, numControlPointsX, numControlPointsY, alpha, mask, distBetweenControlPoints, levelIndex, numLevels);
  
  //release memory
  Fvector2dfree2d(controlPointsIncrement, numControlPointsX);
  Fvector2dfree2d(controlPointsIncrementSmoothed, numControlPointsX);
}

// sub-function 23
void LinearCombinationOfTwoControlPointsMats(Fvector2d **controlPointsNew, Fvector2d **controlPointsOld, Fvector2d **controlPointsUpdated, float weighting, int numControlPointsX, int numControlPointsY)
{
  int i,j;
  
  for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
      {
      controlPointsUpdated[i][j].x = weighting*controlPointsOld[i][j].x + (1-weighting)*controlPointsNew[i][j].x;
      controlPointsUpdated[i][j].y = weighting*controlPointsOld[i][j].y + (1-weighting)*controlPointsNew[i][j].y;
      }
}

// sub-function 24
void Convolution2DFloat(float** input,Ivector2d inputSize,float** mask,Ivector2d maskSize,float** output)
{
  int row,col,ii,jj,iii,jjj;
  Ivector2d halfSize,evenOrNot;
  double sum;

  halfSize.x=(maskSize.x-1)/2;
  halfSize.y=(maskSize.y-1)/2;

  evenOrNot.x=-(maskSize.x-1)%2;
  evenOrNot.y=-(maskSize.y-1)%2;

  for (row=0;row<inputSize.x;row++)
    for (col=0;col<inputSize.y;col++)
      {
      sum=0;
	for (ii=-halfSize.x-evenOrNot.x; ii<=halfSize.x; ii++)
	  {
	  iii = row+ii;
	  if (iii<0)			iii = -iii -1;
	  else if (iii>=inputSize.x)	iii = (inputSize.x-1) - (iii-inputSize.x);

	  for (jj=-halfSize.y-evenOrNot.y; jj<=halfSize.y; jj++)
	    {
	    jjj = col+jj;
	    if (jjj<0) 			jjj = -jjj -1;
	    if (jjj>=inputSize.y)	jjj = (inputSize.y-1) - (jjj-inputSize.y);

	    sum += input[iii][jjj] * mask[ii+halfSize.x+evenOrNot.x][jj+halfSize.y+evenOrNot.y];
	    } // jj
	  } // ii
	output[row][col]=(float)sum;
      }
}


// sub-function 26
void smoothDeformationField(Fvector2d **dfField, Fvector2d **dfFieldSmoothed, Ivector2d dfSize, int levelIndex, int numLevels)
{
  int i,j;
  float **smoothMat;
  int gaussianRadius;
  float gaussianKernel;
  
  switch (levelIndex)
    {
      case 4:  //coarsest level
        gaussianRadius = 1;
        gaussianKernel = 1.0;
        break;
      case 3:  // coarse level
        gaussianRadius = 2;
        gaussianKernel = 1.25;
        break;
      case 2:  // middle level
        gaussianRadius = 3;
        gaussianKernel = 2.2;
        break;
      case 1:
        gaussianRadius = 4;
        gaussianKernel = 3.0;
        break;
    
      default:
        break;
    }
  
  smoothMat = Falloc2d(2*gaussianRadius+1, 2*gaussianRadius+1);
  gauss2D(2*gaussianRadius+1, gaussianKernel, 2*gaussianRadius+1, gaussianKernel, 0, smoothMat);
  
  float **dfFieldX, **dfFieldY;
  float **dfFieldXSmoothed, **dfFieldYSmoothed;
  dfFieldX = Falloc2d(dfSize.x, dfSize.y);
  dfFieldY = Falloc2d(dfSize.x, dfSize.y);
  dfFieldXSmoothed = Falloc2d(dfSize.x, dfSize.y);
  dfFieldYSmoothed = Falloc2d(dfSize.x, dfSize.y);
  for (i=0;i<dfSize.x;i++)
    for (j=0;j<dfSize.y;j++)
      {
      dfFieldX[i][j] = dfField[i][j].x;
      dfFieldY[i][j] = dfField[i][j].y;
      }
      
  Ivector2d smoothMatSize;
  smoothMatSize.x = 2*gaussianRadius+1;
  smoothMatSize.y = 2*gaussianRadius+1;
  
  Convolution2DFloat(dfFieldX, dfSize, smoothMat, smoothMatSize, dfFieldXSmoothed);
  Convolution2DFloat(dfFieldY, dfSize, smoothMat, smoothMatSize, dfFieldYSmoothed);
  
  for (i=0;i<dfSize.x;i++)
    for (j=0;j<dfSize.y;j++)
      {
        dfFieldSmoothed[i][j].x = dfFieldXSmoothed[i][j];
        dfFieldSmoothed[i][j].y = dfFieldYSmoothed[i][j];
      }
      
  // release memory
  Ffree2d(dfFieldX, dfSize.x);
  Ffree2d(dfFieldY, dfSize.x);
  Ffree2d(dfFieldXSmoothed, dfSize.x);
  Ffree2d(dfFieldYSmoothed, dfSize.x);
}

// sub-function 27
void UpsampleDisplacementAtControlPoints(Fvector2d **controlPointsPreviousLevel, Fvector2d **controlPointsThisLevel, int numControlPointsXPreviousLevel, int numControlPointsYPreviousLevel, int numControlPointsXThisLevel, int numControlPointsYThisLevel, int distBetweenControlPoints, unsigned char **maksForControlPointsThisLevel, int levelIndex, int numLevels)
{
  int i,j, ii,jj, s,t;
  int iceil, ifloor, jceil, jfloor;
  Ivector2d realCoordinate;
  float iOffsetFloor, iOffsetCeil, jOffsetFloor, jOffsetCeil;
  float offset;
  Fvector2d **controlPointsThisLevelBackup;
  Fvector2d **controlPointsTemp1, **controlPointsTemp2;
  controlPointsThisLevelBackup = Fvector2dalloc2d(numControlPointsXThisLevel, numControlPointsYThisLevel);
  controlPointsTemp1 = Fvector2dalloc2d(numControlPointsXThisLevel, numControlPointsYThisLevel);
  controlPointsTemp2 = Fvector2dalloc2d(numControlPointsXThisLevel, numControlPointsYThisLevel);
  
  float maxDisplacementAllowed;
  if (distBetweenControlPoints<14)
    maxDisplacementAllowed = distBetweenControlPoints * pow(1.35, numLevels+1-levelIndex);
  else
    maxDisplacementAllowed = distBetweenControlPoints * pow(1.25, numLevels-levelIndex);
    
  
  // upsample and interpolate
  for (i=0;i<MIN(numControlPointsXPreviousLevel*2, numControlPointsXThisLevel-1);i++)
    for (j=0;j<MIN(numControlPointsYPreviousLevel*2, numControlPointsYThisLevel-1);j++)
      {
      realCoordinate.x = i*distBetweenControlPoints;
      realCoordinate.y = j*distBetweenControlPoints;
     
      //printf("(i,j) = (%d, %d)\n",i,j);
      if (i%2==0 & j%2==0) // both even numbers
        {
        controlPointsThisLevelBackup[i][j].x = controlPointsPreviousLevel[i/2][j/2].x;
        controlPointsThisLevelBackup[i][j].y = controlPointsPreviousLevel[i/2][j/2].y;
        }
      else if (i%2==0 & j%2!=0) // i even, j odd
        {
        jfloor = (int)floor((float)j/2.0);
        jceil = (int)ceil((float)j/2.0);
        if (jceil>(numControlPointsYPreviousLevel-1))
          {
          controlPointsThisLevelBackup[i][j].x = controlPointsPreviousLevel[i/2][jfloor].x;
          controlPointsThisLevelBackup[i][j].y = controlPointsPreviousLevel[i/2][jfloor].y;
          }
        else
          {
          jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
          jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
          controlPointsThisLevelBackup[i][j].x = (jOffsetFloor*controlPointsPreviousLevel[i/2][jfloor].x + jOffsetCeil*controlPointsPreviousLevel[i/2][jceil].x);
          controlPointsThisLevelBackup[i][j].y = (jOffsetFloor*controlPointsPreviousLevel[i/2][jfloor].y + jOffsetCeil*controlPointsPreviousLevel[i/2][jceil].y);
          }
        }
      else if (i%2!=0 & j%2==0) // i odd, j even
        {
        ifloor = (int)floor((float)i/2.0);
        iceil = (int)ceil((float)i/2.0);
        if (iceil>(numControlPointsXPreviousLevel-1))
          {
          controlPointsThisLevelBackup[i][j].x = controlPointsPreviousLevel[ifloor][j/2].x;
          controlPointsThisLevelBackup[i][j].y = controlPointsPreviousLevel[ifloor][j/2].y;
          }
        else
          {
          iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
          iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
          controlPointsThisLevelBackup[i][j].x = (iOffsetFloor*controlPointsPreviousLevel[ifloor][j/2].x + iOffsetCeil*controlPointsPreviousLevel[iceil][j/2].x);
          controlPointsThisLevelBackup[i][j].y = (iOffsetFloor*controlPointsPreviousLevel[ifloor][j/2].y + iOffsetCeil*controlPointsPreviousLevel[iceil][j/2].y);
          }
        } //else if (i%2!=0 & j%2==0) // i odd, j even
      else if (i%2!=0 & j%2!=0) // both odd
        {
        ifloor = (int)floor((float)i/2.0);
        iceil = (int)ceil((float)i/2.0);
        jfloor = (int)floor((float)j/2.0);
        jceil = (int)ceil((float)j/2.0);
        if (iceil>(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1))
          {
          controlPointsThisLevelBackup[i][j].x = controlPointsPreviousLevel[ifloor][jfloor].x;
          controlPointsThisLevelBackup[i][j].y = controlPointsPreviousLevel[ifloor][jfloor].y;
          }
        else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1))
          {
          iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
          iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
          controlPointsThisLevelBackup[i][j].x = (iOffsetFloor*controlPointsPreviousLevel[ifloor][jfloor].x + iOffsetCeil*controlPointsPreviousLevel[iceil][jfloor].x);
          controlPointsThisLevelBackup[i][j].y = (iOffsetFloor*controlPointsPreviousLevel[ifloor][jfloor].y + iOffsetCeil*controlPointsPreviousLevel[iceil][jfloor].y);
          }
        else if (iceil>(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1))
          {
          jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
          jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
          controlPointsThisLevelBackup[i][j].x = (jOffsetFloor*controlPointsPreviousLevel[ifloor][jfloor].x + jOffsetCeil*controlPointsPreviousLevel[ifloor][jceil].x);
          controlPointsThisLevelBackup[i][j].y = (jOffsetFloor*controlPointsPreviousLevel[ifloor][jfloor].y + jOffsetCeil*controlPointsPreviousLevel[ifloor][jceil].y);
          }
        else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1))
          {
          iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
          iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
          jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
          jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
          controlPointsThisLevelBackup[i][j].x = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[ifloor][jfloor].x + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[ifloor][jceil].x + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[iceil][jfloor].x + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[iceil][jceil].x);
          controlPointsThisLevelBackup[i][j].y = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[ifloor][jfloor].y + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[ifloor][jceil].y + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[iceil][jfloor].y + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[iceil][jceil].y);
          }
        } // else if (i%2!=0 & j%2!=0) // both odd
        
        
        controlPointsThisLevelBackup[i][j].x *= 2;
        controlPointsThisLevelBackup[i][j].y *= 2;
      }// for each control point this level
  
  // step 1: copy controlPointsBackup to controlPoints
  for(i=0;i<numControlPointsXThisLevel;i++)
    for (j=0;j<numControlPointsYThisLevel;j++)
      {
      controlPointsTemp1[i][j].x = controlPointsThisLevelBackup[i][j].x;
      controlPointsTemp1[i][j].y = controlPointsThisLevelBackup[i][j].y;
      controlPointsTemp2[i][j].x = controlPointsThisLevelBackup[i][j].x;
      controlPointsTemp2[i][j].y = controlPointsThisLevelBackup[i][j].y;
      }

    // step 2: if displacement at a control point exceeds maximum allowed, distribute it along the gradient direction.
    for (i=0;i<numControlPointsXThisLevel;i++)
      for (j=0;j<numControlPointsYThisLevel;j++)
        {
        // x directions, -
        if (fabs(controlPointsTemp1[i][j].x)>maxDisplacementAllowed)
          {
          offset = controlPointsTemp1[i][j].x - sign(controlPointsTemp1[i][j].x)*maxDisplacementAllowed;
          controlPointsTemp1[i][j].x = sign(controlPointsTemp1[i][j].x)*maxDisplacementAllowed;
          
          ii = i-1;
          while (ii>=0 && ii<numControlPointsXThisLevel)
             {
             controlPointsTemp1[ii][j].x += offset*pow(0.75, fabs((float)(ii-i)));
             if (fabs(controlPointsTemp1[ii][j].x)>maxDisplacementAllowed)
               {
               offset = controlPointsTemp1[ii][j].x - sign(controlPointsTemp1[i][j].x)*maxDisplacementAllowed;
               controlPointsTemp1[ii][j].x = sign(controlPointsTemp1[i][j].x)*maxDisplacementAllowed;
               ii += -1;
               } // if(ii,j) 
             else
               break;
             } // while
          }// if (i,j)
          
        // y direction, -
        if (fabs(controlPointsTemp1[i][j].y)>maxDisplacementAllowed)
          {
          offset = controlPointsTemp1[i][j].y - sign(controlPointsTemp1[i][j].y)* maxDisplacementAllowed;
          controlPointsTemp1[i][j].y = sign(controlPointsTemp1[i][j].y)*maxDisplacementAllowed;
          
          jj = j-1;
          while (jj>=0 && jj<numControlPointsYThisLevel)
             {
             controlPointsTemp1[i][jj].y += offset*pow(0.75, fabs((float)(jj-j)));
             if (fabs(controlPointsTemp1[i][jj].y)>maxDisplacementAllowed)
               {
               offset = controlPointsTemp1[i][jj].y - sign(controlPointsTemp1[i][j].y)*maxDisplacementAllowed;
               controlPointsTemp1[i][jj].y = sign(controlPointsTemp1[i][j].y)*maxDisplacementAllowed;
               jj += -1;
               } // if(i,jj) 
             else
               break;
             } // while
          }// if (i,j)
        } // for for
    
    for (i=(numControlPointsXThisLevel-1);i>=0;i--)
      for (j=(numControlPointsYThisLevel-1);j>=0;j--)
        {
        // x directions, +
        if (fabs(controlPointsTemp2[i][j].x)>maxDisplacementAllowed)
          {
          offset = controlPointsTemp2[i][j].x - sign(controlPointsTemp2[i][j].x)*maxDisplacementAllowed;
          controlPointsTemp2[i][j].x = sign(controlPointsTemp2[i][j].x)*maxDisplacementAllowed;
          
          ii = i+1;
          while (ii>=0 && ii<numControlPointsXThisLevel)
             {
             controlPointsTemp2[ii][j].x += offset*pow(0.75, fabs((float)(ii-i)));
             if (fabs(controlPointsTemp2[ii][j].x)>maxDisplacementAllowed)
               {
               offset = controlPointsTemp2[ii][j].x - sign(controlPointsTemp2[i][j].x)*maxDisplacementAllowed;
               controlPointsTemp2[ii][j].x = sign(controlPointsTemp2[i][j].x)*maxDisplacementAllowed;
               ii += 1;
               } // if(ii,j) 
             else
               break;
             } // while
          }// if (i,j)
          
        // y direction, +
        if (fabs(controlPointsTemp2[i][j].y)>maxDisplacementAllowed)
          {
          offset = controlPointsTemp2[i][j].y - sign(controlPointsTemp2[i][j].y)*maxDisplacementAllowed;
          controlPointsTemp2[i][j].y = sign(controlPointsTemp2[i][j].y)*maxDisplacementAllowed;
          
          jj = j+1;
          while (jj>=0 && jj<numControlPointsYThisLevel)
             {
             controlPointsTemp2[i][jj].y += offset*pow(0.75, fabs((float)(jj-j)));
             if (fabs(controlPointsTemp2[i][jj].y)>maxDisplacementAllowed)
               {
               offset = controlPointsTemp2[i][jj].y - sign(controlPointsTemp2[i][j].y)*maxDisplacementAllowed;
               controlPointsTemp2[i][jj].y = sign(controlPointsTemp2[i][j].y)*maxDisplacementAllowed;
               jj += 1;
               } // if(i,jj) 
             else
               break;
             } // while
          }// if (i,j)
       } // for for
    
    for (i=0;i<numControlPointsXThisLevel;i++)
      for (j=0;j<numControlPointsYThisLevel;j++)
        {
        // x
        if (fabs(controlPointsTemp1[i][j].x)>fabs(controlPointsTemp2[i][j].x))
          controlPointsThisLevel[i][j].x = 0.7 * controlPointsTemp1[i][j].x + 0.3* controlPointsTemp2[i][j].x;
        else
          controlPointsThisLevel[i][j].x = 0.3 * controlPointsTemp1[i][j].x + 0.7* controlPointsTemp2[i][j].x;
          
        //y
        if (fabs(controlPointsTemp1[i][j].y)>fabs(controlPointsTemp2[i][j].y))
          controlPointsThisLevel[i][j].y = 0.7 * controlPointsTemp1[i][j].y + 0.3* controlPointsTemp2[i][j].y;
        else
          controlPointsThisLevel[i][j].y = 0.3 * controlPointsTemp1[i][j].y + 0.7* controlPointsTemp2[i][j].y;
        }
    
    Fvector2dfree2d(controlPointsThisLevelBackup, numControlPointsXThisLevel);
    Fvector2dfree2d(controlPointsTemp1, numControlPointsXThisLevel);
    Fvector2dfree2d(controlPointsTemp2, numControlPointsXThisLevel);
}



//sub-function 28
void GenerateMaskForControlPoints(unsigned char **sliceA, unsigned char **sliceB, Ivector2d sliceSize, int distBetweenControlPoints, unsigned char **mask, int levelIndex, int numLevels)
{
   int i,j, x,y;
   int dilateRadius;
   float dilationFactor;
   int s;
      
   printf("\nGenerating a mask for control points ...");
   
   if (distBetweenControlPoints<10)
     dilationFactor = pow(1.6, (numLevels-levelIndex+1));
   else
     dilationFactor = pow(1.45, (numLevels-levelIndex+1));
     
   dilateRadius = (int)(dilationFactor*distBetweenControlPoints);
   for (i=0; i<sliceSize.x; i++)
     for (j=0; j<sliceSize.y; j++)
       {  
        mask[i][j]=0;
        
        if ( sliceA[i][j]>0 || sliceB[i][j]>0 )
          mask[i][j]=255;
        else
          {
          for (s=0;s<32;s++)
            {
            x = (int)(i + dilateRadius*cos(s*G_PI/16));
            y = (int)(j + dilateRadius*sin(s*G_PI/16));
            
            if ( sliceA[MIN((sliceSize.x-1),MAX(x,0))][MIN((sliceSize.y-1), MAX(0,y))]>0 || sliceB[MIN((sliceSize.x-1),MAX(x,0))][MIN((sliceSize.y-1), MAX(0,y))]>0 )
                {
                  mask[i][j]=255;
                  break;
                } //if                
            } // for 
          } // else
       } // for each point
   printf("done!\n");
}




// sub-function 29
void saveDisplacementAtControlPoints2D(Fvector2d **controlPoints, int numControlPointsXThisLevel, int numControlPointsYThisLevel)
{
   int i, j;
   char fileNameX[1024], fileNameY[1024];
   FILE *fp1, *fp2;
   
   sprintf(fileNameX, "%s", "dvXMat.txt");
   sprintf(fileNameY, "%s", "dvYMat.txt");
   
   printf("Saving displacement at control points...\n");
   fp1 = fopen(fileNameX, "w");
   fp2 = fopen(fileNameY, "w");
   for (i=0;i<numControlPointsXThisLevel;i++)
     {
     printf("i=%d\n", i);
     for (j=0;j<numControlPointsYThisLevel;j++)
       {
        fprintf(fp1, "%f ", controlPoints[i][j].y);
        fprintf(fp2, "%f ", controlPoints[i][j].x);
       }
     }
  
  fclose(fp1);
  fclose(fp2);
}



// sub-function 30
void SmoothUCImage2D(unsigned char **slice, unsigned char **sliceSmoothed, Ivector2d sliceSize)
{
  int i,j;
  float **smoothMat;
  Ivector2d inputSize, smoothMatSize;
  float **sliceFloat, **sliceSmoothedFloat;
  int gaussianRadius = 2;
  float gaussianKernel = 1.4;
  
  if (sliceSize.x>100 && sliceSize.x<=200)
    {
    gaussianRadius = 3;
    gaussianKernel = 2;
    }
  else if (sliceSize.x>200)
    {
    gaussianRadius = 4;
    gaussianKernel = 2.5;
    }
  
  inputSize.x = sliceSize.x+2*gaussianRadius;
  inputSize.y = sliceSize.y+2*gaussianRadius;
  smoothMatSize.x = 2*gaussianRadius+1;
  smoothMatSize.y = 2*gaussianRadius+1;
  
  smoothMat = Falloc2d(smoothMatSize.x, smoothMatSize.y);
  gauss2D(smoothMatSize.x, gaussianKernel, smoothMatSize.y, gaussianKernel, 0, smoothMat);
  
  sliceFloat = Falloc2d(inputSize.x, inputSize.y);
  sliceSmoothedFloat = Falloc2d(inputSize.x, inputSize.y);
  
  for (i=0;i<inputSize.x;i++)
    for (j=0;j<inputSize.y;j++)
      {
      if ( (i<gaussianRadius) | (j<gaussianRadius) | (i>=(sliceSize.x+gaussianRadius)) | (j>=(sliceSize.y+gaussianRadius)) )
        {
        sliceFloat[i][j] = 0.0;
        sliceSmoothedFloat[i][j] = 0.0;
        }
      else
        {
        sliceFloat[i][j] = (float)slice[i-gaussianRadius][j-gaussianRadius];
        sliceSmoothedFloat[i][j] = 0.0;
        }
      }
      
  Convolution2DFloat(sliceFloat, inputSize, smoothMat, smoothMatSize, sliceSmoothedFloat);
  
  for (i=0;i<sliceSize.x;i++)
    for (j=0;j<sliceSize.y;j++)
      sliceSmoothed[i][j] = (int)(sliceSmoothedFloat[i+gaussianRadius][j+gaussianRadius]+0.5);
      
  // release momery
  Ffree2d(sliceSmoothedFloat, inputSize.x);
  Ffree2d(sliceFloat, inputSize.x);
  Ffree2d(smoothMat, smoothMatSize.x);
}




// sub-function 31
bool ReadFeaturesFromFeatureList(const char* featureImageListFile,unsigned char ***featureMap, int x_size,int y_size)
{
  FILE* fp=fopen(featureImageListFile,"r");

  int featureNum = 0;
  if (fscanf(fp,"%d",&featureNum) != 1) {
    fprintf(stderr, "Failed to read number of feature images from file %s!\n", featureImageListFile);
    fclose(fp);
    return false;
  }
 
  int    featureIndex = 0;
  Image* featureimage = NULL;
  string featureDir   = ::basis::os::path::dirname(featureImageListFile);
  char   filename[1024];

  for (featureIndex = 0; featureIndex < featureNum; featureIndex++) {
    if (fscanf(fp, "%s", filename) != 1) {
      fprintf(stderr, "Failed to read filename of feature number %d from file %s!\n", featureIndex + 1, featureImageListFile);
      break;
    }
    string filepath = ::basis::os::path::join(::basis::os::path::abspath(featureDir.c_str()), filename);
    featureimage = ReadNiftiImage(filepath.c_str());
    if (featureimage == NULL) {
      fprintf(stderr, "Failed to read feature image %d from file %s!\n", featureIndex + 1, filepath.c_str());
      break;
    }
    if (featureimage->hdr.datatype != DT_UNSIGNED_CHAR) {
      fprintf(stderr, "Feature image %s has datatype %d, but DT_UNSIGNED_CHAR is required!\n", filepath.c_str(), featureimage->hdr.datatype);
      break;
    }
    featureimage->SetFormat(Image::FORMAT_DRAMMS);
    if (featureimage->region.nz > 1) {
      fprintf(stderr, "Feature image %d read from file %s has more than two dimensions!\n", featureIndex + 1, filepath.c_str());
      break;
    }
    if (featureimage->region.nx != x_size || featureimage->region.ny != y_size) {
      fprintf(stderr, "Feature image %d read from file %s has size [%d, %d], but size [%d, %d] is required!\n",
            featureIndex + 1, filepath.c_str(),
            featureimage->region.nx,
            featureimage->region.ny,
            x_size, y_size);
      break;
    }
    for (int i = 0; i < x_size;i++) {
      for (int j = 0; j < y_size;j++) {
        featureMap[featureIndex][i][j] = featureimage->img.uc[0][i][j];
      }
    }
    delete featureimage;
    featureimage = NULL;
  }

  if (featureimage) delete featureimage;
  fclose(fp);

  return (featureIndex == featureNum); // all feature images read successfully
}


// sub-function 32
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(unsigned char ***featureMapA, unsigned char ***featureMapB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defFiled, int numFeatures, float ***featureMapA2BUpdated, float ***featureMapA2BOld, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind, int controlPointIndex1, int controlPointIndex2, float increment, char orientation, int NthLevel, int updateFeaturesOrNot)
{
  int i,j, featureIndex;
  int xcoor, ycoor;
  float deltaEnergySquared;
  float energyChange=0;
  
  for (i=-affectingRadiusAhead;i<affectingRadiusBehind;i++)
    for (j=-affectingRadiusAhead;j<affectingRadiusBehind;j++)
      {
        xcoor = realCoorThisControlPoint.x + j;
        ycoor = realCoorThisControlPoint.y + i;
        
        if (xcoor>=0 && xcoor<sliceSize.x && ycoor>=0 && ycoor<sliceSize.y)
          {
          deltaEnergySquared = 0.0;
          for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
            {
            featureMapA2BUpdated[featureIndex][xcoor][ycoor] = FFDLocalEffect(featureMapA[featureIndex], xcoor, ycoor, controlPoints, defFiled, sliceSize, distBetweenControlPoints, controlPointIndex1, controlPointIndex2, increment, orientation, NNO, featureIndex, numFeatures);
            deltaEnergySquared += (featureMapA2BUpdated[featureIndex][xcoor][ycoor]+featureMapA2BOld[featureIndex][xcoor][ycoor]-2*(float)featureMapB[featureIndex][xcoor][ycoor])*(featureMapA2BUpdated[featureIndex][xcoor][ycoor]-featureMapA2BOld[featureIndex][xcoor][ycoor]);
            
            if (updateFeaturesOrNot==YYES)
              featureMapA2BOld[featureIndex][xcoor][ycoor] = featureMapA2BUpdated[featureIndex][xcoor][ycoor];
            }
            energyChange += deltaEnergySquared;
          }
      } // calculate energy change brought by deltaX
    
  energyChange /= numFeatures;
  
  return energyChange;
}


// sub-function 32 -(2)
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(unsigned char ***featureMapA, unsigned char ***featureMapB, Ivector2d sliceSize, int distBetweenControlPoints, Fvector2d **controlPoints, Fvector2d **defFiled, int numFeatures, float ***featureMapA2BUpdated, float ***featureMapA2BOld, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind, int controlPointIndex1, int controlPointIndex2, float increment, char orientation, int NthLevel)
{
  int i,j, featureIndex;
  int xcoor, ycoor;
  float deltaEnergySquared;
  float energyChange=0;
  
  // to speed up, sample the voxels in the neighborhood of control point (controlPointIndex1, controlPointIndex2)
  int startPoint = -(int)((float)affectingRadiusAhead*pow(0.7, NthLevel));
  int endPoint = (int)((float)affectingRadiusBehind*pow(0.7, NthLevel));
  int interval = (int)ceil((float)(endPoint-startPoint)/5.0);
  int numPointsSampled=0;
  
  for (i=startPoint;i<=endPoint;i=i+interval)
    for (j=startPoint;j<=endPoint;j=j+interval)
      {
        xcoor = realCoorThisControlPoint.x + j;
        ycoor = realCoorThisControlPoint.y + i;
        
        if (xcoor>=0 && xcoor<sliceSize.x && ycoor>=0 && ycoor<sliceSize.y)
          {
          deltaEnergySquared = 0.0;
          for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
            {
            featureMapA2BUpdated[featureIndex][xcoor][ycoor] = FFDLocalEffect(featureMapA[featureIndex], xcoor, ycoor, controlPoints, defFiled, sliceSize, distBetweenControlPoints, controlPointIndex1, controlPointIndex2, increment, orientation, YYES, featureIndex, numFeatures);
            deltaEnergySquared += (featureMapA2BUpdated[featureIndex][xcoor][ycoor]+featureMapA2BOld[featureIndex][xcoor][ycoor]-2*(float)featureMapB[featureIndex][xcoor][ycoor])*(featureMapA2BUpdated[featureIndex][xcoor][ycoor]-featureMapA2BOld[featureIndex][xcoor][ycoor]);
            }
          energyChange += deltaEnergySquared;
          numPointsSampled++;
          }
      } // calculate energy change brought by deltaX

  float samplingRate = pow(0.85,NthLevel)*(float)distBetweenControlPoints*(float)distBetweenControlPoints*16.0/(float)numPointsSampled;   
  energyChange *= (samplingRate/(float)numFeatures);
  
  return energyChange;
}


// sub-function 33
void CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(float ***featureMapOrig, float ***featureMapCopied, Ivector2d sliceSize, int numFeatures, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind)
{
  int i,j, featureIndex;
  int xcoor, ycoor;
  
  for (i=-affectingRadiusAhead;i<affectingRadiusBehind;i++)
    for (j=-affectingRadiusAhead;j<affectingRadiusBehind;j++)
      {
        xcoor = realCoorThisControlPoint.x + j;
        ycoor = realCoorThisControlPoint.y + i;
        
        if (xcoor>=0 && xcoor<sliceSize.x && ycoor>=0 && ycoor<sliceSize.y)
          {
            for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
              featureMapCopied[featureIndex][xcoor][ycoor] = featureMapOrig[featureIndex][xcoor][ycoor];
          }
      }
  
}


// sub-function 34
void CopyDispAtTheBlockInducedByThisControlPoint(Fvector2d **defFieldOrig, Fvector2d **defFieldCopied, Ivector2d sliceSize, Ivector2d realCoorThisControlPoint, int affectingRadiusAhead, int affectingRadiusBehind)
{
  int i,j;
  int xcoor, ycoor;
  
  for (i=-affectingRadiusAhead;i<affectingRadiusBehind;i++)
    for (j=-affectingRadiusAhead;j<affectingRadiusBehind;j++)
      {
        xcoor = realCoorThisControlPoint.x + j;
        ycoor = realCoorThisControlPoint.y + i;
        
        if (xcoor>=0 && xcoor<sliceSize.x && ycoor>=0 && ycoor<sliceSize.y)
            {
            defFieldCopied[xcoor][ycoor].x = defFieldOrig[xcoor][ycoor].x;
            defFieldCopied[xcoor][ycoor].y = defFieldOrig[xcoor][ycoor].y;
            }
      } 
}



// sub-function 35
void GenerateFloatImageFromDeformationField(unsigned char **sliceA, unsigned char **sliceB, float **sliceA2BFloat, Fvector2d **defField, Ivector2d sliceSize)
{
  int x,y;
  
  for (x=0;x<sliceSize.x;x++)
    for (y=0;y<sliceSize.y;y++)
      {
      sliceA2BFloat[x][y] = bilinearInterpolation(sliceA, (float)x+defField[x][y].x, (float)y+defField[x][y].y, sliceSize.x, sliceSize.y);
      
      if (  (sliceA2BFloat[x][y]==0.0) && ( ((x+defField[x][y].x)<0) || ((x+defField[x][y].x)>=(sliceSize.x-1)) || ((y+defField[x][y].y)<0) || ((y+defField[x][y].y)>=(sliceSize.y-1)) )   )
        sliceA2BFloat[x][y] = MAX(0.0, sliceA[x][y]);
      }
      
}



// sub-function 36

float bilinearInterpolation(unsigned char **img, float x, float y, int sizeX, int sizeY)
{
  if ( x<0.0 || y<0.0 || x>(float)(sizeX-1) || y>(float)(sizeY-1) ) return 0.0;

  float intensity=0.0;
  int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
  float distXtoXFloor, distXtoXCeil, distYtoYFloor, distYtoYCeil;

  xFloor = static_cast<int>(floor(x));
  xCeil  = static_cast<int>(ceil(x));
  yFloor = static_cast<int>(floor(y));
  yCeil  = static_cast<int>(ceil(y));


  if ( (xFloor==xCeil)&(xFloor!=(sizeX-1)) )  xCeil+=1;
  if ( (xFloor==xCeil)&(xFloor==(sizeX-1)) )  xFloor-=1;  if (xFloor<0) xFloor=0;
  if ( (yFloor==yCeil)&(yFloor!=(sizeY-1)) )  yCeil+=1;
  if ( (yFloor==yCeil)&(yFloor==(sizeY-1)) )  yFloor-=1;  if (yFloor<0) yFloor=0;


  distXtoXFloor = x-xFloor;
  distXtoXCeil  = xCeil-x;
  distYtoYFloor = y-yFloor;
  distYtoYCeil  = yCeil-y;

  intensity = img[xFloor][yFloor]*distXtoXCeil*distYtoYCeil +
                          img[xFloor][yCeil]*distXtoXCeil*distYtoYFloor +
                          img[xCeil][yFloor]*distXtoXFloor*distYtoYCeil +
                          img[xCeil][yCeil]*distXtoXFloor*distYtoYFloor;

  return intensity;
}

