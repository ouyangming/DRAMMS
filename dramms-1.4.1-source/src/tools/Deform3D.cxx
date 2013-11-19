/**
 * @file  Deform3D.cxx
 * @brief Deformabley register 3D images given a list of feature images in both image spaces.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <Fast_PD.h>

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>

#include <common/matrix.h>  
#include <common/imageio.h>
#include <common/general.h>
#include <common/mvcd.h>
#include <common/cres.h>
#include <common/image.h>
#include <common/utilities.h>

#include "DiscreteOptimizationUtilities.h"
#include "GeneralUtilities.h"

#include <dramms/basis.h> // exename(), print_contact()


// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;


#define SSD  0
#define CC   1
#define MINCPDISTXY 3
#define MINCPDISTZ  1


// ===========================================================================
// help
// ===========================================================================

// ---------------------------------------------------------------------------
void print_help()
{
    string exec_name = exename();
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_imageA> <input_imageB> <prefix_featureFile> <output_imageA2B> <def_field> " << endl;
    cout << endl;
    cout << "Description:" << endl;
    cout << "  This program outputs the registered image and the corresponding deformation" << endl;
    cout << "  field given two 3D images and files listing feature images computed from these." << endl;
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
    cout << "  -b <int>,<int>,<int>  The number of voxels between neigboring control points" << endl;
	cout << "                        in x, y, z directions. If not specified, the optimal spacing" << endl;
    cout << "                        is determined automatically based on the image size." << endl;
    cout << "  -p                    Enable padding of the input images if the image region of the" << endl;
    cout << "                        images is not big enough to fit enough control points." << endl;
    cout << "                        (default: no padding)" << endl;
	cout << "  -I <dir>              Request storage of intermediate images and deformations to" << endl;
    cout << "                        files in the specified directory. (default: none)" << endl;
	cout << "  -f <int>              Method to combine two deformation fields, where" << endl;
	cout << "                        0: addition, 1: composition. (default: 0)" << endl;
	cout << "  -C <int>              Request use of mutual-saliency weighting, where 0: no, >0: yes. (default: 1)" << endl;
	cout << "  -e <int>              Threshold for generating mask from input image A, i.e., the subject image." << endl;
    cout << "                        As the input image is rescaled to the intensity range [0-255], this" << endl;
    cout << "                        threshold value has to be given in this range. (default: 12)" << endl;
	cout << "  -a <file>             Mask image. The specified image file must be defined in the" << endl;
    cout << "                        space of image B, i.e., the template space." << endl;
    cout << "                        If a mask image is specified, option -e is ignored. (default: none)" << endl;
	cout << "  -A                    Cost-function-masking using the input mask image specified by -a option. " << endl;
	cout << "  -i <char>             initial deformation (must be defined in template space)"<< endl;
	cout << "  -u <int>              how much memory to use (from high accuracy to low accuracy, 3-most memory (default); 2-a little less memory; 1-less memory; 0-much less memory). " << endl;
    cout << endl;
	cout << " Parameters for discrete optimization: " << endl;
	cout << "  -F <int>              Request fast approximation of deformation, where 0: no, 1: yes. (default: 0)" << endl;
	cout << "  -k <int>              Maximum number of iterations at each level of the multi-resolution pyramid. (default 10)" << endl;
	cout << "  -s <float>            Convergence factor for global to local search. (default: 0.67) " << endl;
	cout << "  -n <int>              Number of discrete labels, where a larger number means denser sampling" << endl;
	cout << "                        of possible displacement space. (default: 4)" << endl;
	cout << "  -m <int>              Defines the way the cost of pixels projects onto the control points," << endl;
	cout << "                        where 0: nearest neighborhood, 1: trilinear, 2: cubic B-spline. (default: 2) " << endl;
	cout << "  -M <int>              Defines the way the pairwise potential is calculated, where" << endl;
	cout << "                        0: pott's distance, 1: truncated absolute difference," << endl;
    cout << "                        2: truncated quadratic differences. (default: 1)" << endl;
    cout << "  -t <int>              Theshold for maximum cost for the pairwise potential. (default: 1)" << endl;
    cout << "  -w <int>              Request the use of the spatially-varying weight when  projecting" << endl;
	cout << "                        the cost to the control points, where 0: no, 1: yes. (default: 1)" << endl;
	cout << "  -g <float>            Regularization weight controlling deformation aggressiveness. (default: 0.1) " << endl;
	cout << "  -S <int>              Similarity measure on attribute vectors, where 0: vector difference," << endl;
    cout << "                        1: vector correlation coefficent. (default: 0)" << endl;
	cout << endl;
    cout << "Standard arguments:" << endl;
    cout << "  -v                    Increase verbosity of output messages." << endl;
    cout << "  -h                    Print help and exit." << endl;
    cout << endl;
	cout << "Example:" << endl;
    cout << "  " << exec_name << " A.nii.gz B.nii.gz Gabor_ A2B.nii.gz def_A2B.nii.gz -r3 -b7,2 -C1 -g0.2 " << endl;
    cout << endl;
    print_contact();
}



char inputImgNameA[1024],inputImgNameB[1024],outputImgNameA2B[1024];
char featureListPrefix[1024], featureListNameA[1024], featureListNameB[1024];
char outputDeformationFieldName[1024];
char folderNameForIntermediateResults[1024];
char imgAThisLevelName[1024], imgBThisLevelName[1024], imgA2BThisLevelName[1024], dfFieldB2AThisLevelName[1024], maskThisLevelName[1024];
char confidenceMapThisLevelName[1024];
char inputMaskName[1024], initDeformationName[1024];

int i,j,k;
int s,t;
int x,y,z;
int featureIndex;
int iterIndex1, interIndex2;
Ivector3d imageSize;
Ivector3d imageSizeThisLevel;
Ivector3d imageSizePreviousLevel;
int saveIntermediateOrNot;
int numControlPointsX, numControlPointsY, numControlPointsZ, numControlPointsXThisLevel, numControlPointsYThisLevel, numControlPointsZThisLevel;
int numControlPointsXPreviousLevel, numControlPointsYPreviousLevel, numControlPointsZPreviousLevel;
int numFeaturesA, numFeaturesB;

unsigned char ***outputImgA2B, ***maskForControlPoints;
float ***imgA2BThisLevelFloat;
float ***confidenceMapPreviousLevel;
Fvector3d ***dfFieldB2A, ***controlPoints;
Fvector3d ***dfFieldB2APreviousLevel, ***controlPointsThisLevel, ***controlPointsPreviousLevel;
unsigned char ****featureMapA, ****featureMapB;
float ****featureMapA2BInitialized;

// ===========================================================================
// TODO: To improve performance, create local Image instance in each function
//       which uses trilinear interpolation.
// ===========================================================================

// ---------------------------------------------------------------------------
inline float trilinearInterpolation(unsigned char ***Img, float xc, float yc, float zc, Ivector3d imageSize)
{
    Image image(Img, imageSize.x, imageSize.y, imageSize.z, DT_UNSIGNED_CHAR, 1, Image::FORMAT_DRAMMS, false);
    return image.value(xc, yc, zc);
}

// ===========================================================================
// list of sub-functions (details in the end of this file)
// ===========================================================================

void SaveConfidenceMapImage(char* filename,float*** ConfidenceMapImage, unsigned char***mask, int x_size,int y_size, int z_size);
void SaveDeformationField3D(char* filename,Fvector3d ***deformationField, int image_X, int image_Y, int image_Z);
void DownSampleUCImage3D(unsigned char ***inputImg, Ivector3d imageSize, unsigned char ***imgThisLevel, int resolutionLevel);
void convertFloatImage2UCImage(float ***imageFloat, unsigned char ***imageUC, Ivector3d imageSize);
void convertUCImage2FloatImage(unsigned char ***imageUC, float ***imageFloat, Ivector3d imageSize);
void FFD3D(unsigned char ***imageA, unsigned char ***imageB, unsigned char ****featureMapA, unsigned char ****featureMapB, unsigned char ***maskForControlPoints, unsigned char ***maskForImages, float*** confidenceMap, Ivector3d imageSize, int numFeatures, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defB2A, float ***imageA2B, int levelIndex, int numLevels, int numForegroundVoxels);
float Bspline(float u, int level);
float calculateEnergyOnFeatures(unsigned char ****tempFeatureMap, float ****movingFeatureMap, float*** confidenceMap, Ivector3d imageSize, int numFeatures);
void gauss3D(int numPoints1, float std1, int numPoints2, float std2, int numPoints3, float std3, float ***h);
void gauss2D(int numPoints1, float std1, int numPoints2, float std2, float theta, float **h);
float gauss1D(float x, float std);
float searchDisplacementAtThisControlPointBasedOnFeatures(unsigned char ***imageA, unsigned char ***imageB, unsigned char ****featureMapA, unsigned char ****featureMapB, float ****featureMapA2BTemp, float ****featureMapA2BIterative, float ***confidenceMap, Ivector3d imageSize, int numFeatures, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defField, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, int NthIter, int NthLevel, float initEnergyAB, int numForegroundVoxels);
void GenerateFloatImageAndDeformationFieldByFFD(unsigned char ***imageA, unsigned char ***imageB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***imageA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels);
void GenerateFloatImageAndDeformationFieldByFFD2(unsigned char ***imageA, unsigned char ***imageB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***imageA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels, int iterIdx, int AdditionOrComposition, float *weightsX, float *weightsY, float *weightsZ, int method, int inputInitDeformationOrNot, int encourageDiffeomorphismOrNot);
void GenerateDeformationFieldByFFD2(Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***imageA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels, int iterIdx, int AdditionOrComposition, float *weightsX, float *weightsY, float *weightsZ, int method, int inputInitDeformationOrNot);
float GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(unsigned char ****featureMapA, unsigned char ****featureMapB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***confidenceMap, int numFeatures, float ****featureMapA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels);
void GenerateFeaturesAndDeformationFieldByFFD(unsigned char ****featureMapA, unsigned char ****featureMapB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, int numFeatures, float ****featureMapA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels);
float FFDLocalEffect(unsigned char ***imageA, int xcoor, int ycoor, int zcoor, Fvector3d ***controlPoints, Fvector3d ***defField, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, float increment, char orientation, int tempOrNot, int featureIndex, int numFeatures);
void CopyFloatImage3D(float ***imgA, float***imgB, Ivector3d imgSize);
float sign(float x);
void CopyFeatureMapsFloat(float ****featureMapOrig, float ****featureMapCopy, Ivector3d imageSize, int numFeatures);
void CopyDisplacementAtControlPoints3D(Fvector3d ***controlPoints, Fvector3d ***controlPointsBackup, int numControlPointsX, int numControlPointsY, int numcontrolPointsZ);
void CopyDeformationField(Fvector3d ***def, Fvector3d ***defCopy, Ivector3d imageSize);
void CopyConfidenceMap(float ***MSMap, float ***MSMapCopy, Ivector3d imageSize);
void calculateDisplacementIncrementAtControlPoints(Fvector3d ***controlPointsBackup, Fvector3d ***controlPoints, Fvector3d ***controlPointsIncrement, int numControlPointsX, int numControlPointsY, int numControlPointsZ);
void smoothDisplacementAtControlPoints3D(Fvector3d ***controlPoints, Fvector3d ***controlPointsSmoothed, int numControlPointsX, int numControlPointsY, int numControlPointsZ, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int levelIndex, int numLevels, float factor);
void UpdateControlPointsWithSmoothIncrement2(Fvector3d ***controlPointsBackup, Fvector3d ***controlPointsIncrement, Fvector3d ***controlPointsIncrementSmoothed, Fvector3d ***controlPoints, int numControlPointsX, int numControlPointsY, int numControlPointsZ, float alpha, unsigned char ***mask, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int levelIndex, int numLevels);
void UpdateControlPointsWithSmoothIncrement(Fvector3d ***controlPointsBackup, Fvector3d ***controlPoints, Fvector3d ***controlPointsUpdated, int numControlPointsX, int numControlPointsY, int numControlPointsZ, int IterIndex, int maxNumIterInResolution, int levelIndex, int numLevels, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, unsigned char ***mask);
void LinearCombinationOfTwoControlPointsMats(Fvector3d ***controlPointsNew, Fvector3d ***controlPointsOld, Fvector3d ***controlPointsUpdated, float weighting, int numControlPointsX, int numControlPointsY, int numControlPointsZ);
void Convolution3DFloat(float*** input,Ivector3d inputSize, float*** mask,Ivector3d maskSize, float*** output);
void Convolution2DFloat(float** input,Ivector2d inputSize,float** mask,Ivector2d maskSize,float** output);
void smoothDeformationField(Fvector3d ***dfField, Fvector3d ***dfFieldSmoothed, Ivector3d dfSize, int levelIndex, int numLevels);
void UpsampleDisplacementAtControlPoints(Fvector3d ***controlPointsPreviousLevel, Fvector3d ***controlPointsThisLevel, int numControlPointsXPreviousLevel, int numControlPointsYPreviousLevel, int numControlPointsZPreviousLevel, int numControlPointsXThisLevel, int numControlPointsYThisLevel, int numControlPointsZThisLevel, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, unsigned char ***maskForControlPointsThisLevel, int levelIndex, int numLevels);
void UpsampleDeformationField(Fvector3d ***dfFieldB2APreviousLevel, Fvector3d ***dfFieldB2AThisLevel, Ivector3d imageSizePreviousLevel, Ivector3d imageSizeThisLevel, Fvector3d resolutionRatio, int smoothOrNot, int levelIndex, int numLevels, float foregroundRatio);
void UpsampleConfidenceMap(float ***MSMapPreviousLevel, float ***MSMapThisLevel, Ivector3d imageSizePreviousLevel, Ivector3d imageSizeThisLevel);
int GenerateMask(unsigned char ***imageA, unsigned char ***imageB, Ivector3d imageSize, unsigned char ***maskForImages, unsigned char ***maskForControlPoints, int levelIndex, int numLevels, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int foregroundThre);
int DilateImageMaskIntoControlPointMask(unsigned char ***maskForImages, unsigned char ***maskForControlPoints, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int levelIndex, int numLevels);
void saveDisplacementAtControlPoints3D(Fvector3d ***controlPoints, int numControlPointsXThisLevel, int numControlPointsYThisLevel, int numControlPointsZThisLevel);
void SmoothUCImage3D(unsigned char ***image, unsigned char ***imageSmoothed, Ivector3d imageSize);
bool ReadFeaturesFromFeatureList(const char* featureImageListFile,unsigned char ****featureMap, const Image::Region& region);
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(unsigned char ****featureMapA, unsigned char ****featureMapB, float ***confidenceMap, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defField, int numFeatures, float ****featureMapA2BUpdated, float ****featureMapA2BOld, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, float increment, char orientation, int NthLevel, int updateFeaturesOrNot);
float GenerateFeaturesFromDeformationFieldAndCalculateEnergy(unsigned char ****featureMapA, unsigned char ****featureMapB, float ****featureMapA2BFloat, Fvector3d ***defField, float ***confidenceMap, Ivector3d imageSize, int numFeatures);
void GenerateFloatImageFromDeformationField(unsigned char ***imageA, unsigned char ***imageB, float ***imageA2BFloat, Fvector3d ***defField, Ivector3d imageSize);
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(unsigned char ****featureMapA, unsigned char ****featureMapB, float ***confidenceMap, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defField, int numFeatures, float ****featureMapA2BUpdated, float ****featureMapA2BOld, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, float increment, char orientation, int NthLevel);
void CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(float ****featureMapOrig, float ****featureMapCopied, Ivector3d imageSize, int numFeatures, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind);
void CopyDispAtTheBlockInducedByThisControlPoint(Fvector3d ***defFieldOrig, Fvector3d ***defFieldCopied, Ivector3d imageSize, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind);
void initializeConfidenceMapWithOnes(float ***confidenceMap, Ivector3d imageSize);
void UpdateConfidenceMapByInputMask(float ***confidenceMap, unsigned char ***maskForImage, Ivector3d imageSize);
void calculateConfidenceMap(float ****featureMapA2BIterative, unsigned char ****featureMapB, unsigned char ***mask, Ivector3d imageSize, int numFeatures, int levelIndex, int distBetweenControlPointsX, int distBetweenControlPointsY, int SimilarityMeasure, float ***confidenceMap);
void calculateConfidenceMap3D(float ****featureMapA2BIterative, unsigned char ****featureMapB, unsigned char ***mask, Ivector3d imageSize, int numFeatures, int levelIndex, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int SimilarityMeasure, float ***confidenceMap);
float calculateEuclideanDistanceBetweenTwoFeatureVectors(float ****featureMapA2B, unsigned char ****featureMapB, int xA, int yA, int zA, int xB, int yB, int zB, int numFeatures);
float calculateCorrelationCoefficientBetweenTwoFeatureVectors(float ****featureMapA2B, unsigned char ****featureMapB, int xA, int yA, int zA, int xB, int yB, int zB, int numFeatures);
void ApplyDeformationFieldOnImage(unsigned char ***image, Fvector3d ***dfFieldB2A, Ivector3d imageSize, unsigned char ***imageA2B);
void ApplyDeformationFieldOnFeatures(unsigned char ****featureMapA, unsigned char ****featureMapB, int numFeatures, Ivector3d imageSize, Fvector3d ***dfFieldB2A, float ****featureMapA2BInitialized);
void IncorporateInitialDeformation(Fvector3d*** def, Fvector3d*** init_def, Ivector3d defsize, int levelIndex, int numImgLevels);

float discreteOptimization(unsigned char ****SF, unsigned char ****TF, Fvector3d ***defField, int numSamples, Ivector3d imageSize, Fvector3d resolutionRatio, int levelIndex, int method, int distMethod, float threshold, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int numFeatures, int distBasedWeight, float regWeight,  float ***confidenceMap, unsigned char ***mask, Fvector3d ***controlPointDisp, int AdditionOrComposition, int chk, float *weightsX, float *weightsY, float *weightsZ, float label_factor, int iter,  int MaxNumIterInThisResolution, int indexGridLevel, int numGridLevels, int fastApproximationOrNot, int SimilarityMeasure, float initEnergyPrevIter);

// ===========================================================================
// main
// ===========================================================================

int main(int argc,char *argv[])
{
  int numForegroundVoxels;
  float foregroundRatio;
  FILE *fpA, *fpB;
  struct tms startTime,endTime;
  float durTime;
  Fvector3d voxelSize;
  Fvector3d resolutionRatio;
  float minResolution;
  int iter;
  bool sizeok;

  int numImgLevels=3; 
  int distBetweenControlPointsX = 0; // if not specified by user, it is determined below
  int distBetweenControlPointsY = 0; // if not specified by user, it is determined below
  int distBetweenControlPointsZ = 0; // if not specified by user, it is determined below
  saveIntermediateOrNot = NNO;
  int AdditionOrComposition = 0;  // 0 for addition and 1 for composition
  int useConfidenceMapOrNot = 1;  // default: on
  int foregroundThre = 12; // threshold for generating binary mask for control points
  int inputMaskOrNot = NNO;
  int CostFunctionMasking = NNO;
  int fastApproximationOrNot = NNO;
  int encourageDiffeomorphismOrNot = NNO;
  int SimilarityMeasure=SSD; // by default, similarity is measured on difference of vectors
  sprintf(inputMaskName, "NULL");
  int inputInitDeformationOrNot = NNO;
  int useMemory=3;
  
  int numGridLevels, indexGridLevel; // to allow two grid levels in the highest resolution.
  
  int numSamples = 4;
  int numSamplesThisLevel;  //to increase number of samples at lower image resolution
  int method = 2;  //1 for trilinear; 2 for cubic B-splines
  int distMethod = 2;
  float threshold = 0.8; //20;  // threshold for pair-wise distance
  int distBasedWeight = 1;   // strongly suggested to be 1
  float regWeight = 0.1f;	  // default parameters for Aris' discrete optimization
  int MaxNumIterInEachResolution = 10;  // number of iterations in each resolution
  int MaxNumIterInThisResolution;
  int MaxNumIterInThisGridLevelThisResolution;
  float label_factor = 0.67;  // factor to control global-to-local search during iterations in each resolution

  bool enable_padding = false; // whether to allow input images to be padded

  int verbose = 0; // verbosity of output messages

  int c=-1;
  while((c=getopt(argc,argv,"r:b:n:m:M:I:f:t:w:g:u:C:e:a:Ak:s:F:S:D:pi:hv")) != -1)
    {
      switch(c)
		{
			case 'r':
				sscanf(optarg, "%d", &numImgLevels);
				break;
				
			case 'b':
				sscanf(optarg, "%d, %d, %d", &distBetweenControlPointsX,
                                             &distBetweenControlPointsY,
                                             &distBetweenControlPointsZ);
				break;
				
			case 'n':
			    sscanf(optarg, "%d", &numSamples);
				break;
			
			case 'm':
				sscanf(optarg, "%d", &method);
				break;
			
			case 'M':
				sscanf(optarg, "%d", &distMethod);
				break;
		
		    case 'f':
			    sscanf(optarg, "%d", &AdditionOrComposition);
				break;
				
			case 't':
				sscanf(optarg, "%f", &threshold);
				break;
				
			case 'w':
				sscanf(optarg, "%d", &distBasedWeight);
				break;
				
			case 'g':
				sscanf(optarg, "%f", &regWeight);
				break;
			
			case 'u':
				sscanf(optarg, "%d", &useMemory);
				break;
				
			case 'I':
			    saveIntermediateOrNot = YYES;
			    sscanf(optarg, "%s", folderNameForIntermediateResults);
				break;
				
			case 'C':
				sscanf(optarg, "%d",&useConfidenceMapOrNot);
				break;
				
			case 'e':
				sscanf(optarg, "%d", &foregroundThre);
				break;
			
			case 'a':
			    inputMaskOrNot=YYES;
				sscanf(optarg, "%s", inputMaskName);
				break;
			
			case 'A':
				CostFunctionMasking=YYES;
				break;
				
			case 'k':
				sscanf(optarg, "%d", &MaxNumIterInEachResolution);
				break;
				
			case 's':
				sscanf(optarg, "%f", &label_factor);
				break;
			
			case 'F':
				sscanf(optarg, "%d", &fastApproximationOrNot);
				break;
			
			case 'D':
				sscanf(optarg, "%d", &encourageDiffeomorphismOrNot);
				break;
			
			case 'S':
				sscanf(optarg, "%d", &SimilarityMeasure);
				break;
			
			case 'i':
				inputInitDeformationOrNot = YYES;
				sscanf(optarg, "%s", initDeformationName);
				break;

            case 'p':
                enable_padding = true;
                break;

			case 'h':
				print_help();
                exit(0);
				
            case 'v':
                verbose++;
                break;

			default:
                // error message printed by getopt() already
                exit(1);
		}
    }
	
	argc -= optind;
    argv += optind;

    if (argc != 5) {
        print_help();
        exit(1);
    }
	
   if ( (SimilarityMeasure!=SSD)&(SimilarityMeasure!=CC) ) 
		{
		printf("\n\Error: DRAMMS only accepts SSD (-S0, default) or CC (-S1) as similarity measures!\n\n");
		exit(1);
		}

   if (CostFunctionMasking==YYES && inputMaskOrNot==NNO)
		{
		printf("\n\Error: Cost function masking (-A option) requires the input of a binary mask through -a option!\n\n");
		exit(1);
		}
   // ------------------------------------------------------------------------
   // read input images	
   const char* inputImgNameA              = argv[0];
   const char* inputImgNameB              = argv[1];
   const char* featureListPrefix          = argv[2];
   const char* outputImgNameA2B           = argv[3];
   const char* outputDeformationFieldName = argv[4];
   
   Image* inputimageA = ReadImage(inputImgNameA, DT_UNSIGNED_CHAR);
   if (inputimageA == NULL) {
     fprintf(stderr, "Failed to read image A from file %s!\n", inputImgNameA);
     exit(1);
   }
   if (inputimageA->region.nz <= 1) {
     fprintf(stderr, "Image %s has no three dimensions! Use Deform2D to register two-dimensional images.\n", inputImgNameA);
     exit(1);
   }

   Image* inputimageB = ReadImage(inputImgNameB, DT_UNSIGNED_CHAR);
   if (inputimageB == NULL) {
     fprintf(stderr, "Failed to read image B from file %s!\n", inputImgNameB);
     delete inputimageA;
     exit(1);
   }
   if (inputimageB->region.nz <= 1) {
     fprintf(stderr, "Image %s has no three dimensions! Use Deform2D to register two-dimensional images.\n", inputImgNameB);
     delete inputimageA;
     exit(1);
   }

   Image* mask = NULL;
   if (inputMaskOrNot) {
     mask = ReadImage(inputMaskName, DT_UNSIGNED_CHAR);
     if (mask == NULL) {
        fprintf(stderr, "Failed to read mask from file %s!\n", inputMaskName);
        delete inputimageA;
        delete inputimageB;
        exit(1);
     }
     if (mask->region.nx != inputimageB->region.nx ||
            mask->region.ny != inputimageB->region.ny ||
            mask->region.nz != inputimageB->region.nz) {
        fprintf(stderr, "Mask size does not match size of image B!\n");
        delete inputimageA;
        delete inputimageB;
        exit(1);
     }
   } else {
	 mask = new Image(inputimageB->region.nx,
                      inputimageB->region.ny,
                      inputimageB->region.nz,
                      DT_UNSIGNED_CHAR, 1,
                      Image::FORMAT_DRAMMS);
   }
 
   Image* initdeformation = NULL;
   if (inputInitDeformationOrNot) {
	    initdeformation = ReadDeformationField(initDeformationName);
        if (initdeformation == NULL) {
            fprintf(stderr, "Failed to read initial deformation from file %s!\n", initDeformationName);
            delete inputimageA;
            delete inputimageB;
            delete mask;
            exit(1);
        }
	}

   // ------------------------------------------------------------------------
   // determine control point spacing and required image region
   bool spacing_ok = true;
   int maxdist = static_cast<int>(std::max(0.1 * inputimageB->region.nx,
                                           0.1 * inputimageB->region.ny));
   if (distBetweenControlPointsX > maxdist) {
      fprintf(stderr, "Invalid control point spacing! Please choose a value smaller than %d for x.\n", maxdist);
      spacing_ok = false;
   }
   if (distBetweenControlPointsY > maxdist) {
      fprintf(stderr, "Invalid control point spacing! Please choose a value smaller than %d for y.\n", maxdist);
      spacing_ok = false;
   }
   if (distBetweenControlPointsZ > static_cast<int>(inputimageB->region.nz)) {
      fprintf(stderr, "Invalid control point spacing! Please choose a value smaller than %d for z.\n", static_cast<int>(inputimageB->region.nz));
      spacing_ok = false;
   }
   if (!spacing_ok) {
      delete inputimageA;
      delete inputimageB;
      delete mask;
      if (initdeformation) delete initdeformation;
      exit(1);
   }

   Image::Region required_region;
   if (!enable_padding) required_region = inputimageB->region;

   {
      if (verbose > 0) {
         printf("\nDetermine optimal control point spacing and required image region...\n");
         fflush(stdout);
      }

      Image* tmpmask = GenerateMask(inputimageB, 25);
      if (tmpmask == NULL) {
         fprintf(stderr, "Failed to generate temporary mask!\n");
         delete inputimageA;
         delete inputimageB;
         delete mask;
         if (initdeformation) delete initdeformation;
         exit(1);
      }
      RemoveNoiseFromMask(tmpmask);
	  bool useLessMemory=false;
	  if (useMemory==1) useLessMemory=true;
      DetermineControlPointSpacing(tmpmask, distBetweenControlPointsX,
                                            distBetweenControlPointsY,
                                            distBetweenControlPointsZ,
                                            required_region, useLessMemory);
      delete tmpmask;

      if (verbose > 1) {
         printf("\nControl point spacing: %d,%d,%d\n", distBetweenControlPointsX,
                                                     distBetweenControlPointsY,
                                                     distBetweenControlPointsZ);
         printf("Input image region:    offset = (%d,%d,%d), size = (%d,%d,%d)\n",
                 inputimageB->region.ox, inputimageB->region.oy, inputimageB->region.oz,
                 inputimageB->region.nx, inputimageB->region.ny, inputimageB->region.nz);
         printf("Required image region: offset = (%d,%d,%d), size = (%d,%d,%d)\n",
                 required_region.ox, required_region.oy, required_region.oz,
                 required_region.nx, required_region.ny, required_region.nz);
         fflush(stdout);
      }
   }

   // ------------------------------------------------------------------------
   // pad input images if required
   if (required_region != inputimageB->region) {
        bool   ok  = true;
        Image* tmp = NULL;

        if (verbose > 0) {
            printf("\nPadding input images such that image domain covers required region...\n");
            fflush(stdout);
        }

        // image A
        tmp = ResizeImage(inputimageA, required_region, 0, true);
        if (tmp == NULL) ok = false;
        if (ok) {
            delete inputimageA;
            inputimageA = tmp;
        }
        // image B
        if (ok) {
            tmp = ResizeImage(inputimageB, required_region, 0, true);
            if (tmp == NULL) ok = false;
            if (ok) {
                delete inputimageB;
                inputimageB = tmp;
            }
        }
        // mask
        if (ok) {
            tmp = ResizeImage(mask, required_region, 0, true);
            if (tmp == NULL) ok = false;
            if (ok) {
                delete mask;
                mask = tmp;
            }
        }
        // initial deformation
        if (ok && initdeformation) {
            tmp = ResizeImage(initdeformation, required_region, 0, true);
            if (tmp == NULL) ok = false;
            if (ok) {
                delete initdeformation;
                initdeformation = tmp;
            }
        }
        // bail out on error
        if (!ok) {
            fprintf(stderr, "Failed to allocate memory!\n");
            delete inputimageA;
            delete inputimageB;
            delete mask;
            if (initdeformation) delete initdeformation;
            exit(1);
        }
   }
 
   // ------------------------------------------------------------------------
   // common attributes
   imageSize.x = inputimageB->region.nx;
   imageSize.y = inputimageB->region.ny;
   imageSize.z = inputimageB->region.nz;
   voxelSize.x = inputimageB->hdr.pixdim[1];
   voxelSize.y = inputimageB->hdr.pixdim[2];
   voxelSize.z = inputimageB->hdr.pixdim[3];

   numControlPointsX = static_cast<int>(ceil((float)imageSize.x/(float)distBetweenControlPointsX));
   numControlPointsY = static_cast<int>(ceil((float)imageSize.y/(float)distBetweenControlPointsY));
   numControlPointsZ = static_cast<int>(ceil((float)imageSize.z/(float)distBetweenControlPointsZ));
   controlPoints = Fvector3dalloc3d(numControlPointsX, numControlPointsY, numControlPointsZ);
   
   minResolution = MIN(MIN(voxelSize.x, voxelSize.y), voxelSize.z);
   resolutionRatio.x = voxelSize.x/minResolution;
   resolutionRatio.y = voxelSize.y/minResolution;
   resolutionRatio.z = voxelSize.z/minResolution;
   printf("\ninput image A = %s \ninput image B = %s\n", inputImgNameA, inputImgNameB);
   printf("\nfeature list prefix = %s\n", featureListPrefix);
   printf("\nimage size = (%d, %d, %d)\n", imageSize.x, imageSize.y, imageSize.z);
   printf("\nvoxel size = (%f, %f, %f)\n", voxelSize.x, voxelSize.y, voxelSize.z);
   printf("\nresolution ratio = (%f, %f, %f)\n", resolutionRatio.x, resolutionRatio.y, resolutionRatio.z);
   printf("\nsave intermediate results? (1-yes; 0-no) - %d", saveIntermediateOrNot);
   printf("\nHow to combine two deformation fields (0-addition, 1-composition): %d\n ", AdditionOrComposition);
   printf("\nUse mutual-saliency map weighting or not (0-no,1-yes)? : %d\n", useConfidenceMapOrNot);
   printf("inputInitDeformationOrNot (0-no, 1-yes)? : %d\n", inputInitDeformationOrNot); 
   if (saveIntermediateOrNot==YYES)
     printf(", into folder %s\n", folderNameForIntermediateResults);
   else
     printf("\n");
   printf("\nmulti-resolution: number of levels = %d\n\nnumber of voxels between two adjacent control points: x-%d, y-%d, z-%d (in finest resolution)\n", numImgLevels, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ);
   printf("\nencourage diffeomorphism or not = %d (0-no, 1-yes)\n", encourageDiffeomorphismOrNot);

   // ------------------------------------------------------------------------
   // multi-resolution implementation   
   times(&startTime);

   for (i = 0; i < numImgLevels; i++) {
	 int levelIndex = numImgLevels-i;
     int levelRatio = static_cast<int>(pow(2, levelIndex - 1));

	 numGridLevels=1;
	 //if ( i>1 && imageSize.x>220 && imageSize.y>220 && imageSize.z>150 &&regWeight<=5.0)
	 if ( i>1 && MAX(imageSize.x*resolutionRatio.x, imageSize.y*resolutionRatio.y)>150 && (imageSize.z*resolutionRatio.z)>150 && regWeight<=5.5 && inputInitDeformationOrNot==NNO )   // Yangming changed on 03/30/2012 to solve registration errors in adni data, need to check if this works in LONI dataset
	 	numGridLevels=2;       // allow 2 grid levels in the highest resolution if image is big and requirement for accuracy is high
	 if ( useMemory == 0 )
		numGridLevels=1;
		
     imageSizeThisLevel.x = imageSize.x / levelRatio;
     imageSizeThisLevel.y = imageSize.y / levelRatio;
     imageSizeThisLevel.z = imageSize.z / levelRatio;

	 printf("\033[31m\n---------------------------------------------------------------\n\033[m");
	 printf("\033[31m\nSearch in resolution level %d, image size = (%d, %d, %d)\n\033[m", levelIndex, imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);
	 printf("\033[31m\n---------------------------------------------------------------\n\033[m");
	 
	 // Step 1-1: downsample input images at current level
	 printf("step 1-1: downsample input images to 1/%d of their size...\n", levelRatio);

	 Image* imageAthislevel = DownsampleImage(inputimageA, levelRatio);
  	 Image* imageBthislevel = DownsampleImage(inputimageB, levelRatio);
	 Image* maskthislevel   = DownsampleImage(mask,        levelRatio);

     Image* initdeformationthislevel = NULL;
     if (initdeformation != NULL) {
        initdeformationthislevel = DownsampleImage(initdeformation, levelRatio);
     }

     assert(imageBthislevel->region.nx == imageSizeThisLevel.x);
     assert(imageBthislevel->region.ny == imageSizeThisLevel.y);
     assert(imageBthislevel->region.nz == imageSizeThisLevel.z);

	 int distBetweenControlPointsXThisLevel = distBetweenControlPointsX;
	 int distBetweenControlPointsYThisLevel = distBetweenControlPointsY;
	 int distBetweenControlPointsZThisLevel = distBetweenControlPointsZ;

	 // Step 1-2: at each level, allocate memory for the images and deformations
	 printf("step 1-2: allocate required memory for this level...\n");
	 Image* deformationthislevel          = new Image(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z, DT_FLOAT, 3, Image::FORMAT_DRAMMS);
	 Image* confidencemapthislevel        = new Image(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z, DT_FLOAT);
	 Image* imageA2Bthisleveluc           = new Image(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z, DT_UNSIGNED_CHAR);
	 Image* maskforcontrolpointsthislevel = new Image(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z, DT_UNSIGNED_CHAR);

     deformationthislevel         ->CopyRegion(imageBthislevel); // important for cropping as implicitly done by WriteImage()
     confidencemapthislevel       ->CopyRegion(imageBthislevel);
     imageA2Bthisleveluc          ->CopyRegion(imageBthislevel);
     maskforcontrolpointsthislevel->CopyRegion(imageBthislevel);

     deformationthislevel         ->CopyTransform(imageBthislevel); // defines image to be in space B
     confidencemapthislevel       ->CopyTransform(imageBthislevel);
     imageA2Bthisleveluc          ->CopyTransform(imageBthislevel);
     maskforcontrolpointsthislevel->CopyTransform(imageBthislevel);

	 unsigned char*** imgAThisLevel 			    = imageAthislevel->img.uc;
	 unsigned char*** imgBThisLevel 			    = imageBthislevel->img.uc;
	 unsigned char*** maskForImagesThisLevel        = maskthislevel->img.uc;
     Fvector3d***     initDeformationThisLevel      = initdeformationthislevel ? initdeformationthislevel->img.v3 : NULL;
	 float***		  confidenceMapThisLevel 	    = confidencemapthislevel->img.fl;
	 Fvector3d***     dfFieldB2AThisLevel    	    = deformationthislevel->img.v3;
	 unsigned char*** imgA2BThisLevelUC			    = imageA2Bthisleveluc->img.uc;
	 unsigned char*** maskForControlPointsThisLevel = maskforcontrolpointsthislevel->img.uc;

	 imgA2BThisLevelFloat = Falloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z); 

	 // mask in the space of image B
	 printf("dim1=(%d, %d, %d), dim2=(%d, %d, %d)",
            maskthislevel->region.nx,
            maskthislevel->region.ny,
            maskthislevel->region.nz,
            maskforcontrolpointsthislevel->region.nx,
            maskforcontrolpointsthislevel->region.ny,
            maskforcontrolpointsthislevel->region.nz);

	 if (inputMaskOrNot)
		{
		printf("input mask %s read. \nDilate it into mask for control points...", inputMaskName);
		numForegroundVoxels = DilateImageMaskIntoControlPointMask(maskForImagesThisLevel,maskForControlPointsThisLevel, imageSizeThisLevel, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, levelIndex, numImgLevels);
		printf("done!\n");
		foregroundRatio = (float)numForegroundVoxels/(float)(imageSizeThisLevel.x*imageSizeThisLevel.y*imageSizeThisLevel.z);
		}
	 else
		{
		printf("\t generate mask for control points \n");
		printf("\t foreground threshold = %d\n", foregroundThre);
		numForegroundVoxels = GenerateMask(imgAThisLevel, imgBThisLevel, imageSizeThisLevel, maskForImagesThisLevel, maskForControlPointsThisLevel, levelIndex, numImgLevels, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, foregroundThre);
		foregroundRatio = (float)numForegroundVoxels/(float)(imageSizeThisLevel.x*imageSizeThisLevel.y*imageSizeThisLevel.z);
	 	}

	 // Step 2: read features for both image under registration
	 printf("Step 2: Read attribute images at this resolution...\n");
	 sprintf(featureListNameA, "%sA_level%d.lst", featureListPrefix, levelIndex);
	 sprintf(featureListNameB, "%sB_level%d.lst", featureListPrefix, levelIndex);
	 printf("\nfeature list name A = %s\n", featureListNameA);
	 printf("feature list name B = %s\n", featureListNameB);
	 
	 if (NULL==(fpA=fopen(featureListNameA,"rb"))){
       printf("File %s doesn't exist!\n",featureListNameA);
	   exit(1);
	   }
	 if (NULL==(fpB=fopen(featureListNameB,"rb"))){
       printf("File %s doesn't exist!\n",featureListNameB);
	   exit(1);
	   }
	 fscanf(fpA,"%d",&numFeaturesA);  
	 fscanf(fpB,"%d",&numFeaturesB);
	 if (numFeaturesA != numFeaturesB) {
	   printf("number of features for image A and B must match! (numFeatureA = %d, numFeaturesB = %d\n)", numFeaturesA, numFeaturesB);
	   exit(1);
	   }
	 else
	   printf("numFeatureA = %d\nnumFeaturesB = %d\n", numFeaturesA, numFeaturesB);
	   
	 printf("\nReading feature images ...\n");
	 featureMapA = (unsigned char****)malloc(sizeof(unsigned char***)*numFeaturesA);	   
	 featureMapB = (unsigned char****)malloc(sizeof(unsigned char***)*numFeaturesB);	
	 #pragma omp parallel for shared(numFeaturesA, imageSizeThisLevel) private(featureIndex) num_threads(numFeaturesA)
     for (featureIndex=0;featureIndex<numFeaturesA;featureIndex++)
	   {
	   featureMapA[featureIndex] = UCalloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);	 
	   featureMapB[featureIndex] = UCalloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);	 
	   }
	   
	 //omp_set_num_threads(2);
	 #pragma omp parallel sections
	 {
	 #pragma omp section
	 {
	 ReadFeaturesFromFeatureList(featureListNameA, featureMapA, imageBthislevel->region);
	 }
	 #pragma omp section
	 {
	 ReadFeaturesFromFeatureList(featureListNameB, featureMapB, imageBthislevel->region);
	 }
	 }
		 
	 // step 3: initialize displacement at control points
	 for (indexGridLevel=1; indexGridLevel<=numGridLevels; indexGridLevel++)
	 { //for (indexGridLevel=1; indexGridLevel<=numGridLevels; indexGridLevel++)

	 if ( (numGridLevels==2)&&(indexGridLevel==1) )
		{
		distBetweenControlPointsXThisLevel = MAX((int)((float)distBetweenControlPointsX*pow(1.15, 3-levelIndex)+0.5), MINCPDISTXY);
		distBetweenControlPointsYThisLevel = MAX((int)((float)distBetweenControlPointsY*pow(1.15, 3-levelIndex)+0.5), MINCPDISTXY);
		distBetweenControlPointsZThisLevel = MAX((int)((float)distBetweenControlPointsZ*pow(1.15, 3-levelIndex)+0.5), MINCPDISTZ);
		}
	 else 
		{
		distBetweenControlPointsXThisLevel = distBetweenControlPointsX;
		distBetweenControlPointsYThisLevel = distBetweenControlPointsY;
		distBetweenControlPointsZThisLevel = distBetweenControlPointsZ;
		
		if (levelIndex==3)   // lowest image resolution
			{
			if (foregroundRatio>0.40)
				{
				distBetweenControlPointsXThisLevel = (int)(distBetweenControlPointsX*2.0+0.5);
				distBetweenControlPointsYThisLevel = (int)(distBetweenControlPointsY*2.0+0.5);
				distBetweenControlPointsZThisLevel = (int)(distBetweenControlPointsZ*2.0+0.5);
				}
			else if (foregroundRatio>0.25)
				{
				distBetweenControlPointsXThisLevel = (int)(distBetweenControlPointsX*1.5+0.5);
				distBetweenControlPointsYThisLevel = (int)(distBetweenControlPointsY*1.5+0.5);
				distBetweenControlPointsZThisLevel = (int)(distBetweenControlPointsZ*1.5+0.5);
				}
			}
		if (levelIndex==2) // middle image resolution
			{
			if (foregroundRatio>0.25)
				{
				distBetweenControlPointsXThisLevel = distBetweenControlPointsX+1;
				distBetweenControlPointsYThisLevel = distBetweenControlPointsY+1;
				distBetweenControlPointsZThisLevel = distBetweenControlPointsZ+1;
				}
			else
				{
				distBetweenControlPointsXThisLevel = distBetweenControlPointsX;
				distBetweenControlPointsYThisLevel = distBetweenControlPointsY;
				distBetweenControlPointsZThisLevel = distBetweenControlPointsZ;
				}
			}
		if (levelIndex==1)  // highest image resolution
			{
			distBetweenControlPointsXThisLevel = MAX(MINCPDISTXY, distBetweenControlPointsXThisLevel-1);
			distBetweenControlPointsYThisLevel = MAX(MINCPDISTXY, distBetweenControlPointsYThisLevel-1);
			if (distBetweenControlPointsZThisLevel>2)
				distBetweenControlPointsZThisLevel = MAX(MINCPDISTZ, distBetweenControlPointsZThisLevel-1);
			}
		}
	 
	 
	 
	 printf("\n\n\n");
	 printf("Step 3: initialize (sparse) displacement at control points in this resolution\n");
	 printf("        grid level %d of %d\n", indexGridLevel, numGridLevels);
	 
	 while(true) {
	 numControlPointsXThisLevel = (int)ceil((float)imageSizeThisLevel.x/(float)distBetweenControlPointsXThisLevel);
	 numControlPointsYThisLevel = (int)ceil((float)imageSizeThisLevel.y/(float)distBetweenControlPointsYThisLevel);
	 numControlPointsZThisLevel = (int)ceil((float)(imageSizeThisLevel.z-1)/(float)distBetweenControlPointsZThisLevel);
	 
	 sizeok=true;
	 if ( (useMemory==2)&&(numControlPointsXThisLevel>80) ) {
		distBetweenControlPointsXThisLevel++;
		sizeok=false;
		}
	 if ( (useMemory==2)&&(numControlPointsYThisLevel>80) ) {
		distBetweenControlPointsYThisLevel++;
		sizeok=false;
		}
	 if ( (useMemory==2)&&(numControlPointsZThisLevel>100) ) {
		distBetweenControlPointsZThisLevel++;
		sizeok=false;
		}
	 if ( (useMemory==1)&&(numControlPointsXThisLevel>70) ) {
		distBetweenControlPointsXThisLevel++;
		sizeok=false;
		}
	 if ( (useMemory==1)&&(numControlPointsYThisLevel>70) ) {
		distBetweenControlPointsYThisLevel++;
		sizeok=false;
		}
	 if ( (useMemory==1)&&(numControlPointsZThisLevel>90) ) {
		distBetweenControlPointsZThisLevel++;
		sizeok=false;
		}
		
	 if(sizeok)
		break;
	 }
	 printf(" \n*** number of control points at this level = %d (%d*%d*%d)\n\n", (numControlPointsXThisLevel*numControlPointsYThisLevel*numControlPointsZThisLevel), numControlPointsXThisLevel, numControlPointsYThisLevel, numControlPointsZThisLevel);
	 controlPointsThisLevel = Fvector3dalloc3d(numControlPointsXThisLevel, numControlPointsYThisLevel, numControlPointsZThisLevel);
	 
	 // for levels other than the coarsest level, take the results in the previous level as initilization
	 if ( (numImgLevels>1) && (i!=0) && (indexGridLevel==1) )
	   {
	   //1. Initialize the dense deformation 
	   UpsampleDeformationField(dfFieldB2APreviousLevel, dfFieldB2AThisLevel, imageSizePreviousLevel, imageSizeThisLevel, resolutionRatio, NNO, levelIndex, numImgLevels, foregroundRatio);
	   Fvector3dfree3d(dfFieldB2APreviousLevel, imageSizePreviousLevel.z, imageSizePreviousLevel.x);
	   if (inputInitDeformationOrNot){
		  printf("incorporate initial deformation...\n");
		  IncorporateInitialDeformation(dfFieldB2AThisLevel, initDeformationThisLevel, imageSizeThisLevel, levelIndex, numImgLevels);
		  }
	   
	   // 2. initialize the confidence map
	   printf("initial mutual-saliency map...\n");
	   if (useConfidenceMapOrNot>0)
	     {
			if (levelIndex<numImgLevels) // in the high and middle resolutions, upsample the mutual-saliency maps calculated from previous resolution, to balance between accuracy and computation.
			{
			UpsampleConfidenceMap(confidenceMapPreviousLevel, confidenceMapThisLevel, imageSizePreviousLevel, imageSizeThisLevel);
			Ffree3d(confidenceMapPreviousLevel, imageSizePreviousLevel.z, imageSizePreviousLevel.x);
			}
		 } // if use confidence map
	   else
			initializeConfidenceMapWithOnes(confidenceMapThisLevel, imageSizeThisLevel);
	   // Ou added on 05/09/2009 (end)
	   }
	 else   // for the coarsest level
	   {
	   // 1. initialize the displacements at control points
	   for (z=0;z<numControlPointsZThisLevel;z++)
	    for (x=0;x<numControlPointsXThisLevel;x++)
	     for (y=0;y<numControlPointsYThisLevel;y++)
		   {
		   controlPointsThisLevel[z][x][y].x=0;
		   controlPointsThisLevel[z][x][y].y=0;
		   controlPointsThisLevel[z][x][y].z=0;
		   }
	   // 2. initialize the confidence map
	   initializeConfidenceMapWithOnes(confidenceMapThisLevel, imageSizeThisLevel);
	   
	   // 3. initialize the dense deformation if there is input initial deformation
	   if (inputInitDeformationOrNot) {
		  printf("incorporate initial deformation...\n");
		  IncorporateInitialDeformation(dfFieldB2AThisLevel, initDeformationThisLevel, imageSizeThisLevel, levelIndex, numImgLevels);
		  }
	   }
	  
	  
	 //////////////////
	 // Yangming added on June 28, 2013, to correctly use the input mask for cost-function-masking (do not utilize image/feature information in the background of input mask if there is input mask given)
	 /////////////////
	 if (CostFunctionMasking)
		UpdateConfidenceMapByInputMask(confidenceMapThisLevel, maskForImagesThisLevel, imageSizeThisLevel);
	 /////////////////
	 
	 
     if ( (useConfidenceMapOrNot>0)&&(indexGridLevel==1)&&(saveIntermediateOrNot==YYES) )
		{
		sprintf(confidenceMapThisLevelName, "%s/MutualSaliencyMap_level%d.nii.gz", folderNameForIntermediateResults, levelIndex);
		printf("\n** save confidence map into %s **\n\n", confidenceMapThisLevelName);
		WriteImage(confidenceMapThisLevelName, confidencemapthislevel, maskthislevel);
		}
	 
	 
	 // step 4: save initilized transformed image
	 if ( (saveIntermediateOrNot==YYES)&&(indexGridLevel==1) )
	   {
	   printf("Step 4: Save initialized results:\n");
	   #pragma omp parallel sections num_threads(2)
	   {
	   #pragma omp section
	   {
	   sprintf(imgA2BThisLevelName, "%s/A2B_level%d_init.nii.gz", folderNameForIntermediateResults, levelIndex);
	   ApplyDeformationFieldOnImage(imgAThisLevel, dfFieldB2AThisLevel, imageSizeThisLevel, imgA2BThisLevelUC);
	   WriteImage(imgA2BThisLevelName, imageA2Bthisleveluc);
	   printf("\t%s\n", imgA2BThisLevelName);
	   }
	   
	   // save initial deformation field
	   #pragma omp section
	   {
	   printf("Saving initial deformation field in level %d\n", levelIndex);
	   sprintf(dfFieldB2AThisLevelName, "%s/DField_level%d_init.nii.gz", folderNameForIntermediateResults, levelIndex);
	   WriteImage(dfFieldB2AThisLevelName, deformationthislevel);
	   printf("\t%s\n", dfFieldB2AThisLevelName);
	   }
	   }
	   }
	 
	 // Step 5: Optimization:  search displacements at control points in this level, then use the displacement at control points to generate warped image and deformation field
	 printf("\nStep 5: Now start calculating displacement at control points:\n");

	 // pre-computed FFD weights
	 float *weightsXThisLevel, *weightsYThisLevel, *weightsZThisLevel;
	 float *weightsXThisLevel_BSpline, *weightsYThisLevel_BSpline, *weightsZThisLevel_BSpline;
	 printf("weighting scheme = %d (0 for nearest neighbor; 1 for trilinear, 2 for cubic B-spline\n", method);
	 printf("distBasedWeight=%d\n",distBasedWeight);
	 
	 //check if a weighting scheme based on the distance from the nodes will be used and in that case precompute the weights
	 int  chk=0;
	 if(distBasedWeight == 1)
	 {
		//nearest neighbor weighting scheme
		if(method == 0)
		{
			weightsXThisLevel = NULL; weightsYThisLevel = NULL; weightsZThisLevel = NULL;
			chk=1;
		}
		//linear weighting scheme
		else if(method == 1)   
		{
			weightsXThisLevel = (float *)malloc(2*imageSizeThisLevel.x*sizeof(float));
			weightsYThisLevel = (float *)malloc(2*imageSizeThisLevel.y*sizeof(float));
			weightsZThisLevel = (float *)malloc(2*imageSizeThisLevel.z*sizeof(float));		
			precomputeWeights(imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, resolutionRatio, method, weightsXThisLevel, weightsYThisLevel, weightsZThisLevel);
			
			weightsXThisLevel_BSpline = (float *)malloc(4*imageSizeThisLevel.x*sizeof(float));
			weightsYThisLevel_BSpline = (float *)malloc(4*imageSizeThisLevel.y*sizeof(float));
			weightsZThisLevel_BSpline = (float *)malloc(4*imageSizeThisLevel.z*sizeof(float));
			precomputeWeights(imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, resolutionRatio, 2, weightsXThisLevel_BSpline, weightsYThisLevel_BSpline, weightsZThisLevel_BSpline);
			chk=1;
		}
		// cubic B spline weighting scheme
		else if(method == 2)
		{
			weightsXThisLevel = (float *)malloc(4*imageSizeThisLevel.x*sizeof(float));
			weightsYThisLevel = (float *)malloc(4*imageSizeThisLevel.y*sizeof(float));
			weightsZThisLevel = (float *)malloc(4*imageSizeThisLevel.z*sizeof(float));
			precomputeWeights(imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, resolutionRatio, method, weightsXThisLevel, weightsYThisLevel, weightsZThisLevel);
			chk=1;
		}
	 }

	 // hierarchically increase the number of samples in the lower resolution, This will provides better initilization in the higher level
	 numSamplesThisLevel = numSamples;
	 
     // input tentative deformation field "dfFieldB2AThisLevel", output displacement at all control points "controlPointsThisLevel"
     printf("numSamples = %d at this level\n", numSamplesThisLevel);
     printf("imageSize = (%d, %d, %d)\n", imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);
     printf("method = %d\n", method);
     printf("distMethod = %d\n", distMethod);
	 printf("threshold = %f\n", threshold);
     printf("dist between two control points = (%d,%d,%d)\n", distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel);
     printf("num of features = %d\n", numFeaturesA);
     printf("distBasedWeight = %d\n", distBasedWeight);
	 printf("reg weight = %f\n", regWeight);
	 
	 MaxNumIterInThisResolution = MAX(MaxNumIterInEachResolution - (levelIndex-1), 1);
	 
	 if (numGridLevels==1)
		MaxNumIterInThisGridLevelThisResolution = MaxNumIterInThisResolution;
	 else
		MaxNumIterInThisGridLevelThisResolution = MAX((int)(4*MaxNumIterInThisResolution/(numGridLevels*3)+0.5), 1);
	 printf("num iterations in this level of grid in this resoltion = %d\n", MaxNumIterInThisGridLevelThisResolution);
	 printf("label_factor = %f\n", label_factor);
	 
	 
	
	 printf("Start discrete optimization in this level...\n");
     
	 float initEnergyPrevIter = 9999999999.0;
	 float initEnergyThisIter = 0.0;
	 iter=0;
	 while(iter<MaxNumIterInThisGridLevelThisResolution)
	 {
	 printf("************\n");
	 printf("iteration #%d (maximum %d iterations)\n", iter+1, MaxNumIterInThisGridLevelThisResolution);
	 printf("************\n");

	 initEnergyThisIter = discreteOptimization(featureMapA, featureMapB, dfFieldB2AThisLevel, numSamplesThisLevel, imageSizeThisLevel, resolutionRatio, levelIndex, method, distMethod, threshold, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, numFeaturesA, distBasedWeight, regWeight*pow(60.0, (numGridLevels-indexGridLevel)), confidenceMapThisLevel, maskForControlPointsThisLevel, controlPointsThisLevel, AdditionOrComposition, chk, weightsXThisLevel, weightsYThisLevel, weightsZThisLevel, label_factor, iter,  MaxNumIterInThisGridLevelThisResolution, indexGridLevel, numGridLevels, fastApproximationOrNot, SimilarityMeasure, initEnergyPrevIter);

	 // update deformation in this iteration
	 if ( (iter<(MaxNumIterInThisGridLevelThisResolution-1))&&(initEnergyThisIter<initEnergyPrevIter) )
	    {
		if (method==1) 
			GenerateDeformationFieldByFFD2(imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, controlPointsThisLevel, imgA2BThisLevelFloat, dfFieldB2AThisLevel, levelIndex, numImgLevels, iter, AdditionOrComposition, weightsXThisLevel_BSpline, weightsYThisLevel_BSpline, weightsZThisLevel_BSpline, 2, inputInitDeformationOrNot);	 
		else
			GenerateDeformationFieldByFFD2(imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, controlPointsThisLevel, imgA2BThisLevelFloat, dfFieldB2AThisLevel, levelIndex, numImgLevels, iter, AdditionOrComposition, weightsXThisLevel, weightsYThisLevel, weightsZThisLevel, method, inputInitDeformationOrNot);	
		initEnergyPrevIter = initEnergyThisIter;
		iter++;
		}
	 else
		{
		iter = MaxNumIterInThisGridLevelThisResolution;
		// step 6: generate dense deformation field and warp image by the "optimal" displacement at all control points obtained in step 5
		printf("\nStep 6: Generating dense deformation and warp image by optimal displacements at control points\n\n");
		if (method==1)
			GenerateFloatImageAndDeformationFieldByFFD2(imgAThisLevel, imgBThisLevel, imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, controlPointsThisLevel, imgA2BThisLevelFloat, dfFieldB2AThisLevel, levelIndex, numImgLevels, iter, AdditionOrComposition, weightsXThisLevel_BSpline, weightsYThisLevel_BSpline, weightsZThisLevel_BSpline, 2, inputInitDeformationOrNot, encourageDiffeomorphismOrNot);
		else
			GenerateFloatImageAndDeformationFieldByFFD2(imgAThisLevel, imgBThisLevel, imageSizeThisLevel, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, distBetweenControlPointsZThisLevel, controlPointsThisLevel, imgA2BThisLevelFloat, dfFieldB2AThisLevel, levelIndex, numImgLevels, iter, AdditionOrComposition, weightsXThisLevel, weightsYThisLevel, weightsZThisLevel, method, inputInitDeformationOrNot, encourageDiffeomorphismOrNot);
		}
	 } // while
	 
	 Fvector3dfree3d(controlPointsThisLevel, numControlPointsZThisLevel, numControlPointsXThisLevel);
	   if (chk==1)
		{
		free(weightsXThisLevel);
		free(weightsYThisLevel);
		free(weightsZThisLevel);
		
		if (method==1)
			{
			free(weightsXThisLevel_BSpline);
			free(weightsYThisLevel_BSpline);
			free(weightsZThisLevel_BSpline);
			}
		}
		
     }  //for (indexGridLevel=1; indexGridLevel<=numGridLevels; indexGridLevel++)
	 // till here, finished deformation at grid level #indexGridLevel, in resolution #levelIndex

	 // save intermediate results in this level
	 if ( saveIntermediateOrNot==YYES )
	   {
	    printf("saving intermediate warped images and deformation fields into folder %s\n", folderNameForIntermediateResults);
	    sprintf(imgAThisLevelName, "%s/A_level%d.nii.gz", folderNameForIntermediateResults, levelIndex); printf("%s\n", imgAThisLevelName);
		sprintf(imgBThisLevelName, "%s/B_level%d.nii.gz", folderNameForIntermediateResults, levelIndex); printf("%s\n", imgBThisLevelName);
		sprintf(imgA2BThisLevelName, "%s/A2B_level%d.nii.gz", folderNameForIntermediateResults, levelIndex); printf("%s\n", imgA2BThisLevelName);
		sprintf(dfFieldB2AThisLevelName, "%s/DField_level%d.nii.gz", folderNameForIntermediateResults, levelIndex); printf("%s\n", dfFieldB2AThisLevelName);
		sprintf(maskThisLevelName, "%s/mask_level%d.nii.gz", folderNameForIntermediateResults, levelIndex); printf("%s\n", maskThisLevelName);

		WriteImage(imgAThisLevelName, imageAthislevel);
		WriteImage(imgBThisLevelName, imageBthislevel);
		convertFloatImage2UCImage(imgA2BThisLevelFloat, imgA2BThisLevelUC, imageSizeThisLevel);
		WriteImage(imgA2BThisLevelName, imageA2Bthisleveluc);
		WriteImage(dfFieldB2AThisLevelName, deformationthislevel);
		WriteImage(maskThisLevelName, maskforcontrolpointsthislevel);
		}//if (saveIntermediateOrNot==YYES)
	 
	 
	 if ( i==(numImgLevels-1) ) // at the highest resolution
	   {
		printf("\n\nAt the finest resolution,  save final results: \n");
		if (saveIntermediateOrNot==NNO)
			convertFloatImage2UCImage(imgA2BThisLevelFloat, imgA2BThisLevelUC, imageSizeThisLevel);
		WriteImage(outputImgNameA2B, imageA2Bthisleveluc);
		printf("... %s ", outputImgNameA2B);
		WriteImage(outputDeformationFieldName, deformationthislevel);
		printf(" ... %s !\n\n", outputDeformationFieldName);
	   }
	 
	 
	 // 1. save displacements of this control points in this level, as the initialization for displacement at control points of the next level
	 // 2. save confidence map in this level, as the initialization for confidence map of the next level
	 if ( (i!=(numImgLevels-1)) & (numImgLevels>1) ) // coarse levels in the multi-level implementation
	   {
	   dfFieldB2APreviousLevel = Fvector3dalloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);
	   CopyDeformationField(dfFieldB2AThisLevel, dfFieldB2APreviousLevel, imageSizeThisLevel);
	   imageSizePreviousLevel.x = imageSizeThisLevel.x;
	   imageSizePreviousLevel.y = imageSizeThisLevel.y;
	   imageSizePreviousLevel.z = imageSizeThisLevel.z;
	  
	  
	   // calculate confidence maps to be used in the next resolution (level)
	   if (useConfidenceMapOrNot>0)
		  {
			featureMapA2BInitialized = (float****)malloc(sizeof(float***)*numFeaturesA);	
			#pragma omp parallel for shared(featureMapA2BInitialized) private(featureIndex) num_threads(100)
			for (featureIndex=0;featureIndex<numFeaturesA;featureIndex++)
				featureMapA2BInitialized[featureIndex] = Falloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);	
			printf("    first warp the feature images by tentative deformations...\n");			
			ApplyDeformationFieldOnFeatures(featureMapA, featureMapB, numFeaturesA, imageSizeThisLevel, dfFieldB2AThisLevel, featureMapA2BInitialized);   // generation of confidence map requires the update of the feature images after the initial deformation
			printf("    then calculate mutual-saliency map...\n");
			calculateConfidenceMap(featureMapA2BInitialized, featureMapB, maskForControlPointsThisLevel, imageSizeThisLevel, numFeaturesA, levelIndex, distBetweenControlPointsXThisLevel, distBetweenControlPointsYThisLevel, SimilarityMeasure, confidenceMapThisLevel);
		 
			// clear memory
			#pragma omp parallel for shared(featureMapA2BInitialized) private(featureIndex) num_threads(100)
			for (featureIndex=0;featureIndex<numFeaturesA;featureIndex++)
				Ffree3d(featureMapA2BInitialized[featureIndex], imageSizeThisLevel.z, imageSizeThisLevel.x);
			free(featureMapA2BInitialized);

		    // store confidence map
			confidenceMapPreviousLevel = Falloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);
			CopyConfidenceMap(confidenceMapThisLevel, confidenceMapPreviousLevel, imageSizeThisLevel);
		  } // if use mutual-saliency map
	   }

	 // release memory at each level
	 delete imageAthislevel;
	 delete imageBthislevel;
	 delete imageA2Bthisleveluc;
	 Ffree3d(imgA2BThisLevelFloat, imageSizeThisLevel.z, imageSizeThisLevel.x);
	 delete maskthislevel;
	 delete maskforcontrolpointsthislevel;
	 delete confidencemapthislevel;
	 delete deformationthislevel;
	 if (inputInitDeformationOrNot) delete initdeformationthislevel;
	 for (featureIndex=0;featureIndex<numFeaturesA;featureIndex++)
       {
	   UCfree3d(featureMapA[featureIndex], imageSizeThisLevel.z, imageSizeThisLevel.x);
	   UCfree3d(featureMapB[featureIndex], imageSizeThisLevel.z, imageSizeThisLevel.x);
	   }
	 free(featureMapA);
	 free(featureMapB);
	 
	 } //for (i=0; i<numImgLevels; i++)

   times(&endTime);
   durTime=((double)endTime.tms_utime-(double)startTime.tms_utime);
   printf("\n\nDuration=%.2f minutes\n",durTime/(100.0*60.0));

   // release memory at the end
   delete inputimageA;
   delete inputimageB;
   delete mask;
   if (inputInitDeformationOrNot) delete initdeformation;   
}

// ---------------------------------------------------------------------------
void SaveConfidenceMapImage(char* filename,float*** ConfidenceMapImage, unsigned char***mask, int x_size,int y_size, int z_size)
{
  int i,j,k;
  FILE * fp;
  float ***ConfidenceMapToBeSaved;
  ConfidenceMapToBeSaved = Falloc3d(x_size, y_size, z_size);

  if (NULL==(fp=fopen(filename,"w")))
      {
	printf("Unable to open %s!\n",filename);
	exit(0);
      }
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
	  {
	  for (j=0;j<y_size;j++)
	    {
		if (mask[k][i][j]==0)
		  ConfidenceMapToBeSaved[k][i][j]=0.0f;
		else
		  ConfidenceMapToBeSaved[k][i][j]=ConfidenceMapImage[k][i][j];
		}
      fwrite(ConfidenceMapToBeSaved[k][i],sizeof(float),y_size,fp);  
	  }
  fclose(fp);   
  Ffree3d(ConfidenceMapToBeSaved, z_size, x_size);
}

// ---------------------------------------------------------------------------
void DownSampleUCImage3D(unsigned char ***inputImg, Ivector3d imageSize, unsigned char ***imgThisLevel, int resolutionLevel)
{
  Ivector3d imageSizeThisLeve;
  int o,p,q;
  
  imageSizeThisLevel.x = (int)(imageSize.x/(int)pow(2,resolutionLevel));
  imageSizeThisLevel.y = (int)(imageSize.y/(int)pow(2,resolutionLevel));
  imageSizeThisLevel.z = (int)(imageSize.z/(int)pow(2,resolutionLevel));
  
  #pragma omp parallel for private(o,p,q) num_threads(100)
  for (o=0;o<imageSizeThisLevel.z;o++)
   for (p=0;p<imageSizeThisLevel.x;p++)
    for (q=0;q<imageSizeThisLevel.y;q++)
	  imgThisLevel[o][p][q]=inputImg[o*(int)pow(2,resolutionLevel)][p*(int)pow(2,resolutionLevel)][q*(int)pow(2,resolutionLevel)];
}

// ---------------------------------------------------------------------------
void convertFloatImage2UCImage(float ***imageFloat, unsigned char ***imageUC, Ivector3d imageSize)
{
  int x,y,z;
  
  for (z=0;z<imageSize.z;z++)
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  imageUC[z][x][y]=(int)(imageFloat[z][x][y]+0.5);
}

// ---------------------------------------------------------------------------
void convertUCImage2FloatImage(unsigned char ***imageUC, float ***imageFloat, Ivector3d imageSize)
{
  int x,y,z;
  
  for (z=0;z<imageSize.z;z++)
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  imageFloat[z][x][y]=(float)imageUC[z][x][y];
}

// ---------------------------------------------------------------------------
void FFD3D(unsigned char ***imageA, unsigned char ***imageB, unsigned char ****featureMapA, unsigned char ****featureMapB, unsigned char ***maskForControlPoints, unsigned char ***maskForImages, float ***confidenceMap, Ivector3d imageSize, int numFeatures, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defB2A, float ***imageA2B, int levelIndex, int numLevels, int numForegroundVoxels)
{
  int maxNumIterInResolution=5;
  int i;
  int IterIndex, controlPointIndex1, controlPointIndex2, controlPointIndex3, featureIndex;
  Ivector3d realCoorThisControlPoint;
  char filename[1024];
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  
  int incrementControlPoint=1;
  if (levelIndex==1)
    incrementControlPoint=2;   // in the finest resolution, update control point at every 2 control points
  int startControlPoint;

  Fvector3d ***controlPointsBackup;
  controlPointsBackup = Fvector3dalloc3d(numControlPointsX, numControlPointsY, numControlPointsZ);
  
  // get initial enerygy
  float ****featureMapA2BTemp, ****featureMapA2BIterative;
  featureMapA2BTemp = (float****)malloc(sizeof(float***)*numFeatures);	   
  featureMapA2BIterative = (float****)malloc(sizeof(float***)*numFeatures);
  for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
    {
     featureMapA2BTemp[featureIndex]=Falloc3d(imageSize.x, imageSize.y, imageSize.z); 
	 featureMapA2BIterative[featureIndex]=Falloc3d(imageSize.x, imageSize.y, imageSize.z); 
	}
  
  
  float energyAB[maxNumIterInResolution+1];
  for (i=0;i<=maxNumIterInResolution;i++)
    energyAB[i]=0;
  
  
  energyAB[0] = GenerateFeaturesFromDeformationFieldAndCalculateEnergy(featureMapA, featureMapB, featureMapA2BTemp, defB2A, confidenceMap, imageSize, numFeatures);
  CopyFeatureMapsFloat(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures);

  float diffRatio, energyABIterative, energyABTemp;
  energyABIterative = energyAB[0];
  energyABTemp = energyAB[0];
  printf("energy = %f\n", energyAB[0]);
  
  
  for (IterIndex=0;IterIndex<maxNumIterInResolution;IterIndex++)
    {
	  CopyDisplacementAtControlPoints3D(controlPoints, controlPointsBackup, numControlPointsX, numControlPointsY, numControlPointsZ);  // copy: controlPoints -> controlPointsBackup
	  
	  if (IterIndex%2==0)
	    startControlPoint = 1;
	  else
	    startControlPoint = 0;  //update those control points in the even and odd positions in turn
		
		
	  //----------------------------
	  // Calculate displacement at all control points
	  //----------------------------
	  for (controlPointIndex3=0;controlPointIndex3<numControlPointsZ;controlPointIndex3++)
	   for (controlPointIndex1=startControlPoint;controlPointIndex1<numControlPointsX;controlPointIndex1+=incrementControlPoint)
	    for (controlPointIndex2=startControlPoint;controlPointIndex2<numControlPointsY;controlPointIndex2+=incrementControlPoint)
		  {
		  // check if this control point falls into the mask
		  realCoorThisControlPoint.x = controlPointIndex1*distBetweenControlPointsX;
		  realCoorThisControlPoint.y = controlPointIndex2*distBetweenControlPointsY;
		  realCoorThisControlPoint.z = controlPointIndex3*distBetweenControlPointsZ+1;
		  
		  if (maskForControlPoints[realCoorThisControlPoint.z][realCoorThisControlPoint.x][realCoorThisControlPoint.y]>0)
		    {
			printf("\033[33m\n\nUpdate control point (%d, %d, %d) => the real coordinate (%d, %d, %d)\n\033[m", controlPointIndex1, controlPointIndex2, controlPointIndex3, realCoorThisControlPoint.x, realCoorThisControlPoint.y, realCoorThisControlPoint.z);
		    energyABIterative = searchDisplacementAtThisControlPointBasedOnFeatures(imageA, imageB, featureMapA, featureMapB, featureMapA2BTemp, featureMapA2BIterative, confidenceMap, imageSize, numFeaturesA, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defB2A, controlPointIndex1, controlPointIndex2, controlPointIndex3, IterIndex, levelIndex, energyABTemp, numForegroundVoxels); // "controlPoints" will be updated at this control point
			energyABTemp = energyABIterative;
		    } // if
		  else
		    {
		    controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x=0;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y=0;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z=0;
			} // else
		  } // search displacement at each control point
	  
	  //----------------------------------------------------------
	  // Check: after one iteration of search displacement at all control points, check the total energy. If it is relatively stable, no need for further iterations in this resolution.
	  //----------------------------------------------------------
	  printf("\nSmooth displacements at all control points...image level = %d...", levelIndex);
	  UpdateControlPointsWithSmoothIncrement(controlPointsBackup, controlPoints, controlPoints, numControlPointsX, numControlPointsY, numControlPointsZ, IterIndex, maxNumIterInResolution, levelIndex, numLevels, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, maskForControlPoints);
	  smoothDisplacementAtControlPoints3D(controlPoints, controlPoints, numControlPointsX, numControlPointsY, numControlPointsZ, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, levelIndex, numLevels, 0.04*pow(1.5,(float)(3-levelIndex)));  // this factor decreases the smooth kernel to its 2.5% 
	  printf("done!\n");
		
	  printf("Re-generate the *dense* deformation field using now the smoothed displacement at control points...");
	  energyAB[IterIndex+1] = GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(featureMapA, featureMapB, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, confidenceMap, numFeatures, featureMapA2BIterative, defB2A, levelIndex, numLevels);
	  printf("done!\n");

	  diffRatio = (energyAB[IterIndex+1]-energyAB[IterIndex])/energyAB[IterIndex];
	  energyABTemp = energyAB[IterIndex+1];
	  
	  printf("\033[35m\n** %d-th iteration in this level\n\033[m", IterIndex+1);;
	  printf("\033[35m** energy: \033[m");
	  for (i=0;i<=maxNumIterInResolution;i++)
	    printf("\033[35m %f  \033[m", energyAB[i]);
	  printf("\033[35m\n**\n\n\n\033[m");

	  printf("\ndiffRatio = %f, thre = %f\n", diffRatio, 0.10*exp(-pow((float)(numLevels-levelIndex),2.0)/pow((float)numLevels,2.0)));
	  // when total energy doesn't decrease too much, it becomes stable, then terminate iterations at this resolution/level, this requirement gets stricker as resolutions gets finer
	  if ( (fabs(diffRatio)<0.10*exp(-pow((float)(numLevels-levelIndex),2.0)/pow((float)numLevels,2.0)))  || (IterIndex==(maxNumIterInResolution-1)) )
	    {
		GenerateFloatImageFromDeformationField(imageA, imageB, imageA2B, defB2A, imageSize);
	    break;
		}
	  //when total energy increases too much, go back to the previous and terminate iteration at this resolution/level; This requirement gets stricker as resolution gets finer
	  else if ( diffRatio>0.01*exp(-pow((numLevels-levelIndex),2.0)/(float)numLevels) )
		{
		LinearCombinationOfTwoControlPointsMats(controlPoints, controlPointsBackup, controlPoints, 0.85, numControlPointsX, numControlPointsY, numControlPointsZ);
		GenerateFloatImageAndDeformationFieldByFFD(imageA, imageB, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, imageA2B, defB2A, levelIndex, numLevels);
		break;
		}
		
	} // for (IterIndex=0;IterIndex<maxNumIterInResolution;IterIndex++)
	
  Fvector3dfree3d(controlPointsBackup, numControlPointsZ, numControlPointsX);
  for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
    {
	Ffree3d(featureMapA2BTemp[featureIndex], imageSize.z, imageSize.x);
	Ffree3d(featureMapA2BIterative[featureIndex], imageSize.z, imageSize.x);
	}
  free(featureMapA2BTemp);
  free(featureMapA2BIterative);
}

// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
float calculateEnergyOnFeatures(unsigned char ****tempFeatureMap, float ****movingFeatureMap, float ***confidenceMap, Ivector3d imageSize, int numFeatures)
{
  int x,y,z;
  int featureIndex;
  float totalEnergy=0.0;
  float energySquaredAtThisPoint;
  
  for (z=0;z<imageSize.z;z++)
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  {
	  energySquaredAtThisPoint = 0;
	  for (featureIndex=0; featureIndex<numFeatures; featureIndex++) {
	    energySquaredAtThisPoint += pow( ((float)tempFeatureMap[featureIndex][z][x][y] - movingFeatureMap[featureIndex][z][x][y]), 2.0 );
		}
		
	  totalEnergy += energySquaredAtThisPoint * confidenceMap[z][x][y];  // confidence map involved
	  }
	
  totalEnergy /= numFeatures;
  
  return totalEnergy;
}

// ---------------------------------------------------------------------------
void gauss3D(int numPoints1, float std1, int numPoints2, float std2, int numPoints3, float std3, float ***h)
{
  float u[3];
  int i,j,k;
  float sumH=0.0;
  
  for (k=0;k<numPoints3;k++)
   for (i=0;i<numPoints1;i++)
    for (j=0;j<numPoints2;j++)
	  {
	    u[0]=i-(numPoints1-1)/2;
		u[1]=j-(numPoints2-1)/2;
		u[2]=k-(numPoints3-1)/2;
		
		h[k][i][j]=gauss1D(u[0],std1)*gauss1D(u[1],std2)*gauss1D(u[2],std3);
		sumH += h[k][i][j];
	  }
       
  for (k=0;k<numPoints3;k++)
   for (i=0;i<numPoints1;i++)
    for (j=0;j<numPoints2;j++)
	  {
		h[k][i][j]/= sumH;
	  }
}

// ---------------------------------------------------------------------------
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
		sumH += h[i][j];
	  }
        
  for (i=0;i<numPoints2;i++)
    for (j=0;j<numPoints1;j++)
	  {
		h[i][j]/= sumH;
	  }
}

// ---------------------------------------------------------------------------
float gauss1D(float x, float std)
{
  float value;
  value = exp(-pow(x,2.0)/(2.0*pow(std,2.0))) / (std*sqrt(2*G_PI));
  
  return value;
}

// ---------------------------------------------------------------------------
float searchDisplacementAtThisControlPointBasedOnFeatures(unsigned char ***imageA, unsigned char ***imageB, unsigned char ****featureMapA, unsigned char ****featureMapB, float ****featureMapA2BTemp, float ****featureMapA2BIterative, float ***confidenceMap, Ivector3d imageSize, int numFeatures, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defField, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, int NthIter, int NthLevel, float initEnergyAB, int numForegroundVoxels)
{
  int maxIter1=1;
  int maxIter2=1;
  int iterIndex1, iterIndex2;
  int featureIndex;
  int i,j,k;
  int xcoor, ycoor, zcoor;
  float rate;
    
  float maxDisplacementAllowedX;
  if (distBetweenControlPointsX<14)
    maxDisplacementAllowedX = 1.05*distBetweenControlPointsX * pow(1.2, 4-NthLevel);
  else 
    maxDisplacementAllowedX = 0.9*distBetweenControlPointsX * pow(1.15, 3-NthLevel);
	
  float maxDisplacementAllowedY;
  if (distBetweenControlPointsY<14)
    maxDisplacementAllowedY = 1.05*distBetweenControlPointsY * pow(1.2, 4-NthLevel);
  else 
    maxDisplacementAllowedY = 0.9*distBetweenControlPointsY * pow(1.15, 3-NthLevel);
	
  float maxDisplacementAllowedZ;
  if (distBetweenControlPointsZ<4)
    maxDisplacementAllowedZ = 0.5*distBetweenControlPointsZ * pow(1.2, 4-NthLevel);
  else 
    maxDisplacementAllowedZ = 0.65*distBetweenControlPointsZ * pow(1.15, 3-NthLevel);
	
  
  float deltaX = 1.5*0.3 / exp(-pow(NthIter,2.0)/25.0) * (float)distBetweenControlPointsX/3;
  float deltaY = 1.5*0.3 / exp(-pow(NthIter,2.0)/25.0) * (float)distBetweenControlPointsY/3;
  float deltaZ = 0.05 / exp(-pow(NthIter,2.0)/25.0) * (float)distBetweenControlPointsZ/1.5;
  
  float minIncrementX = 1.2*0.5*pow(1.25,(3-NthLevel))*exp((float)distBetweenControlPointsX/50.0);
  float minIncrementY = 1.2*0.5*pow(1.25,(3-NthLevel))*exp((float)distBetweenControlPointsY/50.0);
  float minIncrementZ = 1.0*0.05*pow(1.1,(3-NthLevel))*exp((float)distBetweenControlPointsZ/10.0);
  float maxIncrementX = maxDisplacementAllowedX*0.6;
  float maxIncrementY = maxDisplacementAllowedY*0.6;
  float maxIncrementZ = maxDisplacementAllowedZ*0.6;
  
  printf("\n(deltaX, deltaY, deltaZ) = (%f, %f, %f)\n",deltaX, deltaY, deltaZ);
  printf("minimum increment required in each step = (%f, %f, %f)\n", minIncrementX, minIncrementY, minIncrementZ);
  printf("maximum increment allowed in each step = (%f, %f, %f)\n", maxIncrementX, maxIncrementY, maxIncrementZ);
   
  Ivector3d realCoorThisControlPoint;
  realCoorThisControlPoint.x = controlPointIndex1*distBetweenControlPointsX;
  realCoorThisControlPoint.y = controlPointIndex2*distBetweenControlPointsY;
  realCoorThisControlPoint.z = controlPointIndex3*distBetweenControlPointsZ+1;
  
  float dxIterative=controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x;
  float dyIterative=controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y;
  float dzIterative=controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z;
  float diffEnergyX, diffEnergyY, diffEnergyZ, deltaEnergyX, deltaEnergyY, deltaEnergyZ;
  float deltaEnergyXSquared, deltaEnergyYSquared, deltaEnergyZSquared;
  float gradientX, gradientY, gradientZ;
  float incrementX, incrementY, incrementZ;
  float actualIncrementX, actualIncrementY, actualIncrementZ;
  float previousIncrementX=0, previousIncrementY=0, previousIncrementZ=0;
  float eta;
  
  int affectingRadiusXAhead = distBetweenControlPointsX*2;
  int affectingRadiusXBehind = distBetweenControlPointsX*2;
  int affectingRadiusYAhead = distBetweenControlPointsY*2;
  int affectingRadiusYBehind = distBetweenControlPointsY*2;
  int affectingRadiusZAhead = distBetweenControlPointsZ*2;
  int affectingRadiusZBehind = distBetweenControlPointsZ*2;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  
  Fvector3d ***defFieldTemp;
  defFieldTemp = Fvector3dalloc3d(imageSize.x, imageSize.y, imageSize.z);
  for (k=-affectingRadiusZAhead;k<affectingRadiusZBehind;k++)
   for (i=-affectingRadiusXAhead;i<affectingRadiusXBehind;i++)
    for (j=-affectingRadiusYAhead;j<affectingRadiusYBehind;j++)
	  {
	  xcoor = realCoorThisControlPoint.x + i;
	  ycoor = realCoorThisControlPoint.y + j;
	  zcoor = realCoorThisControlPoint.z + k;
	  if ( (xcoor>=0)&(xcoor<imageSize.x)&(ycoor>=0)&(ycoor<imageSize.y)&(zcoor>=0)&(zcoor<imageSize.z) )
	    {
		defFieldTemp[zcoor][xcoor][ycoor].x = defField[zcoor][xcoor][ycoor].x;
		defFieldTemp[zcoor][xcoor][ycoor].y = defField[zcoor][xcoor][ycoor].y;
		defFieldTemp[zcoor][xcoor][ycoor].z = defField[zcoor][xcoor][ycoor].z;
		}
	  }
  
  float energyAB = initEnergyAB;
  float previousEnergyAB = energyAB;
  float ENERGY[maxIter1+1];
  for (i=0;i<=maxIter1;i++)
    ENERGY[i] = 0;
  
  
  for (iterIndex1=0;iterIndex1<maxIter1;iterIndex1++)
    {
	ENERGY[iterIndex1]=energyAB;
	previousEnergyAB = energyAB;
	
	//------------------------------------
	// search in x direction
	//-----------------------------------
	for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
	  {
	    printf("X: Iter #%d, %d\n", iterIndex1, iterIndex2);
		controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x += deltaX;

		diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, deltaX, 'X', NthLevel);
		// gradient
		gradientX = diffEnergyX/deltaX;
		printf("\tdiffEnergyX/energyAB = %f, gradientX = %f\n", diffEnergyX/energyAB, gradientX);  // for debug
		
		// if change is too small, no need to update displacement at this control point, break;
		if ( fabs(diffEnergyX/energyAB) <= 0.0000025/ ((float)distBetweenControlPointsX/3.0)* exp(-pow(NthIter,2.0)/5.0) )
          {
		   controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x = dxIterative;
		   printf("\tGradientX is too small!\n");
		   printf("\tdx = %f, dy = %f, dz=%f, diffEnergyX = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyX);
           break;
		  } // if gradient is too small, break;
        else // else, proceed along the gradient direction
		  {
		   // get a proper eta
		   eta = 200*exp(-pow(iterIndex1,2.0)/100.0)*exp(-pow(iterIndex2,2.0)/100.0)* pow(10, -8)*(energyAB/(float)(numForegroundVoxels*distBetweenControlPointsX*distBetweenControlPointsY*distBetweenControlPointsZ)) * pow(4.0, (float)(3-NthLevel));
            
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
		   printf("\teta = %f, gradient = %f\n",eta, gradientX);
           printf("\ttentative incrementX = %f", incrementX);
		   
		   
		   // contraint 1 on increment : it can not exceed maximum displacement allowed
           dxIterative = dxIterative + incrementX;
           if (fabs(dxIterative)>maxDisplacementAllowedX)
             dxIterative = sign(dxIterative)*maxDisplacementAllowedX;
			 
		   // constraint 2 on increment: it must guarantee that dx proceeds in the right direction; otherwise, increment=0.1*increment
		   actualIncrementX = dxIterative - controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x + deltaX;
		   printf(", actual incrementX = %f\n", actualIncrementX);
           controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x = dxIterative;
        
		
		   if ( ((incrementX==-previousIncrementX)&(fabs(incrementX)==maxIncrementX)) || ((incrementX==-previousIncrementX)&(fabs(incrementX)==minIncrementX)) )
		    {// abondon this move
			printf("\tAbondon this move\n");
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x -= incrementX;
			dxIterative -= incrementX;
			printf("\tdx = %f, dy = %f, dz = %f\n\tenergy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			break;
			}

		   CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind,affectingRadiusZAhead, affectingRadiusZBehind);
		   diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementX, 'X', NthLevel, NNO);
		   energyAB = previousEnergyAB + diffEnergyX;
		   rate = diffEnergyX/previousEnergyAB;
		   printf("\tenergyAB = %f, \n\tpreviousEnergyAB = %f, \n\tchangeRate = %f\n", energyAB, previousEnergyAB, rate );   // for debug
		   
		   
		   if ( (fabs(rate)<pow(10,-5.0)*pow(1.5,(NthLevel-3))) || ((fabs(dxIterative)==maxDisplacementAllowedX)&&(rate<0.005)) )
		    {
			// too little change to be continued or the displacement has reached the maximum allowed
			printf("\tToo little change or the displacement has reached maximum allowed -- no need to move in this direction at this stage!\n");
			printf("\tdx = %f, dy = %f, dz = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z);
			previousEnergyAB = energyAB;
			previousIncrementX = 0;
			CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
		    break;
			}
		   else if ( rate > 0.005 * pow(0.8, (3-NthLevel)) )
		    {
			// moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
			dxIterative = dxIterative - incrementX + 0.2*incrementX*0.005* pow(0.8, (3-NthLevel))/rate;
			actualIncrementX = 0.2*incrementX*0.005* pow(0.8, (3-NthLevel))/rate;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x = dxIterative;
			CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementX, 'X', NthLevel, NNO);
			energyAB = previousEnergyAB + diffEnergyX;
			printf("\thas moved too much.. Aho, move back!\n");
			printf("\tdx = %f, dy = %f, dz=%f, diffEnergyX = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyX, energyAB);
			if (energyAB>(1+0.005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  dxIterative = dxIterative - 0.2*incrementX*0.005* pow(0.8, (3-NthLevel))/rate;
			  actualIncrementX = 0;
			  controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x = dxIterative;
			  printf("\thas moved too much.. Aho, move back!\n");
			  printf("\tdx = %f, dy = %f, dz=%f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			  diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementX, 'X', NthLevel, YYES);
			  break;
			  }
			else if (energyAB>0.9995*previousEnergyAB && energyAB<=(1+0.0005*pow(0.8,(3-NthLevel)))*previousEnergyAB) 
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;			
			  }
			else
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  }
			}
		   else if ( (rate > 0.000025 * pow(0.8, (3-NthLevel))) && (rate <= 0.005 * pow(0.8, (3-NthLevel))) )
		    {
			// moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
			dxIterative = dxIterative - 0.5*incrementX;
			actualIncrementX = 0.5*incrementX;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x = dxIterative;
			CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementX, 'X', NthLevel, NNO);
			energyAB = previousEnergyAB + diffEnergyX;
			printf("\thas moved too much.. Aho, move back!\n");
			printf("\tdx = %f, dy = %f, dz=%f, diffEnergyX = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyX, energyAB);
			if (energyAB>(1+0.005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  dxIterative = dxIterative - 0.5*incrementX;
			  actualIncrementX = 0;
			  controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x = dxIterative;
			  printf("\thas moved too much.. Aho, move back!\n");
			  printf("\tdx = %f, dy = %f, dz=%f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			  diffEnergyX = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementX, 'X', NthLevel, YYES);
			  break;
			  }
			else if (energyAB>0.9995*previousEnergyAB && energyAB<=(1+0.0005*pow(0.8,(3-NthLevel)))*previousEnergyAB) 
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;			
			  }
			else
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;
			  }
			}
		   else
		    {
			previousEnergyAB = energyAB; 
			previousIncrementX = incrementX;
			CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
            } // else, increment is proper
		   printf("\tdx = %f, dy = %f, dz=%f, diffEnergyX = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyX);
          } // else, proceed along the gradient direction
      } //for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)

	//-----------------------------------
	// search in y direction
	//-----------------------------------
	for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
	  {
	    printf("Y: Iter #%d, %d\n", iterIndex1, iterIndex2);
		controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y += deltaY;

		diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, deltaY, 'Y', NthLevel);
		// gradient
		gradientY = diffEnergyY/deltaY;
		printf("\tdiffEnergyY/energyAB = %f, gradientY = %f\n", diffEnergyY/energyAB, gradientY);  // for debug
		
		if ( fabs(diffEnergyY/energyAB) <= 0.0000025/ ((float)distBetweenControlPointsY/3)* exp(-pow(NthIter,2.0)/5.0) )
          {
		   controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y = dyIterative;
		   printf("\tGradientY is too small!\n");
		   printf("\tdx = %f, dy = %f, dz = %f, diffEnergyY = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyY);
           break;
		  } // if gradient is too small, break;
        else // else, proceed along the gradient direction
		  {
		   // get a proper eta
		   eta = 200*exp(-pow(iterIndex1,2.0)/100.0)*exp(-pow(iterIndex2,2.0)/100.0)* pow(10, -8)*(energyAB/(float)(numForegroundVoxels*distBetweenControlPointsX*distBetweenControlPointsY*distBetweenControlPointsZ)) * pow(4.0, (float)(3-NthLevel));

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
             incrementY = sign(incrementY)*maxIncrementY;  // do not move too much in one iteration
		   printf("\teta = %f, gradient = %f\n",eta, gradientY);
		   printf("\ttentative incrementY = %f", incrementY);
           
		   // contraint 1 on increment : it can not exceed 90% of the distance between two adjacent control points
           dyIterative = dyIterative + incrementY;
           if (fabs(dyIterative)>maxDisplacementAllowedY)
             dyIterative = sign(dyIterative)*maxDisplacementAllowedY;
             
           
		   // constraint 2 on increment: it must guarantee that dx proceeds in the right direction; otherwise, increment=0.1*increment
           actualIncrementY = dyIterative - controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y + deltaY;
		   printf(", actual incrementY = %f\n", actualIncrementY);
		   controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y = dyIterative;
		   
		   
		   if ( ((incrementY==-previousIncrementY)&(fabs(incrementY)==maxIncrementY)) || ((incrementY==-previousIncrementY)&(fabs(incrementY)==minIncrementY)) )
		    {// abondon this move
			printf("\tAbondon this move\n");
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y -= incrementY;
			dyIterative -= incrementY;
			printf("\tdx = %f, dy = %f, dz = %f\n\tenergy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			break;
			} 

		   CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
		   diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementY, 'Y', NthLevel, NNO);
		   
		   energyAB = previousEnergyAB + diffEnergyY;
		   rate = diffEnergyY/previousEnergyAB;
		   printf("\tenergyAB = %f, \n\tpreviousEnergyAB = %f, \n\tchangeRate = %f\n", energyAB, previousEnergyAB, rate );   // for debug
		   
		   
		   if ( (fabs(rate)<pow(10,-5.0)*pow(1.5,(NthLevel-3))) || ((fabs(dyIterative)==maxDisplacementAllowedY)&&(rate<0.005)) )
		    {
			// too little change to be continued or the displacement has reached the maximum allowed
			printf("\tToo little change or the displacement has reached maximum allowed -- no need to move in this direction at this stage!\n");
			printf("\tdx = %f, dy = %f, dz=%f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z);
			previousEnergyAB = energyAB;
			previousIncrementY = 0;
			CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
		    break;
			}
		   else if ( rate > 0.005 * pow(0.8, (3-NthLevel)) )
		    {
			// moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
			dyIterative = dyIterative - incrementY + 0.2*incrementY*0.005* pow(0.8, (3-NthLevel))/rate;
			actualIncrementY = 0.2*incrementY*0.005* pow(0.8, (3-NthLevel))/rate;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y = dyIterative;
			CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementY, 'Y', NthLevel, NNO);
			energyAB = previousEnergyAB + diffEnergyY;
			printf("\thas moved too much.. Aho, move back!\n");
			printf("\tdx = %f, dy = %f, dz=%f, diffEnergyY = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyY, energyAB);
			if (energyAB>(1+0.005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  dyIterative = dyIterative - 0.2*incrementY*0.005* pow(0.8, (3-NthLevel))/rate;
			  actualIncrementY = 0;
			  controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y = dyIterative;
			  printf("\thas moved too much.. Aho, move back!\n");
			  printf("\tdx = %f, dy = %f, dz=%f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			  diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementY, 'Y', NthLevel, YYES);
			  break;
			  }
			else if (energyAB>0.9995*previousEnergyAB && energyAB<=(1+0.0005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;			
			  }
			else
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  }
			}
		   else if ( (rate > 0.000025 * pow(0.8, (3-NthLevel))) && (rate > 0.005 * pow(0.8, (3-NthLevel))) )
		    {
			// moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
			dyIterative = dyIterative - 0.5*incrementY;
			actualIncrementY = 0.5*incrementY;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y = dyIterative;
			CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementY, 'Y', NthLevel, NNO);
			energyAB = previousEnergyAB + diffEnergyY;
			printf("\thas moved too much.. Aho, move back!\n");
			printf("\tdx = %f, dy = %f, dz=%f, diffEnergyY = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyY, energyAB);
			if (energyAB>(1+0.005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  dyIterative = dyIterative - 0.5*incrementY;
			  actualIncrementY = 0;
			  controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y = dyIterative;
			  printf("\thas moved too much.. Aho, move back!\n");
			  printf("\tdx = %f, dy = %f, dz=%f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			  diffEnergyY = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementY, 'Y', NthLevel, YYES);
			  break;
			  }
			else if (energyAB>0.9995*previousEnergyAB && energyAB<=(1+0.0005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;			
			  }
			else
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;
			  }		
			}
		   else
		    {
			previousEnergyAB = energyAB; 
			previousIncrementY = incrementY;
			CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
            } // else, increment is proper
		   printf("\tdx = %f, dy = %f, dz=%f, diffEnergyY = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyY);
		  } // else, proceed along the gradient direction
      } //for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)

	//-----------------------------------
	// search in z direction
	//-----------------------------------
	controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z=0;
	
	for (iterIndex2=0; iterIndex2<maxIter2; iterIndex2++)
	  {
	    printf("Z: Iter #%d, %d\n", iterIndex1, iterIndex2);
		controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z += deltaZ;
		
		diffEnergyZ = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, deltaZ, 'Z', NthLevel);
		// gradient
		gradientZ = diffEnergyZ/deltaZ;
		printf("\tdiffEnergyZ/energyAB = %f, gradientZ = %f\n", diffEnergyZ/energyAB, gradientZ);  // for debug
		
		if ( fabs(diffEnergyZ/energyAB) <= 0.0000025/ ((float)distBetweenControlPointsZ/3.0)* exp(-pow(NthIter,2.0)/5.0) )
          {
		   controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z = dzIterative;
		   printf("\tGradientZ is too small!\n");
		   printf("\tdx = %f, dy = %f, dz = %f, diffEnergyZ = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyZ);
           break;
		  } // if gradient is too small, break;
        else // else, proceed along the gradient direction
		  {
		   // get a proper eta
		   eta = 70*exp(-pow(iterIndex1,2.0)/100.0)*exp(-pow(iterIndex2,2.0)/100.0)* pow(10, -8)*(energyAB/(float)(numForegroundVoxels*distBetweenControlPointsX*distBetweenControlPointsY*distBetweenControlPointsZ)) * pow(3.5, (float)(3-NthLevel));

           if (diffEnergyZ/energyAB > 0.5)   
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*0.1*diffEnergyZ/energyAB;   // in case of big gradient, do not behave too aggressively
           else if ( (diffEnergyZ/energyAB>0.2) & (diffEnergyZ/energyAB<=0.5) )
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*0.5*diffEnergyZ/energyAB;
           else if ( (diffEnergyZ/energyAB > 0.1) & (diffEnergyZ/energyAB<=0.2) )
             eta = exp(-pow(NthIter,2.0)/2.0)*exp(-pow(iterIndex1,2.0)/10.0) *eta*diffEnergyZ/energyAB;
           
		   // get a proper increment 
           incrementZ = - eta*gradientZ;
           		   
           if (fabs(incrementZ)<minIncrementZ)
			 incrementZ = sign(incrementZ)*minIncrementZ;
           else if (fabs(incrementZ)>maxIncrementZ)
             incrementZ = sign(incrementZ)*maxIncrementZ;  // do not move too much in one iteration
		   printf("\teta = %f, gradient = %f\n",eta, gradientZ);
		   printf("\ttentative incrementZ = %f", incrementZ);
           
		   // contraint 1 on increment : it can not exceed 90% of the distance between two adjacent control points
           dzIterative = dzIterative + incrementZ;
           if (fabs(dzIterative)>maxDisplacementAllowedZ)
             dzIterative = sign(dzIterative)*maxDisplacementAllowedZ;
             
           
		   // constraint 2 on increment: it must guarantee that dx proceeds in the right direction; otherwise, increment=0.1*increment
           actualIncrementZ = dzIterative - controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z + deltaZ;
		   printf(", actual incrementZ = %f\n", actualIncrementZ);
		   controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z = dzIterative;
		   
		   
		   if ( ((incrementZ==-previousIncrementZ)&(fabs(incrementZ)==maxIncrementZ)) || ((incrementZ==-previousIncrementZ)&(fabs(incrementZ)==minIncrementZ)) )
		    {// abondon this move
			printf("\tAbondon this move\n");
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z -= incrementZ;
			dzIterative -= incrementZ;
			printf("\tdx = %f, dy = %f, dz = %f\n\tenergy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			break;
			} 
			
		   CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
		   diffEnergyZ = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementZ, 'Z', NthLevel, NNO);
		   
		   energyAB = previousEnergyAB + diffEnergyZ;
		   rate = diffEnergyZ/previousEnergyAB;
		   printf("\tenergyAB = %f, \n\tpreviousEnergyAB = %f, \n\tchangeRate = %f\n", energyAB, previousEnergyAB, rate );   // for debug
		   
		   
		   if ( (fabs(rate)<pow(10,-5.0)*pow(1.5,(NthLevel-3))) || ((fabs(dzIterative)==maxDisplacementAllowedZ)&&(rate<0.005)) )
		    {
			// too little change to be continued or the displacement has reached the maximum allowed
			printf("\tToo little change or the displacement has reached maximum allowed -- no need to move in this direction at this stage!\n");
			printf("\tdx = %f, dy = %f, dz=%f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z);
			previousEnergyAB = energyAB;
			previousIncrementZ = 0;
			CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
		    break;
			}
		   else if ( rate > 0.005 * pow(0.8, (3-NthLevel)) )
		    {
			// moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
			dzIterative = dzIterative - incrementZ + 0.2*incrementZ*0.005* pow(0.8, (3-NthLevel))/rate;
			actualIncrementZ = 0.2*incrementZ*0.005* pow(0.8, (3-NthLevel))/rate;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z = dzIterative;
			CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			diffEnergyZ = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementZ, 'Z', NthLevel, NNO);
			energyAB = previousEnergyAB + diffEnergyZ;
			printf("\thas moved too much.. Aho, move back!\n");
			printf("\tdx = %f, dy = %f, dz=%f, diffEnergyZ = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyZ, energyAB);
			if (energyAB>(1+0.005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  dzIterative = dzIterative - 0.2*incrementZ*0.005* pow(0.8, (3-NthLevel))/rate;
			  actualIncrementZ = 0;
			  controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z = dzIterative;
			  printf("\thas moved too much.. Aho, move back!\n");
			  printf("\tdx = %f, dy = %f, dz=%f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			  diffEnergyZ = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementZ, 'Z', NthLevel, YYES);
			  break;
			  }
			else if (energyAB>0.9995*previousEnergyAB && energyAB<=(1+0.0005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;			
			  }
			else
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  }
			}
		   else if ( (rate > 0.000025 * pow(0.8, (3-NthLevel))) && (rate > 0.005 * pow(0.8, (3-NthLevel))) )
		    {
			// moved too much in the gradient direction so that the total energy started to increase. Therefore, need to move back
			dzIterative = dzIterative - 0.5*incrementZ;
			actualIncrementZ = 0.5*incrementZ;
			controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z = dzIterative;
			CopyDispAtTheBlockInducedByThisControlPoint(defField, defFieldTemp, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			diffEnergyZ = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defFieldTemp, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementZ, 'Z', NthLevel, NNO);
			energyAB = previousEnergyAB + diffEnergyZ;
			printf("\thas moved too much.. Aho, move back!\n");
			printf("\tdx = %f, dy = %f, dz=%f, diffEnergyY = %f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyY, energyAB);
			if (energyAB>(1+0.005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  dzIterative = dzIterative - 0.5*incrementZ;
			  actualIncrementZ = 0;
			  controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z = dzIterative;
			  printf("\thas moved too much.. Aho, move back!\n");
			  printf("\tdx = %f, dy = %f, dz=%f\n\tupdated energy = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, previousEnergyAB);
			  diffEnergyZ = UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(featureMapA, featureMapB, confidenceMap, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPoints, defField, numFeatures, featureMapA2BTemp, featureMapA2BIterative, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind, controlPointIndex1, controlPointIndex2, controlPointIndex3, actualIncrementZ, 'Z', NthLevel, YYES);
			  break;
			  }
			else if (energyAB>0.9995*previousEnergyAB && energyAB<=(1+0.0005*pow(0.8,(3-NthLevel)))*previousEnergyAB)
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;			
			  }
			else
			  {
			  previousEnergyAB = energyAB;
			  CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			  break;
			  }		
			}
		   else
		    {
			previousEnergyAB = energyAB; 
			previousIncrementZ = incrementZ;
			CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(featureMapA2BTemp, featureMapA2BIterative, imageSize, numFeatures, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
			CopyDispAtTheBlockInducedByThisControlPoint(defFieldTemp, defField, imageSize, realCoorThisControlPoint, affectingRadiusXAhead, affectingRadiusXBehind, affectingRadiusYAhead, affectingRadiusYBehind, affectingRadiusZAhead, affectingRadiusZBehind);
            } // else, increment is proper
		   printf("\tdx = %f, dy = %f, dz=%f, diffEnergyZ = %f\n\n", controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].x, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].y, controlPoints[controlPointIndex3][controlPointIndex1][controlPointIndex2].z, diffEnergyZ);
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
  Fvector3dfree3d(defFieldTemp, imageSize.z, imageSize.x);
  
  return energyAB;
}

// ---------------------------------------------------------------------------
void GenerateFloatImageAndDeformationFieldByFFD(unsigned char ***imageA, unsigned char ***imageB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***imageA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels)
{
  int x,y,z;
  int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
  int l, m,n;
  float dx, dy, dz;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
  float u,v,w;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  
 
  for (z=0;z<imageSize.z;z++)
   { // for each slice
   printf("\tslice #%d...\n",z);
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  { //for each point
		dx=0;
		dy=0;
		dz=0;
		
		// find out effecting control points
		indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPointsX)-1;
		indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPointsY)-1;
		indexFirstEffectiveControlPointZ = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ)-1;
		u = (float)x/(float)distBetweenControlPointsX - floor((float)x/(float)distBetweenControlPointsX);
		v = (float)y/(float)distBetweenControlPointsY - floor((float)y/(float)distBetweenControlPointsY);
		w = (float)(z-1)/(float)distBetweenControlPointsZ - floor((float)(z-1)/(float)distBetweenControlPointsZ);
		
		// calculate displacement at this voxel (x,y) under the influence of all effective control points
		for (l=0;l<=3;l++)
		 for (m=0;m<=3;m++)
		  for (n=0;n<=3;n++)
		    {
			indexEffectiveControlPointX = l+indexFirstEffectiveControlPointX;
			indexEffectiveControlPointY = m+indexFirstEffectiveControlPointY;
			indexEffectiveControlPointZ = n+indexFirstEffectiveControlPointZ;
			
			if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0 & indexEffectiveControlPointZ<numControlPointsZ) )
			  {
				dx += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
			  }
			}
		
		// save deformation field
		defField[z][x][y].x = dx;
		defField[z][x][y].y = dy;
		defField[z][x][y].z = dz;
		
		//interpolate intensity (float)
		imageA2BFloat[z][x][y] = trilinearInterpolation(imageA, (float)x+dx, (float)y+dy, (float)z+dz, imageSize);	
	  } // for each point on this slice
  } // for each slice
}

// ---------------------------------------------------------------------------
// speed-up version of GenerateFloatImageAndDeformationFieldByFFD
void GenerateFloatImageAndDeformationFieldByFFD2(unsigned char ***imageA, unsigned char ***imageB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***imageA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels, int iterIdx, int AdditionOrComposition, float *weightsX, float *weightsY, float *weightsZ, int method, int inputInitDeformationOrNot, int encourageDiffeomorphismOrNot)
{
  int x,y,z;
  float xx,yy,zz;
  int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
  float w1,w2,w3,w4,w5,w6,w7,w8;
  int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
  int l, m,n;
  float dx, dy, dz;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
  int u,v,w;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  int numPoints = numControlPointsX * numControlPointsY * numControlPointsZ;
  Fvector3d ***defFieldTemp;
  float weight, sumWeight;
  
  // pre-calculate the weights for B-spline, added on June 8, 2009 (begin)
  int maxSizeXY = (int)MAX((float)imageSize.x, (float)imageSize.y);
  int maxSizeXYZ = (int)MAX((float)maxSizeXY, (float)imageSize.z);
 
  if ( ~(levelIndex==numLevels && iterIdx==0) && (AdditionOrComposition==1) )  // 0 for addition and 1 for composition
    defFieldTemp = Fvector3dalloc3d(imageSize.x, imageSize.y, imageSize.z);  
	
	
  #pragma omp parallel for shared(method,distBetweenControlPointsX,distBetweenControlPointsY,distBetweenControlPointsZ,controlPoints,weightsX,weightsY,weightsZ) private(x,y,z,u,v,w,indexEffectiveControlPointX,indexEffectiveControlPointY,indexEffectiveControlPointZ,indexFirstEffectiveControlPointX,indexFirstEffectiveControlPointY,indexFirstEffectiveControlPointZ,dx,dy,dz,l,m,n,weight) num_threads(100)
  for (z=0;z<imageSize.z;z++)
   { // for each slice
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  { //for each point
	  dx=0.0;
	  dy=0.0;
	  dz=0.0;

	  u = (int)floor((float)x/(float)distBetweenControlPointsX);
      v = (int)floor((float)y/(float)distBetweenControlPointsY);
      w = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ);

	  if (method==0)  // nearest neighborhood weighting mechanism
		{ // if method==0
		indexEffectiveControlPointX = u;
		indexEffectiveControlPointY = v;
		indexEffectiveControlPointZ = w;
		
        dx = controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
		dy = controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
		dz = controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
		} // if method==0
	  else if (method==1)  // trilinear weighting mechanism
	    { // if method==1
		sumWeight = 0.0;
		for (l=0;l<=1;l++)
		  for (m=0;m<=1;m++)
		    for (n=0;n<=1;n++)
			  {
			  indexEffectiveControlPointX = u+l;
			  indexEffectiveControlPointY = v+m;
			  indexEffectiveControlPointZ = w+n;
			  
			  if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				{
				weight = weightsX[x+l*imageSize.x] * weightsY[y+m*imageSize.y] * weightsZ[z+n*imageSize.z];
				sumWeight += weight;
				dx += weight * controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += weight * controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += weight * controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
				}
			  }
		
		dx /= sumWeight;
		dy /= sumWeight;
		dz /= sumWeight;
		} // if method==1
	  else if (method==2)
		{  // if method==2
		// find out effecting control points
		indexFirstEffectiveControlPointX = u-1;
		indexFirstEffectiveControlPointY = v-1;
		indexFirstEffectiveControlPointZ = w-1;

		// calculate displacement at this voxel (x,y) under the influence of all effective control points
		sumWeight = 0.0;
		for (l=0;l<=3;l++)
		 for (m=0;m<=3;m++)
		  for (n=0;n<=3;n++)
		    {
			indexEffectiveControlPointX = l+indexFirstEffectiveControlPointX;
			indexEffectiveControlPointY = m+indexFirstEffectiveControlPointY;
			indexEffectiveControlPointZ = n+indexFirstEffectiveControlPointZ;

			if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0 & indexEffectiveControlPointZ<numControlPointsZ) )
			  {
				weight = weightsX[x+l*imageSize.x]*weightsY[y+m*imageSize.y]*weightsZ[z+n*imageSize.z];
				sumWeight += weight;
				dx += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
			  }
			}
		
		dx /= sumWeight;
		dy /= sumWeight;
		dz /= sumWeight;
		} // if method==2
	  
	  
	  // deformation field
	  if ( levelIndex==numLevels && iterIdx==0 && inputInitDeformationOrNot==NNO )  // in the coarsest resolution
		{
		defField[z][x][y].x = dx;
		defField[z][x][y].y = dy;
		defField[z][x][y].z = dz;
		}
      else if ( AdditionOrComposition==0 )  // 0 for addition and 1 for composition
		{ 
		defField[z][x][y].x += dx;
		defField[z][x][y].y += dy;
		defField[z][x][y].z += dz;		  
		}
	  else if ( AdditionOrComposition==1 )  // 0 for addition and 1 for composition
		{
		defFieldTemp[z][x][y].x = defField[z][x][y].x;
		defFieldTemp[z][x][y].y = defField[z][x][y].y;
		defFieldTemp[z][x][y].z = defField[z][x][y].z;
		defField[z][x][y].x = dx;
		defField[z][x][y].y = dy;
		defField[z][x][y].z = dz;
		}
		  
		
		                						
	  if ( AdditionOrComposition==0)
		{
		//interpolate intensity (float)
		imageA2BFloat[z][x][y] = trilinearInterpolation(imageA, (float)x+defField[z][x][y].x, (float)y+defField[z][x][y].y, (float)z+defField[z][x][y].z, imageSize);	
		}
		
	  } // for each point on this slice
   } // for each slice

   // in case we need to composite deformation field
   if ( ~(levelIndex==numLevels && iterIdx==0) && (AdditionOrComposition==1) ) // composition
     { // if composition
	 printf("\tComposing deformations...\n");
	 
	 #pragma omp parallel for private(x,y,z,xx,yy,zz,xFloor,xCeil,yFloor,yCeil,zFloor,zCeil,dx,dy,dz,w1,w2,w3,w4,w5,w6,w7,w8) num_threads(100)
	 for (z=0;z<imageSize.z;z++)
       { // for each slice
       for (x=0;x<imageSize.x;x++)
         for (y=0;y<imageSize.y;y++)
		   { // for each point on the slice
		   xx = MAX(MIN((float)x+defField[z][x][y].x, (float)(imageSize.x-1)), 0);
		   yy = MAX(MIN((float)y+defField[z][x][y].y, (float)(imageSize.y-1)), 0);
		   zz = MAX(MIN((float)z+defField[z][x][y].z, (float)(imageSize.z-1)), 0);
		   
		   xFloor = (int)floor(xx);
		   xCeil  = (int)ceil(xx);
		   yFloor = (int)floor(yy);
		   yCeil  = (int)ceil(yy);
	       zFloor = (int)floor(zz);
		   zCeil  = (int)ceil(zz);
		   
		   if ( (xFloor==xCeil)&(xFloor!=(imageSize.x-1)) )  xCeil += 1;
		   if ( (xFloor==xCeil)&(xFloor==(imageSize.x-1)) )  xFloor -= 1;
		   if ( (yFloor==yCeil)&(yFloor!=(imageSize.y-1)) )  yCeil += 1;
		   if ( (yFloor==yCeil)&(yFloor==(imageSize.y-1)) )  yFloor -= 1;
		   if ( (zFloor==zCeil)&(zFloor!=(imageSize.z-1)) )  zCeil += 1;
		   if ( (zFloor==zCeil)&(zFloor==(imageSize.z-1)) )  zFloor -= 1;
		   
		   w1 = (x-xFloor)*(y-yFloor)*(z-zFloor);
		   w2 = (x-xFloor)*(y-yFloor)*(zCeil-z);
		   w3 = (x-xFloor)*(yCeil-y)*(z-zFloor);
		   w4 = (x-xFloor)*(yCeil-y)*(zCeil-z);
		   w5 = (xCeil-x)*(y-yFloor)*(z-zFloor);
		   w6 = (xCeil-x)*(y-yFloor)*(zCeil-z);
		   w7 = (xCeil-x)*(yCeil-y)*(z-zFloor);
		   w8 = (xCeil-x)*(yCeil-y)*(zCeil-z);
		   
		   // incremental displacement at (x,y,z)
		   dx =  w1*defFieldTemp[zCeil][xCeil][yCeil].x
		       + w2*defFieldTemp[zFloor][xCeil][yCeil].x
		       + w3*defFieldTemp[zCeil][xCeil][yFloor].x
			   + w4*defFieldTemp[zFloor][xCeil][yFloor].x
			   + w5*defFieldTemp[zCeil][xFloor][yCeil].x
		       + w6*defFieldTemp[zFloor][xFloor][yCeil].x
		       + w7*defFieldTemp[zCeil][xFloor][yFloor].x
			   + w8*defFieldTemp[zFloor][xFloor][yFloor].x;
			   
		   dy =  w1*defFieldTemp[zCeil][xCeil][yCeil].y
		       + w2*defFieldTemp[zFloor][xCeil][yCeil].y
		       + w3*defFieldTemp[zCeil][xCeil][yFloor].y
			   + w4*defFieldTemp[zFloor][xCeil][yFloor].y
			   + w5*defFieldTemp[zCeil][xFloor][yCeil].y
		       + w6*defFieldTemp[zFloor][xFloor][yCeil].y
		       + w7*defFieldTemp[zCeil][xFloor][yFloor].y
			   + w8*defFieldTemp[zFloor][xFloor][yFloor].y;
			   
		   dz =  w1*defFieldTemp[zCeil][xCeil][yCeil].z
		       + w2*defFieldTemp[zFloor][xCeil][yCeil].z
		       + w3*defFieldTemp[zCeil][xCeil][yFloor].z
			   + w4*defFieldTemp[zFloor][xCeil][yFloor].z
			   + w5*defFieldTemp[zCeil][xFloor][yCeil].z
		       + w6*defFieldTemp[zFloor][xFloor][yCeil].z
		       + w7*defFieldTemp[zCeil][xFloor][yFloor].z
			   + w8*defFieldTemp[zFloor][xFloor][yFloor].z;

		   defField[z][x][y].x += dx;
		   defField[z][x][y].y += dy;
		   defField[z][x][y].z += dz;

		   //interpolate intensity (float)
		   imageA2BFloat[z][x][y] = trilinearInterpolation(imageA, (float)x+defField[z][x][y].x, (float)y+defField[z][x][y].y, (float)z+defField[z][x][y].z, imageSize);	

		   } // for each point on the slice
	   }// for each slice
	 }  // if composition

   if ( ~(levelIndex==numLevels && iterIdx==0) && (AdditionOrComposition==1) )  // 0 for addition and 1 for composition
     Fvector3dfree3d(defFieldTemp, imageSize.z, imageSize.x);
	
	
   // ***********************
   // below are added on 01/24/2011 to modulate deformation properties.
   // ***********************
   if (encourageDiffeomorphismOrNot==YYES)
	{ //if (encourageDiffeomorphismOrNot==YYES)
	float ***Jacobian;
	Fvector3d ***defTemp;
	Matrix *Jaco ;
	int extremeJacobian=NNO;
	float chuzi;
	int it;
	int xx, yy, zz, xxL, xxR, yyL, yyR, zzL, zzR, tx,ty,tz;

    Jacobian = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
	defTemp = Fvector3dalloc3d(imageSize.x, imageSize.y, imageSize.z);  
	CreateMatrix(&Jaco, 3, 3);
	
	for (it=0;it<2;it++)
	{ // two iterations to reduce extreme jacobians, at a cost of little accuray
	if ( (it==0)||(it>0&&extremeJacobian==YYES) )
	{
	// backup
	for (z=0;z<imageSize.z;z++)
	  for (x=0;x<imageSize.x;x++)
		for (y=0;y<imageSize.y;y++)
		  {
		  defTemp[z][x][y].x = defField[z][x][y].x;
		  defTemp[z][x][y].y = defField[z][x][y].y;
		  defTemp[z][x][y].z = defField[z][x][y].z;
		  }
		  
	
	// calculate jacobian
	for (z=0;z<imageSize.z;z++)
	  for (x=0;x<imageSize.x;x++)
		for (y=0;y<imageSize.y;y++)
		  { // for...for...for...
		  // calculate Jacobian
		  Jacobian[z][x][y]=0.0;
		  
		  if ( x==0 )
				{
				Jaco->data[1][1] = 1 + (defTemp[z][x+1][y].x - defTemp[z][x][y].x);   
				Jaco->data[0][1] = 0 + (defTemp[z][x+1][y].y - defTemp[z][x][y].y);
				Jaco->data[2][1] = 0 + (defTemp[z][x+1][y].z - defTemp[z][x][y].z);	
				}
			else if ( x==(imageSize.x-1) )
				{
				Jaco->data[1][1] = 1 + (defTemp[z][x][y].x - defTemp[z][x-1][y].x);   
				Jaco->data[0][1] = 0 + (defTemp[z][x][y].y - defTemp[z][x-1][y].y);
				Jaco->data[2][1] = 0 + (defTemp[z][x][y].z - defTemp[z][x-1][y].z);
				}
			else 
				{
				Jaco->data[1][1] = 1 + (defTemp[z][x+1][k].x - defTemp[z][x-1][y].x)/2.0;
				Jaco->data[0][1] = 0 + (defTemp[z][x+1][k].y - defTemp[z][x-1][y].y)/2.0;
				Jaco->data[2][1] = 0 + (defTemp[z][x+1][k].z - defTemp[z][x-1][y].z)/2.0;
				}
		   
		   
		   
		   
			if ( y==0 )
				{
				Jaco->data[1][0] = 0 + (defTemp[z][x][y+1].x - defTemp[z][x][y].x);
				Jaco->data[0][0] = 1 + (defTemp[z][x][y+1].y - defTemp[z][x][y].y);
				Jaco->data[2][0] = 0 + (defTemp[z][x][y+1].z - defTemp[z][x][y].z);
				}
			else if ( y==(imageSize.y-1) )
				{
				Jaco->data[1][0] = 0 + (defTemp[z][x][y].x - defTemp[z][x][y-1].x);
				Jaco->data[0][0] = 1 + (defTemp[z][x][y].y - defTemp[z][x][y-1].y);
				Jaco->data[2][0] = 0 + (defTemp[z][x][y].z - defTemp[z][x][y-1].z);
				}
			else
				{
				Jaco->data[1][0] = 0 + (defTemp[z][x][y+1].x - defTemp[z][x][y-1].x)/2.0;
				Jaco->data[0][0] = 1 + (defTemp[z][x][y+1].y - defTemp[z][x][y-1].y)/2.0;
				Jaco->data[2][0] = 0 + (defTemp[z][x][y+1].z - defTemp[z][x][y-1].z)/2.0;
				}
		
		
		   
			if ( z==0 )
				{
				Jaco->data[1][2] = 0 + (defTemp[z+1][x][y].x - defTemp[z][x][y].x);
				Jaco->data[0][2] = 0 + (defTemp[z+1][x][y].y - defTemp[z][x][y].y);
				Jaco->data[2][2] = 1 + (defTemp[z+1][x][y].z - defTemp[z][x][y].z);
				}
			else if ( z==(imageSize.z-1) )
				{
				Jaco->data[1][2] = 0 + (defTemp[z][x][y].x - defTemp[z-1][x][y].x);
				Jaco->data[0][2] = 0 + (defTemp[z][x][y].y - defTemp[z-1][x][y].y);
				Jaco->data[2][2] = 1 + (defTemp[z][x][y].z - defTemp[z-1][x][y].z);
				}
			else
				{
				Jaco->data[1][2] = 0 + (defTemp[z+1][x][y].x - defTemp[z-1][x][y].x)/2.0;
				Jaco->data[0][2] = 0 + (defTemp[z+1][x][y].y - defTemp[z-1][x][y].y)/2.0;
				Jaco->data[2][2] = 1 + (defTemp[z+1][x][y].z - defTemp[z-1][x][y].z)/2.0;
				}
				
				
			Jacobian[z][x][y] = Jaco->data[0][0]*(Jaco->data[1][1]*Jaco->data[2][2]-Jaco->data[2][1]*Jaco->data[1][2]) - Jaco->data[1][0]*(Jaco->data[0][1]*Jaco->data[2][2]-Jaco->data[2][1]*Jaco->data[0][2]) + Jaco->data[2][0]*(Jaco->data[0][1]*Jaco->data[1][2]-Jaco->data[1][1]*Jaco->data[0][2]);
			
			
			if ( Jacobian[z][x][y]<-1.5 || Jacobian[z][x][y]>10.0 )
			    {
				extremeJacobian = YYES;
				
				// method 1: reduce displacement at this voxel at some ratio<1
				// chuzi = 1.0-Jacobian[z][x][y]*(it*5+1)/40.0;
				// //printf("at (%d,%d,%d), jacobian=%f\n",x,y,z,Jacobian[z][x][y]);
				// chuzi = MIN(chuzi, 2.0);
				
				// defField[z][x][y].x = defTemp[z][x][y].x/chuzi;
				// defField[z][x][y].y = defTemp[z][x][y].y/chuzi;
				// defField[z][x][y].z = defTemp[z][x][y].z/chuzi;

				// method 2: smooth displacement at this voxel
				chuzi=1.0;
				for (tz=-1;tz<=1;tz++)
				  for (tx=-1;tx<=1;tx++)
				    for (ty=-1;ty<=1;ty++)
					  {
						xx = MIN(MAX(x+tx, 0), imageSize.x-1);
						yy = MIN(MAX(y+ty, 0), imageSize.y-1);
						zz = MIN(MAX(z+tz, 0), imageSize.z-1);
						
						xxL = MAX(xx-1, 0);
						xxR = MIN(xx+1, imageSize.x-1);
						yyL = MAX(yy-1, 0);
						yyR = MIN(yy+1, imageSize.y-1);
						zzL = MAX(zz-1, 0);
						zzR = MIN(zz+1, imageSize.z-1);
				
				
						defField[zz][xx][yy].x = (0.36* defTemp[zz][xx][yy].x +
									0.11* defTemp[zz][xxL][yy].x  +
									0.11* defTemp[zz][xxR][yy].x  +
									0.11* defTemp[zz][xx][yyL].x  +
									0.11* defTemp[zz][xx][yyR].x  +
									0.05* defTemp[zz][xxL][yyL].x +
									0.05* defTemp[zz][xxL][yyR].x +
									0.05* defTemp[zz][xxR][yyL].x +
									0.05* defTemp[zz][xxR][yyR].x)/chuzi ;
						defField[zz][xx][yy].y = (0.36* defTemp[zz][xx][yy].y +
									0.11* defTemp[zz][xxL][yy].y  +
									0.11* defTemp[zz][xxR][yy].y  +
									0.11* defTemp[zz][xx][yyL].y  +
									0.11* defTemp[zz][xx][yyR].y  +
									0.05* defTemp[zz][xxL][yyL].y +
									0.05* defTemp[zz][xxL][yyR].y +
									0.05* defTemp[zz][xxR][yyL].y +
									0.05* defTemp[zz][xxR][yyR].y)/chuzi ;
						defField[zz][xx][yy].z = (0.36* defTemp[zz][xx][yy].z +
									0.11* defTemp[zz][xxL][yy].z  +
									0.11* defTemp[zz][xxR][yy].z  +
									0.11* defTemp[zz][xx][yyL].z  +
									0.11* defTemp[zz][xx][yyR].z  +
									0.05* defTemp[zz][xxL][yyL].z +
									0.05* defTemp[zz][xxL][yyR].z +
									0.05* defTemp[zz][xxR][yyL].z +
									0.05* defTemp[zz][xxR][yyR].z)/chuzi ;
								
						imageA2BFloat[zz][xx][yy] = trilinearInterpolation(imageA, (float)xx+defField[zz][xx][yy].x, (float)yy+defField[zz][xx][yy].y, (float)zz+defField[zz][xx][yy].z, imageSize);	
					  }
				}
				
			
		  }// for...for...for...
	}
	} // two iterations to remove extreme (negative or positive) jacobians, at a cost of little accuray
		  
	Ffree3d(Jacobian, imageSize.z, imageSize.x);
	FreeMatrix(Jaco);
	Fvector3dfree3d(defTemp, imageSize.z, imageSize.x);
	}//if (encourageDiffeomorphism==YYES)
}

// ---------------------------------------------------------------------------
void GenerateDeformationFieldByFFD2(Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***imageA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels, int iterIdx, int AdditionOrComposition, float *weightsX, float *weightsY, float *weightsZ, int method, int inputInitDeformationOrNot)
{
  int x,y,z;
  float xx,yy,zz;
  int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
  float w1,w2,w3,w4,w5,w6,w7,w8;
  int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
  int l, m,n;
  float dx, dy, dz;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
  int u,v,w;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  int numPoints = numControlPointsX * numControlPointsY * numControlPointsZ;
  Fvector3d ***defFieldTemp;
  float weight, sumWeight;

  // pre-calculate the weights for B-spline, added on June 8, 2009 (begin)
  int maxSizeXY = (int)MAX((float)imageSize.x, (float)imageSize.y);
  int maxSizeXYZ = (int)MAX((float)maxSizeXY, (float)imageSize.z);
 
  if ( ~(levelIndex==numLevels && iterIdx==0) && (AdditionOrComposition==1) )  // 0 for addition and 1 for composition
    defFieldTemp = Fvector3dalloc3d(imageSize.x, imageSize.y, imageSize.z);  
	
  #pragma omp parallel for shared(method,distBetweenControlPointsX,distBetweenControlPointsY,distBetweenControlPointsZ,controlPoints,weightsX,weightsY,weightsZ) private(x,y,z,u,v,w,indexEffectiveControlPointX,indexEffectiveControlPointY,indexEffectiveControlPointZ,indexFirstEffectiveControlPointX,indexFirstEffectiveControlPointY,indexFirstEffectiveControlPointZ,dx,dy,dz,l,m,n,weight) num_threads(100)
  for (z=0;z<imageSize.z;z++)
   { // for each slice
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  { //for each point
	  dx=0.0;
	  dy=0.0;
	  dz=0.0;
		
	  u = (int)floor((float)x/(float)distBetweenControlPointsX);
      v = (int)floor((float)y/(float)distBetweenControlPointsY);
      w = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ);
		
	  if (method==0)  // nearest neighborhood weighting mechanism
		{ // if method==0
		indexEffectiveControlPointX = u;
		indexEffectiveControlPointY = v;
		indexEffectiveControlPointZ = w;
		
        dx = controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
		dy = controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
		dz = controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
		} // if method==0
	  else if (method==1)  // trilinear weighting mechanism
	    { // if method==1
		sumWeight=0.0;
		for (l=0;l<=1;l++)
		  for (m=0;m<=1;m++)
		    for (n=0;n<=1;n++)
			  {
			  indexEffectiveControlPointX = u+l;
			  indexEffectiveControlPointY = v+m;
			  indexEffectiveControlPointZ = w+n;
			  
			  
			  if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				{
				weight = weightsX[x+l*imageSize.x] * weightsY[y+m*imageSize.y] * weightsZ[z+n*imageSize.z];
				sumWeight += weight;
				dx += weight * controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += weight * controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += weight * controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
				}
			  }
		dx /= sumWeight;
		dy /= sumWeight;
		dz /= sumWeight;
		} // if method==1
	  else if (method==2)
		{  // if method==2
		// find out effecting control points
		indexFirstEffectiveControlPointX = u-1;
		indexFirstEffectiveControlPointY = v-1;
		indexFirstEffectiveControlPointZ = w-1;

		// calculate displacement at this voxel (x,y) under the influence of all effective control points
		sumWeight=0.0;
		for (l=0;l<=3;l++)
		 for (m=0;m<=3;m++)
		  for (n=0;n<=3;n++)
		    {
			indexEffectiveControlPointX = l+indexFirstEffectiveControlPointX;
			indexEffectiveControlPointY = m+indexFirstEffectiveControlPointY;
			indexEffectiveControlPointZ = n+indexFirstEffectiveControlPointZ;
			
			if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0 & indexEffectiveControlPointZ<numControlPointsZ) )
			  {
				weight = weightsX[x+l*imageSize.x]*weightsY[y+m*imageSize.y]*weightsZ[z+n*imageSize.z];
				sumWeight += weight;
				dx += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
			  }
			}
			
		dx /= sumWeight;
		dy /= sumWeight;
		dz /= sumWeight;
		} // if method==2
	  
	  
	  // deformation field
	  if ( levelIndex==numLevels && iterIdx==0 && inputInitDeformationOrNot==NNO )  // in the coarsest resolution
		{
		defField[z][x][y].x = dx;
		defField[z][x][y].y = dy;
		defField[z][x][y].z = dz;
		}
      else if ( AdditionOrComposition==0 )  // 0 for addition and 1 for composition
		{ 
		defField[z][x][y].x += dx;
		defField[z][x][y].y += dy;
		defField[z][x][y].z += dz;		  
		}
	  else if ( AdditionOrComposition==1 )  // 0 for addition and 1 for composition
		{
		defFieldTemp[z][x][y].x = defField[z][x][y].x;
		defFieldTemp[z][x][y].y = defField[z][x][y].y;
		defFieldTemp[z][x][y].z = defField[z][x][y].z;
		defField[z][x][y].x = dx;
		defField[z][x][y].y = dy;
		defField[z][x][y].z = dz;
		}		  
		
		
	  } // for each point on this slice
   } // for each slice

   // in case we need to composite deformation field
   if ( ~(levelIndex==numLevels && iterIdx==0) && (AdditionOrComposition==1) ) // composition
     { // if composition
	 printf("\tComposing deformations...\n");

	 #pragma omp parallel for private(x,y,z,xx,yy,zz,xFloor,xCeil,yFloor,yCeil,zFloor,zCeil,dx,dy,dz,w1,w2,w3,w4,w5,w6,w7,w8) num_threads(100)
	 for (z=0;z<imageSize.z;z++)
       { // for each slice
       //printf("\tslice #%d...\n",z);
       for (x=0;x<imageSize.x;x++)
         for (y=0;y<imageSize.y;y++)
		   { // for each point on the slice
		   xx = MAX(MIN((float)x+defField[z][x][y].x, (float)(imageSize.x-1)), 0);
		   yy = MAX(MIN((float)y+defField[z][x][y].y, (float)(imageSize.y-1)), 0);
		   zz = MAX(MIN((float)z+defField[z][x][y].z, (float)(imageSize.z-1)), 0);
		   
		   xFloor = (int)floor(xx);
		   xCeil  = (int)ceil(xx);
		   yFloor = (int)floor(yy);
		   yCeil  = (int)ceil(yy);
	       zFloor = (int)floor(zz);
		   zCeil  = (int)ceil(zz);
		   
		   if ( (xFloor==xCeil)&(xFloor!=(imageSize.x-1)) )  xCeil += 1;
		   if ( (xFloor==xCeil)&(xFloor==(imageSize.x-1)) )  xFloor -= 1;
		   if ( (yFloor==yCeil)&(yFloor!=(imageSize.y-1)) )  yCeil += 1;
		   if ( (yFloor==yCeil)&(yFloor==(imageSize.y-1)) )  yFloor -= 1;
		   if ( (zFloor==zCeil)&(zFloor!=(imageSize.z-1)) )  zCeil += 1;
		   if ( (zFloor==zCeil)&(zFloor==(imageSize.z-1)) )  zFloor -= 1;
		   
		   w1 = (x-xFloor)*(y-yFloor)*(z-zFloor);
		   w2 = (x-xFloor)*(y-yFloor)*(zCeil-z);
		   w3 = (x-xFloor)*(yCeil-y)*(z-zFloor);
		   w4 = (x-xFloor)*(yCeil-y)*(zCeil-z);
		   w5 = (xCeil-x)*(y-yFloor)*(z-zFloor);
		   w6 = (xCeil-x)*(y-yFloor)*(zCeil-z);
		   w7 = (xCeil-x)*(yCeil-y)*(z-zFloor);
		   w8 = (xCeil-x)*(yCeil-y)*(zCeil-z);
		   
		   // incremental displacement at (x,y,z)
		   dx =  w1*defFieldTemp[zCeil][xCeil][yCeil].x
		       + w2*defFieldTemp[zFloor][xCeil][yCeil].x
		       + w3*defFieldTemp[zCeil][xCeil][yFloor].x
			   + w4*defFieldTemp[zFloor][xCeil][yFloor].x
			   + w5*defFieldTemp[zCeil][xFloor][yCeil].x
		       + w6*defFieldTemp[zFloor][xFloor][yCeil].x
		       + w7*defFieldTemp[zCeil][xFloor][yFloor].x
			   + w8*defFieldTemp[zFloor][xFloor][yFloor].x;
			   
		   dy =  w1*defFieldTemp[zCeil][xCeil][yCeil].y
		       + w2*defFieldTemp[zFloor][xCeil][yCeil].y
		       + w3*defFieldTemp[zCeil][xCeil][yFloor].y
			   + w4*defFieldTemp[zFloor][xCeil][yFloor].y
			   + w5*defFieldTemp[zCeil][xFloor][yCeil].y
		       + w6*defFieldTemp[zFloor][xFloor][yCeil].y
		       + w7*defFieldTemp[zCeil][xFloor][yFloor].y
			   + w8*defFieldTemp[zFloor][xFloor][yFloor].y;
			   
		   dz =  w1*defFieldTemp[zCeil][xCeil][yCeil].z
		       + w2*defFieldTemp[zFloor][xCeil][yCeil].z
		       + w3*defFieldTemp[zCeil][xCeil][yFloor].z
			   + w4*defFieldTemp[zFloor][xCeil][yFloor].z
			   + w5*defFieldTemp[zCeil][xFloor][yCeil].z
		       + w6*defFieldTemp[zFloor][xFloor][yCeil].z
		       + w7*defFieldTemp[zCeil][xFloor][yFloor].z
			   + w8*defFieldTemp[zFloor][xFloor][yFloor].z;
			
		   defField[z][x][y].x += dx;
		   defField[z][x][y].y += dy;
		   defField[z][x][y].z += dz;
   
		   } // for each point on the slice
	   }// for each slice
	 }  // if composition
	 
   
   
   if ( ~(levelIndex==numLevels && iterIdx==0) && (AdditionOrComposition==1) )  // 0 for addition and 1 for composition
     Fvector3dfree3d(defFieldTemp, imageSize.z, imageSize.x);
	
}

// ---------------------------------------------------------------------------
float GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy(unsigned char ****featureMapA, unsigned char ****featureMapB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, float ***confidenceMap, int numFeatures, float ****featureMapA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels)
{
  int x,y,z;
  int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
  int l,m,n;
  int featureIndex;
  float feature;
  float dx,dy,dz;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
  float u,v,w;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  
  float totalEnergy=0.0;
  float energySquaredAtThisPoint;
  
 
  for (z=0;z<imageSize.z;z++)
   { // for each slice
   printf("\tslice #%d...\n",z);
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  { //for each point
		dx=0;
		dy=0;
		dz=0;
		
		energySquaredAtThisPoint = 0.0;
		
		// find out effecting control points
		indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPointsX)-1;
		indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPointsY)-1;
		indexFirstEffectiveControlPointZ = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ)-1;
		u = (float)x/(float)distBetweenControlPointsX - floor((float)x/(float)distBetweenControlPointsX);
		v = (float)y/(float)distBetweenControlPointsY - floor((float)y/(float)distBetweenControlPointsY);
		w = (float)(z-1)/(float)distBetweenControlPointsZ - floor((float)(z-1)/(float)distBetweenControlPointsZ);
	    
		
		// calculate displacement at this voxel (x,y) under the influence of all effective control points
		for (l=0;l<=3;l++)
		 for (m=0;m<=3;m++)
		  for (n=0;n<=3;n++)
		    {
			indexEffectiveControlPointX = l+indexFirstEffectiveControlPointX;
			indexEffectiveControlPointY = m+indexFirstEffectiveControlPointY;
			indexEffectiveControlPointZ = n+indexFirstEffectiveControlPointZ;
			
			if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0 & indexEffectiveControlPointZ<numControlPointsZ) )
			  {
			    dx += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
			  }
			}
		
		if (levelIndex==numLevels)  // in the coarsest resolution
		  {
		  defField[z][x][y].x = dx;
		  defField[z][x][y].y = dy;
		  defField[z][x][y].z = dz;
		  }
		else
		  { 
		  defField[z][x][y].x += dx;
		  defField[z][x][y].y += dy;
		  defField[z][x][y].z += dz;
		  }		
		
		//interpolate features (float)
		for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
		  {
		  feature = trilinearInterpolation(featureMapA[featureIndex], (float)x+defField[z][x][y].x, (float)y+defField[z][x][y].y, (float)z+defField[z][x][y].z, imageSize);

		  featureMapA2BFloat[featureIndex][z][x][y] = feature;
		  energySquaredAtThisPoint += pow( (featureMapA2BFloat[featureIndex][z][x][y] - (float)featureMapB[featureIndex][z][x][y]), 2.0 );
		  }
		totalEnergy += energySquaredAtThisPoint * confidenceMap[z][x][y];   // confidence map involved
	  } // for each point on this slice
  } // for each slice
	  
  // smooth (twice) deformation in each direction
  smoothDeformationField(defField, defField, imageSize, levelIndex, numLevels);
  totalEnergy /= (float)numFeatures;
  
  return totalEnergy;
}

// ---------------------------------------------------------------------------
// speedup version of sub-function GenerateFeaturesAndDeformationFieldByFFDAndCalculateEnergy
void GenerateFeaturesAndDeformationFieldByFFD(unsigned char ****featureMapA, unsigned char ****featureMapB, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, int numFeatures, float ****featureMapA2BFloat, Fvector3d ***defField, int levelIndex, int numLevels)
{
  int x,y,z;
  int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
  int l,m,n;
  int featureIndex;
  float feature;
  float dx,dy,dz;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
  float u,v,w;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  
  int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
  float x_src, y_src, z_src;
  float w1, w2, w3, w4, w5, w6, w7, w8;
  float distXtoXFloor,distXtoXCeil,distYtoYFloor,distYtoYCeil,distZtoZFloor,distZtoZCeil;
 
  // pre-calculate the weights for B-spline, added on June 8, 2009 (begin)
  int maxSizeXY = (int)MAX((float)imageSize.x, (float)imageSize.y);
  int maxSizeXYZ = (int)MAX((float)maxSizeXY, (float)imageSize.z);
  float preWeightsXY[maxSizeXY][4];
  float preWeightsZ[imageSize.z][4];
  int int1, int2;
  float weight;
  for (int1=0;int1<maxSizeXYZ;int1++)
    for (int2=0;int2<4;int2++)
	  {
	  if (int1<maxSizeXY)
	    preWeightsXY[int1][int2] = Bspline(  ( (float)int1/(float)distBetweenControlPointsX-floor((float)int1/(float)distBetweenControlPointsY) ), int2);
	  if (int1<imageSize.z)
	    preWeightsZ[int1][int2] = Bspline( ((float)(int1-1)/(float)distBetweenControlPointsZ - floor((float)(int1-1)/(float)distBetweenControlPointsZ) ), int2);
	  }
	  
  #pragma omp parallel for private(x,y,z,dx,dy,dz,indexFirstEffectiveControlPointX,indexFirstEffectiveControlPointY,indexFirstEffectiveControlPointZ,l,m,n,indexEffectiveControlPointX,indexEffectiveControlPointY,indexEffectiveControlPointZ,weight,xFloor,xCeil,yFloor,yCeil,zFloor,zCeil,x_src,y_src,z_src,distXtoXFloor,distXtoXCeil,distYtoYFloor,distYtoYCeil,distZtoZFloor,distZtoZCeil,w1,w2,w3,w4,w5,w6,w7,w8,featureIndex,feature) num_threads(100)
  for (z=0;z<imageSize.z;z++)
   { // for each slice
   printf("\tslice #%d...\n",z);
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  { //for each point
		dx=0;
		dy=0;
		dz=0;

		// find out effecting control points
		indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPointsX)-1;
		indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPointsY)-1;
		indexFirstEffectiveControlPointZ = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ)-1;

		// calculate displacement at this voxel (x,y) under the influence of all effective control points
		for (l=0;l<=3;l++)
		 for (m=0;m<=3;m++)
		  for (n=0;n<=3;n++)
		    {
			indexEffectiveControlPointX = l+indexFirstEffectiveControlPointX;
			indexEffectiveControlPointY = m+indexFirstEffectiveControlPointY;
			indexEffectiveControlPointZ = n+indexFirstEffectiveControlPointZ;
			
			//printf("index = (%d, %d, %d)\n", indexEffectiveControlPointX,indexEffectiveControlPointY,indexEffectiveControlPointZ);
			if ( (indexEffectiveControlPointX>=0 & indexEffectiveControlPointX<numControlPointsX)&&(indexEffectiveControlPointY>=0 & indexEffectiveControlPointY<numControlPointsY)&&(indexEffectiveControlPointZ>=0 & indexEffectiveControlPointZ<numControlPointsZ) )
			  {
			    weight = preWeightsXY[x][l]*preWeightsXY[y][m]*preWeightsZ[z][n];
				dx += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].x;
				dy += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].y;
				dz += weight*controlPoints[indexEffectiveControlPointZ][indexEffectiveControlPointX][indexEffectiveControlPointY].z;
			  }
			}
		
		//save deformation field
		if (levelIndex==numLevels)  // in the coarsest resolution
		  {
		  defField[z][x][y].x = dx;
		  defField[z][x][y].y = dy;
		  defField[z][x][y].z = dz;
		  }
		else
		  { 
		  defField[z][x][y].x += dx;
		  defField[z][x][y].y += dy;
		  defField[z][x][y].z += dz;
		  }

		//interpolate features (float)
		x_src = (float)x+defField[z][x][y].x;
		y_src = (float)y+defField[z][x][y].y;
		z_src = (float)z+defField[z][x][y].z;
	
		if (x_src<0)  x_src=0.0;
		if (y_src<0)  y_src=0.0;
		if (z_src<0)  z_src=0.0;
		if (x_src>imageSize.x-1) x_src=(float)(imageSize.x-1);
		if (y_src>imageSize.y-1) y_src=(float)(imageSize.y-1);
		if (z_src>imageSize.z-1) z_src=(float)(imageSize.z-1);
	
		xFloor = (int)floor(x_src);
		xCeil  = (int)ceil(x_src);
		yFloor = (int)floor(y_src);
		yCeil  = (int)ceil(y_src);
		zFloor = (int)floor(z_src);
		zCeil  = (int)ceil(z_src);
	
		if ( (xFloor==xCeil)&(xFloor!=(imageSize.x-1)) )  xCeil+=1;
		if ( (xFloor==xCeil)&(xFloor==(imageSize.x-1)) )  xFloor-=1;
		if ( (yFloor==yCeil)&(yFloor!=(imageSize.y-1)) )  yCeil+=1;
		if ( (yFloor==yCeil)&(yFloor==(imageSize.y-1)) )  yFloor-=1;
		if ( (zFloor==zCeil)&(zFloor!=(imageSize.z-1)) )  zCeil+=1;
		if ( (zFloor==zCeil)&(zFloor==(imageSize.z-1)) )  zFloor-=1;
	
		distXtoXFloor = x_src-xFloor;
		distXtoXCeil  = xCeil-x_src;
		distYtoYFloor = y_src-yFloor;
		distYtoYCeil  = yCeil-y_src;
		distZtoZFloor = z_src-zFloor;
		distZtoZCeil  = zCeil-z_src;
	
		w1 = distXtoXFloor*distYtoYFloor*distZtoZFloor;
		w2 = distXtoXFloor*distYtoYFloor*distZtoZCeil;
		w3 = distXtoXFloor*distYtoYCeil*distZtoZFloor;
		w4 = distXtoXFloor*distYtoYCeil*distZtoZCeil;
		w5 = distXtoXCeil*distYtoYFloor*distZtoZFloor;
		w6 = distXtoXCeil*distYtoYFloor*distZtoZCeil;
		w7 = distXtoXCeil*distYtoYCeil*distZtoZFloor;
		w8 = distXtoXCeil*distYtoYCeil*distZtoZCeil;
		
		for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
		  {
		  feature = (float)featureMapA[featureIndex][zCeil][xCeil][yCeil]   * w1 +
					(float)featureMapA[featureIndex][zFloor][xCeil][yCeil]  * w2 +
					(float)featureMapA[featureIndex][zCeil][xCeil][yFloor]  * w3 +
					(float)featureMapA[featureIndex][zFloor][xCeil][yFloor] * w4 +
					(float)featureMapA[featureIndex][zCeil][xFloor][yCeil]  * w5 +
					(float)featureMapA[featureIndex][zFloor][xFloor][yCeil] * w6 +
					(float)featureMapA[featureIndex][zCeil][xFloor][yFloor] * w7 +
					(float)featureMapA[featureIndex][zFloor][xFloor][yFloor]* w8 ;
		  if ( (feature==0.0)&&(x_src<0 || x_src>(imageSize.x-1) || y_src<0 || y_src>(imageSize.y-1) || z_src<0 || z_src>(imageSize.z-1)) )
		    feature = MAX(0.0, (float)featureMapA[featureIndex][z][x][y]);
			
		  featureMapA2BFloat[featureIndex][z][x][y] = feature;
		  }
	  } // for each point on this slice
  } // for each slice
}

// ---------------------------------------------------------------------------
float GenerateFeaturesFromDeformationFieldAndCalculateEnergy(unsigned char ****featureMapA, unsigned char ****featureMapB, float ****featureMapA2BFloat, Fvector3d ***defField, float ***confidenceMap, Ivector3d imageSize, int numFeatures)
{
  int x,y,z;
  int featureIndex;
  float totalEnergy=0.0;
  float energySquaredAtThisPoint;
  
  for (z=0;z<imageSize.z;z++)
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  { //for each point
	  energySquaredAtThisPoint = 0.0;
	  
	  for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
	    {
		featureMapA2BFloat[featureIndex][z][x][y] = trilinearInterpolation(featureMapA[featureIndex], (float)x+defField[z][x][y].x, (float)y+defField[z][x][y].y, (float)z+defField[z][x][y].z, imageSize);
		energySquaredAtThisPoint += pow( (featureMapA2BFloat[featureIndex][z][x][y] - (float)featureMapB[featureIndex][z][x][y]), 2.0 );
		}  // for each feature at this voxel
	  
	  totalEnergy += energySquaredAtThisPoint * confidenceMap[z][x][y];   // confidence map involved
	  }// for each voxel
	 
  totalEnergy /= (float)numFeatures;
  
  return totalEnergy;
}

// ---------------------------------------------------------------------------
void GenerateFloatImageFromDeformationField(unsigned char ***imageA, unsigned char ***imageB, float ***imageA2BFloat, Fvector3d ***defField, Ivector3d imageSize)
{
  int x,y,z;
  
  for (z=0;z<imageSize.z;z++)
   for (x=0;x<imageSize.x;x++)
    for (y=0;y<imageSize.y;y++)
	  {
	  imageA2BFloat[z][x][y] = trilinearInterpolation(imageA, (float)x+defField[z][x][y].x, (float)y+defField[z][x][y].y, (float)z+defField[z][x][y].z, imageSize);
	  }
	  
}

// ---------------------------------------------------------------------------
float FFDLocalEffect(unsigned char ***imageA, int xcoor, int ycoor, int zcoor, Fvector3d ***controlPoints, Fvector3d ***defField, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, float increment, char orientation, int tempOrNot, int featureIndex, int numFeatures)
{
  float intensity;
  int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
  int l, m,n;
  float dx, dy, dz;
  int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
  float u,v,w;
  int numControlPointsX = (int)ceil((float)imageSize.x/distBetweenControlPointsX);
  int numControlPointsY = (int)ceil((float)imageSize.y/distBetweenControlPointsY);
  int numControlPointsZ = (int)ceil((float)(imageSize.z-1)/distBetweenControlPointsZ);
  
  //calculate displacement at this voxel (xcoor, ycoor) under the influence of all effective control points
  dx=defField[zcoor][xcoor][ycoor].x;
  dy=defField[zcoor][xcoor][ycoor].y;
  dz=defField[zcoor][xcoor][ycoor].z;
 
  // find out effecting control points
  indexFirstEffectiveControlPointX = (int)floor((float)xcoor/(float)distBetweenControlPointsX)-1;
  indexFirstEffectiveControlPointY = (int)floor((float)ycoor/(float)distBetweenControlPointsY)-1;
  indexFirstEffectiveControlPointZ = (int)floor((float)(zcoor-1)/(float)distBetweenControlPointsZ)-1;
  u = (float)xcoor/(float)distBetweenControlPointsX - floor((float)xcoor/(float)distBetweenControlPointsX);
  v = (float)ycoor/(float)distBetweenControlPointsY - floor((float)ycoor/(float)distBetweenControlPointsY);
  w = (float)(zcoor-1)/(float)distBetweenControlPointsZ - floor((float)(zcoor-1)/(float)distBetweenControlPointsZ);
  
	
  l = controlPointIndex1-indexFirstEffectiveControlPointX;
  m = controlPointIndex2-indexFirstEffectiveControlPointY;
  n = controlPointIndex3-indexFirstEffectiveControlPointZ;
 

  switch (orientation)
    {
	case 'X':  
	  dx += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*increment;
	  break;
	  
	case 'Y':  
	  dy += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*increment;
	  break;
	
	case 'Z':
	  dz += Bspline(u,l)*Bspline(v,m)*Bspline(w,n)*increment;
	  break;
	  
	default:
	  break;
	}
  
  
  if (tempOrNot==NNO && featureIndex==(numFeatures-1))  // 1) not temp, need to update defField 2) update deformation at this point only once, not numFeatures times !!!
    {
	defField[zcoor][xcoor][ycoor].x = dx;
	defField[zcoor][xcoor][ycoor].y = dy;
	defField[zcoor][xcoor][ycoor].z = dz;
	}
 
  //interpolate intensity (float)
  intensity = trilinearInterpolation(imageA, (float)xcoor+dx, (float)ycoor+dy, (float)zcoor+dz, imageSize);

  return intensity;
}

// ---------------------------------------------------------------------------
void CopyFloatImage3D(float ***imgA, float ***imgB, Ivector3d imgSize)    // copy A to B
{
  int i,j,k;
  
  for (k=0;k<imgSize.z;k++)
   for (i=0;i<imgSize.x;i++)
    for (j=0;j<imgSize.y;j++)
	  {
	  imgB[k][i][j] = imgA[k][i][j];
	  }
}

// ---------------------------------------------------------------------------
void CopyFeatureMapsFloat(float ****featureMapOrig, float ****featureMapCopy, Ivector3d imageSize, int numFeatures)
{
  int i,j,k;
  int featureIndex;
  
  for (k=0;k<imageSize.z;k++)
   for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
	  for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
	    {
		featureMapCopy[featureIndex][k][i][j] = featureMapOrig[featureIndex][k][i][j];
		}
}

// ---------------------------------------------------------------------------
float sign(float x)
{
  float signX;
  signX = (float)((x>0)-(x<0));
  return signX;
}

// ---------------------------------------------------------------------------
void CopyDisplacementAtControlPoints3D(Fvector3d ***controlPoints, Fvector3d ***controlPointsBackup, int numControlPointsX, int numControlPointsY, int numControlPointsZ)
{
  int i,j,k;
  for (k=0;k<numControlPointsZ;k++)
   for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
	  controlPointsBackup[k][i][j].x=controlPoints[k][i][j].x;
	  controlPointsBackup[k][i][j].y=controlPoints[k][i][j].y;
	  controlPointsBackup[k][i][j].z=controlPoints[k][i][j].z;
	  }
}

// ---------------------------------------------------------------------------
void CopyDeformationField(Fvector3d ***def, Fvector3d ***defCopy, Ivector3d imageSize)
{
  int i,j,k;
  for (k=0;k<imageSize.z;k++)
   for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
	  {
	  defCopy[k][i][j].x=def[k][i][j].x;
	  defCopy[k][i][j].y=def[k][i][j].y;
	  defCopy[k][i][j].z=def[k][i][j].z;
	  }
}

// ---------------------------------------------------------------------------
void CopyConfidenceMap(float ***MSMap, float ***MSMapCopy, Ivector3d imageSize)
{
  int i,j,k;
  for (k=0;k<imageSize.z;k++)
   for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
	  {
	  MSMapCopy[k][i][j]=MSMap[k][i][j];
	  MSMapCopy[k][i][j]=MSMap[k][i][j];
	  MSMapCopy[k][i][j]=MSMap[k][i][j];
	  }
}

// ---------------------------------------------------------------------------
void calculateDisplacementIncrementAtControlPoints(Fvector3d ***controlPointsBackup, Fvector3d ***controlPoints, Fvector3d ***controlPointsIncrement, int numControlPointsX, int numControlPointsY, int numControlPointsZ)
{
  int i,j,k;
  for (k=0;k<numControlPointsZ;k++)
   for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
	  controlPointsIncrement[k][i][j].x = controlPoints[k][i][j].x - controlPointsBackup[k][i][j].x;
	  controlPointsIncrement[k][i][j].y = controlPoints[k][i][j].y - controlPointsBackup[k][i][j].y;
	  controlPointsIncrement[k][i][j].z = controlPoints[k][i][j].z - controlPointsBackup[k][i][j].z;
	  }
}

// ---------------------------------------------------------------------------
void smoothDisplacementAtControlPoints3D(Fvector3d ***controlPoints, Fvector3d ***controlPointsSmoothed, int numControlPointsX, int numControlPointsY, int numControlPointsZ, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int levelIndex, int numLevels, float factor)
{
  int i,j,k;
  float ***smoothMat;
  float gaussianKernelXY, gaussianKernelZ;
  int gaussianRadiusXY, gaussianRadiusZ;
  float ***controlPointsX, ***controlPointsY, ***controlPointsZ;
  float ***controlPointsXSmoothed, ***controlPointsYSmoothed, ***controlPointsZSmoothed;
  Ivector3d inputSize, smoothMatSize;

  // determine gaussian kernel and gaussian radius in xy plane
  if (numControlPointsX<=7 && numControlPointsX*distBetweenControlPointsX<80)
   {
   gaussianRadiusXY = 2;
   gaussianKernelXY = 1.4 * factor;
   gaussianRadiusZ = 1;
   gaussianKernelZ = 0.6 * factor; 
   }
  else if (numControlPointsX<=15 && numControlPointsX*distBetweenControlPointsX<120)
   {
   gaussianRadiusXY = 2;
   gaussianKernelXY = 2.0 * factor;
   gaussianRadiusZ = 1;
   gaussianKernelZ = 0.8*factor;
   }
  else if (numControlPointsX>15 || numControlPointsX*distBetweenControlPointsX>200)
   {
   gaussianRadiusXY = 4;
   gaussianKernelXY = 3.8 * factor;
   gaussianRadiusZ = 2;
   gaussianKernelZ = 1.6 * factor;
   }
  else
   {
    gaussianRadiusXY = 4;
    gaussianKernelXY = 4.2 * factor;  // the smoothing kernel is adaptive to the distance between two adjacent control points  
	gaussianRadiusZ = 2;
	gaussianKernelZ = 1.6 * factor;
    if (gaussianKernelXY<0.35)
      gaussianKernelXY = 0.35;  // if kernel is less than 0.35, the gaussian template will be 1 at the center but 0 everywhere else. This should generally be avoided.
   }
  
  // determine guassian kernel and gaussian radius in z direction
  inputSize.x = numControlPointsX+2*gaussianRadiusXY;
  inputSize.y = numControlPointsY+2*gaussianRadiusXY;
  inputSize.z = numControlPointsZ+2*gaussianRadiusZ;
  smoothMatSize.x = 2*gaussianRadiusXY+1;
  smoothMatSize.y = 2*gaussianRadiusXY+1;
  smoothMatSize.z = 2*gaussianRadiusZ+1;
  
  smoothMat = Falloc3d(smoothMatSize.x, smoothMatSize.y, smoothMatSize.z);
  gauss3D(smoothMatSize.x, gaussianKernelXY, smoothMatSize.y, gaussianKernelXY, smoothMatSize.z, gaussianKernelZ, smoothMat);

  controlPointsX = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  controlPointsY = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  controlPointsZ = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  controlPointsXSmoothed = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  controlPointsYSmoothed = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  controlPointsZSmoothed = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  
  for (k=0;k<inputSize.z;k++)
   for (i=0;i<inputSize.x;i++)
    for (j=0;j<inputSize.y;j++)
	  {
	  if ( (i<gaussianRadiusXY) | (j<gaussianRadiusXY) | (k<gaussianRadiusZ) | (i>=(numControlPointsX+gaussianRadiusXY)) | (j>=(numControlPointsY+gaussianRadiusXY)) | (k>=(numControlPointsZ+gaussianRadiusZ)) )
	    {
		controlPointsX[k][i][j] = 0.0;
		controlPointsY[k][i][j] = 0.0;
		controlPointsZ[k][i][j] = 0.0;
		controlPointsXSmoothed[k][i][j] = 0.0;
		controlPointsYSmoothed[k][i][j] = 0.0;
		controlPointsYSmoothed[k][i][j] = 0.0;
		}
	  else
	    {
		controlPointsX[k][i][j] = controlPoints[k-gaussianRadiusZ][i-gaussianRadiusXY][j-gaussianRadiusXY].x;
		controlPointsY[k][i][j] = controlPoints[k-gaussianRadiusZ][i-gaussianRadiusXY][j-gaussianRadiusXY].y;
		controlPointsZ[k][i][j] = controlPoints[k-gaussianRadiusZ][i-gaussianRadiusXY][j-gaussianRadiusXY].z;
		controlPointsXSmoothed[k][i][j] = 0.0;
		controlPointsYSmoothed[k][i][j] = 0.0;
		controlPointsZSmoothed[k][i][j] = 0.0;
	    }
	  }
  Convolution3DFloat(controlPointsX, inputSize, smoothMat, smoothMatSize, controlPointsXSmoothed);
  Convolution3DFloat(controlPointsY, inputSize, smoothMat, smoothMatSize, controlPointsYSmoothed);
  Convolution3DFloat(controlPointsZ, inputSize, smoothMat, smoothMatSize, controlPointsZSmoothed);
  
  for (k=0;k<numControlPointsZ;k++)
   for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
		controlPointsSmoothed[k][i][j].x = controlPointsXSmoothed[k+gaussianRadiusZ][i+gaussianRadiusXY][j+gaussianRadiusXY];
		controlPointsSmoothed[k][i][j].y = controlPointsYSmoothed[k+gaussianRadiusZ][i+gaussianRadiusXY][j+gaussianRadiusXY];
		controlPointsSmoothed[k][i][j].z = controlPointsZSmoothed[k+gaussianRadiusZ][i+gaussianRadiusXY][j+gaussianRadiusXY];
	  }
	  
  // release momery
  Ffree3d(controlPointsX, inputSize.z, inputSize.x);
  Ffree3d(controlPointsY, inputSize.z, inputSize.x);
  Ffree3d(controlPointsZ, inputSize.z, inputSize.x);
  Ffree3d(controlPointsXSmoothed, inputSize.z, inputSize.x);
  Ffree3d(controlPointsYSmoothed, inputSize.z, inputSize.x);
  Ffree3d(controlPointsZSmoothed, inputSize.z, inputSize.x);
  Ffree3d(smoothMat, smoothMatSize.z, smoothMatSize.x);
}

// ---------------------------------------------------------------------------
void UpdateControlPointsWithSmoothIncrement2(Fvector3d ***controlPointsBackup, Fvector3d ***controlPointsIncrement, Fvector3d ***controlPointsIncrementSmoothed, Fvector3d ***controlPoints, int numControlPointsX, int numControlPointsY, int numControlPointsZ, float alpha, unsigned char ***mask, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int levelIndex, int numLevels)
{
  int i,j,k, ii,jj,kk;
  Ivector3d realCoordinate;
  float maxDisplacementAllowedX, maxDisplacementAllowedY, maxDisplacementAllowedZ;
  float offset;
  
  if (distBetweenControlPointsX<14)
    maxDisplacementAllowedX = 1.05*distBetweenControlPointsX * pow(0.95, numLevels+1-levelIndex);
  else
    maxDisplacementAllowedX = 0.9*distBetweenControlPointsX * pow(0.90, numLevels-levelIndex);
	
  if (distBetweenControlPointsY<14)
    maxDisplacementAllowedY = 1.05*distBetweenControlPointsY * pow(0.95, numLevels+1-levelIndex);
  else
    maxDisplacementAllowedY = 0.9*distBetweenControlPointsY * pow(0.90, numLevels-levelIndex);
	 
  if (distBetweenControlPointsZ<4)
    maxDisplacementAllowedZ = 0.5*distBetweenControlPointsZ * pow(0.95, numLevels+1-levelIndex);
  else
    maxDisplacementAllowedZ = 0.65*distBetweenControlPointsZ * pow(0.90, numLevels-levelIndex);
	
	
  printf("alpha = %f, maxum displacement allowed = (%f, %f, %f)\n\n", alpha, maxDisplacementAllowedX, maxDisplacementAllowedY, maxDisplacementAllowedZ);
  
  for (k=0;k<numControlPointsZ;k++)
   for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
	  realCoordinate.x = i*distBetweenControlPointsX;
	  realCoordinate.y = j*distBetweenControlPointsY;
	  realCoordinate.z = k*distBetweenControlPointsZ;
	  
	  controlPoints[k][i][j].x = controlPointsBackup[k][i][j].x + ((1-alpha)*controlPointsIncrement[k][i][j].x + alpha*controlPointsIncrementSmoothed[k][i][j].x);
	  controlPoints[k][i][j].y = controlPointsBackup[k][i][j].y + ((1-alpha)*controlPointsIncrement[k][i][j].y + alpha*controlPointsIncrementSmoothed[k][i][j].y);
	  controlPoints[k][i][j].z = controlPointsBackup[k][i][j].z + ((1-alpha)*controlPointsIncrement[k][i][j].z + alpha*controlPointsIncrementSmoothed[k][i][j].z);
	  }
	
  for (k=0;k<numControlPointsZ;k++)
   for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
	  // last step: check if displacement at a control point has exceeded the maximum allowed
	  if (fabs(controlPoints[k][i][j].x) > maxDisplacementAllowedX)
	    {
		offset = controlPoints[k][i][j].x - sign(controlPoints[k][i][j].x) * maxDisplacementAllowedX;
		controlPoints[k][i][j].x = sign(controlPoints[k][i][j].x) * maxDisplacementAllowedX;
		
		ii = i+(int)sign(controlPoints[k][i][j].x);
		while (ii>=0 && ii<numControlPointsX)
		  {
		  controlPoints[k][ii][j].x += offset*pow(0.8, fabs((float)(ii-i)));
		  if (fabs(controlPoints[k][ii][j].x) > maxDisplacementAllowedX)
		    {
			offset = controlPoints[k][ii][j].x - sign(controlPoints[k][ii][j].x)*maxDisplacementAllowedX;
		    controlPoints[k][ii][j].x = sign(controlPoints[k][ii][j].x) * maxDisplacementAllowedX;
			ii += (int)sign(controlPoints[k][ii][j].x);
			}
		  else 
		    break;
		  }//while
		} // if
	  if (fabs(controlPoints[k][i][j].y) > maxDisplacementAllowedY)
	    {
		offset = controlPoints[k][i][j].y - sign(controlPoints[k][i][j].y)*maxDisplacementAllowedY;
		controlPoints[k][i][j].y = sign(controlPoints[k][i][j].y) * maxDisplacementAllowedY;
				
		jj = j+(int)sign(controlPoints[k][i][j].y);
		while (jj>=0 && jj<numControlPointsY)
		  {
		  controlPoints[k][i][jj].y += offset*pow(0.8, fabs((float)(jj-j)));
		  if (fabs(controlPoints[k][i][jj].y) > maxDisplacementAllowedY)
		    {
			offset = controlPoints[k][i][jj].y - sign(controlPoints[k][i][jj].y)*maxDisplacementAllowedY;
		    controlPoints[k][i][jj].y = sign(controlPoints[k][i][jj].y) * maxDisplacementAllowedY;
			jj += (int)sign(controlPoints[k][i][jj].y);
			}
		  else 
		    break;
		  }//while
		} //if
	  if (fabs(controlPoints[k][i][j].z) > maxDisplacementAllowedZ)
	    {
		offset = controlPoints[k][i][j].z - sign(controlPoints[k][i][j].z)*maxDisplacementAllowedZ;
		controlPoints[k][i][j].z = sign(controlPoints[k][i][j].z) * maxDisplacementAllowedZ;
				
		kk = k+(int)sign(controlPoints[k][i][j].z);
		while (kk>=0 && kk<numControlPointsZ)
		  {
		  controlPoints[kk][i][j].z += offset*pow(0.8, fabs((float)(kk-k)));
		  if (fabs(controlPoints[kk][i][j].z) > maxDisplacementAllowedZ)
		    {
			offset = controlPoints[kk][i][j].z - sign(controlPoints[kk][i][j].z)*maxDisplacementAllowedZ;
		    controlPoints[kk][i][j].z = sign(controlPoints[kk][i][j].z) * maxDisplacementAllowedZ;
			kk += (int)sign(controlPoints[kk][i][j].z);
			}
		  else 
		    break;
		  }//while
		} //if
	  }
}

// ---------------------------------------------------------------------------
void UpdateControlPointsWithSmoothIncrement(Fvector3d ***controlPointsBackup, Fvector3d ***controlPoints, Fvector3d ***controlPointsUpdated, int numControlPointsX, int numControlPointsY, int numControlPointsZ, int IterIndex, int maxNumIterInResolution, int levelIndex, int numLevels, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, unsigned char ***mask)
{
  // idea: controlPointsUpdated = controlPointsBackup + (1-alpha)*(controlPoints-controlPointsBackup) + alpha*smooth(controlPoints-controlPointsBackup)
  Fvector3d ***controlPointsIncrement, ***controlPointsIncrementSmoothed;
  controlPointsIncrement = Fvector3dalloc3d(numControlPointsX, numControlPointsY, numControlPointsZ);
  controlPointsIncrementSmoothed = Fvector3dalloc3d(numControlPointsX, numControlPointsY, numControlPointsZ);
  float alpha;   // smooth weighting
  
  calculateDisplacementIncrementAtControlPoints(controlPointsBackup, controlPoints, controlPointsIncrement, numControlPointsX, numControlPointsY, numControlPointsZ);  // calculate increment: controlPoints - controlPointsBackup
  smoothDisplacementAtControlPoints3D(controlPointsIncrement, controlPointsIncrementSmoothed, numControlPointsX, numControlPointsY, numControlPointsZ, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, levelIndex, numLevels, 1.0); // smooth: controlPointsIncrement -> controlPointsIncrementSmoothed
  alpha=0.55*exp(-pow(IterIndex,2.0)/(2.0*pow(maxNumIterInResolution,2.0)))*pow(0.8, (levelIndex-1)); // Parameter "alpha" controls how much to smooth: the bigger alpha is, the more smoothness there will be. At coarse stage, smooth more; as displacement gets more and more accurate at each control point, smooth less
  UpdateControlPointsWithSmoothIncrement2(controlPointsBackup, controlPointsIncrement, controlPointsIncrementSmoothed, controlPointsUpdated, numControlPointsX, numControlPointsY, numControlPointsZ, alpha, mask, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, levelIndex, numLevels);
  
  //release memory
  Fvector3dfree3d(controlPointsIncrement, numControlPointsZ, numControlPointsX);
  Fvector3dfree3d(controlPointsIncrementSmoothed, numControlPointsZ, numControlPointsX);
}

// ---------------------------------------------------------------------------
void LinearCombinationOfTwoControlPointsMats(Fvector3d ***controlPointsNew, Fvector3d ***controlPointsOld, Fvector3d ***controlPointsUpdated, float weighting, int numControlPointsX, int numControlPointsY, int numControlPointsZ)
{
  int i,j,k;
  
  for (k=0;k<numControlPointsZ;k++)
   for (i=0;i<numControlPointsX;i++)
    for (j=0;j<numControlPointsY;j++)
	  {
	  controlPointsUpdated[k][i][j].x = weighting*controlPointsOld[k][i][j].x + (1-weighting)*controlPointsNew[k][i][j].x;
	  controlPointsUpdated[k][i][j].y = weighting*controlPointsOld[k][i][j].y + (1-weighting)*controlPointsNew[k][i][j].y;
	  controlPointsUpdated[k][i][j].z = weighting*controlPointsOld[k][i][j].z + (1-weighting)*controlPointsNew[k][i][j].z;
	  }
}

// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
void Convolution3DFloat(float*** input,Ivector3d inputSize, float*** mask,Ivector3d maskSize, float*** output)
{
  int row,col,z, ii,jj,kk, iii,jjj,kkk;
  Ivector3d halfSize,evenOrNot;
  double sum;

  halfSize.x=(maskSize.x-1)/2;
  halfSize.y=(maskSize.y-1)/2;
  halfSize.z=(maskSize.z-1)/2;

  evenOrNot.x=-(maskSize.x-1)%2;
  evenOrNot.y=-(maskSize.y-1)%2;
  evenOrNot.z=-(maskSize.z-1)%2;

  for (z=0;z<inputSize.z;z++)
   for (row=0;row<inputSize.x;row++)
    for (col=0;col<inputSize.y;col++)
      {
      sum=0;  // initialize
      
      for (kk=-halfSize.z-evenOrNot.z; kk<=halfSize.z; kk++)
	{
	kkk = z+kk;
	if (kkk<0)  			kkk = - kkk -1;
	else if (kkk>=inputSize.z)	kkk = (inputSize.z-1) - (kkk-inputSize.z);

	for (ii=-halfSize.x-evenOrNot.x; ii<=halfSize.x; ii++)
	  {
	  iii = row+ii;
	  if (iii<0)			iii = - iii - 1;
	  else if (iii>=inputSize.x)	iii = (inputSize.x-1) - (iii-inputSize.x);

	  for (jj=-halfSize.y-evenOrNot.y; jj<=halfSize.y; jj++)
	    {
	    jjj = col+jj;
	    if (jjj<0)			jjj = - jjj - 1;
	    else if (jjj>=inputSize.y)	jjj = (inputSize.y-1) - (jjj-inputSize.y);

	    sum += (float)input[kkk][iii][jjj] * mask[kk+halfSize.z+evenOrNot.z][ii+halfSize.x+evenOrNot.x][jj+halfSize.y+evenOrNot.y];
	    } // jj
	  } // ii
	} //kk

      output[z][row][col] = (float)sum;
      }
}

// ---------------------------------------------------------------------------
void smoothDeformationField(Fvector3d ***dfField, Fvector3d ***dfFieldSmoothed, Ivector3d dfSize, int levelIndex, int numLevels)
{
  int i,j,k;
  float ***smoothMat;
  int gaussianRadiusXY, gaussianRadiusZ;
  float gaussianKernelXY, gaussianKernelZ;
  
  switch (levelIndex)
    {
	  case 4:  //coarsest level
	    gaussianRadiusXY = 1;
		gaussianKernelXY = 0.2;
		gaussianRadiusZ = 1;
		gaussianKernelZ = 0.2;
		break;
	  case 3:  // coarse level
	    gaussianRadiusXY = 2;
		gaussianKernelXY = 0.1 * ((float)dfSize.x/64.0);
		gaussianRadiusZ = 1;
		gaussianKernelZ = 0.075 * ((float)dfSize.z/30.0);
		break;
	  case 2:  // middle level
	    gaussianRadiusXY = 2;
		gaussianKernelXY = 0.2 * ((float)dfSize.x/128.0);
		gaussianRadiusZ = 1;
		gaussianKernelZ = 0.15 * ((float)dfSize.z/60.0);
		break;
	  case 1:
	    gaussianRadiusXY = 4;
		gaussianKernelXY = 0.8 * ((float)dfSize.x/256.0);
		gaussianRadiusZ = 2;
		gaussianKernelZ = 0.3 * ((float)dfSize.z/120.0);
		break;
	
	  default:
	    break;
	}
  
  smoothMat = Falloc3d(2*gaussianRadiusXY+1, 2*gaussianRadiusXY+1, 2*gaussianRadiusZ+1);
  gauss3D(2*gaussianRadiusXY+1, gaussianKernelXY, 2*gaussianRadiusXY+1, gaussianKernelXY, 2*gaussianRadiusZ+1, gaussianKernelZ, smoothMat);
  
  float ***dfFieldX, ***dfFieldY, ***dfFieldZ;
  float ***dfFieldXSmoothed, ***dfFieldYSmoothed, ***dfFieldZSmoothed;
  dfFieldX = Falloc3d(dfSize.x, dfSize.y, dfSize.z);
  dfFieldY = Falloc3d(dfSize.x, dfSize.y, dfSize.z);
  dfFieldZ = Falloc3d(dfSize.x, dfSize.y, dfSize.z);
  dfFieldXSmoothed = Falloc3d(dfSize.x, dfSize.y, dfSize.z);
  dfFieldYSmoothed = Falloc3d(dfSize.x, dfSize.y, dfSize.z);
  dfFieldZSmoothed = Falloc3d(dfSize.x, dfSize.y, dfSize.z);
  for (k=0;k<dfSize.z;k++)
   for (i=0;i<dfSize.x;i++)
    for (j=0;j<dfSize.y;j++)
	  {
	  dfFieldX[k][i][j] = dfField[k][i][j].x;
	  dfFieldY[k][i][j] = dfField[k][i][j].y;
	  dfFieldZ[k][i][j] = dfField[k][i][j].z;
	  }
	  
  Ivector3d smoothMatSize;
  smoothMatSize.x = 2*gaussianRadiusXY+1;
  smoothMatSize.y = 2*gaussianRadiusXY+1;
  smoothMatSize.z = 2*gaussianRadiusZ+1;
  
  Convolution3DFloat(dfFieldX, dfSize, smoothMat, smoothMatSize, dfFieldXSmoothed);
  Convolution3DFloat(dfFieldY, dfSize, smoothMat, smoothMatSize, dfFieldYSmoothed);
  Convolution3DFloat(dfFieldZ, dfSize, smoothMat, smoothMatSize, dfFieldZSmoothed);
  
  for (k=0;k<dfSize.z;k++)
   for (i=0;i<dfSize.x;i++)
    for (j=0;j<dfSize.y;j++)
	  {
		dfFieldSmoothed[k][i][j].x = dfFieldXSmoothed[k][i][j];
		dfFieldSmoothed[k][i][j].y = dfFieldYSmoothed[k][i][j];
		dfFieldSmoothed[k][i][j].z = dfFieldZSmoothed[k][i][j];
	  }
	  
  // release memory
  Ffree3d(dfFieldX, dfSize.z, dfSize.x);
  Ffree3d(dfFieldY, dfSize.z, dfSize.x);
  Ffree3d(dfFieldZ, dfSize.z, dfSize.x);
  Ffree3d(dfFieldXSmoothed, dfSize.z, dfSize.x);
  Ffree3d(dfFieldYSmoothed, dfSize.z, dfSize.x);
  Ffree3d(dfFieldZSmoothed, dfSize.z, dfSize.x);
}

// ---------------------------------------------------------------------------
void UpsampleDisplacementAtControlPoints(Fvector3d ***controlPointsPreviousLevel, Fvector3d ***controlPointsThisLevel, int numControlPointsXPreviousLevel, int numControlPointsYPreviousLevel, int numControlPointsZPreviousLevel, int numControlPointsXThisLevel, int numControlPointsYThisLevel, int numControlPointsZThisLevel, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, unsigned char ***maksForControlPointsThisLevel, int levelIndex, int numLevels)
{
  int i,j,k, ii,jj,kk, r,s,t;
  int iceil, ifloor, jceil, jfloor, kceil, kfloor;
  Ivector3d realCoordinate;
  float iOffsetFloor, iOffsetCeil, jOffsetFloor, jOffsetCeil, kOffsetFloor, kOffsetCeil;
  float offset;
  Fvector3d ***controlPointsThisLevelBackup;
  Fvector3d ***controlPointsTemp1, ***controlPointsTemp2;
  controlPointsThisLevelBackup = Fvector3dalloc3d(numControlPointsXThisLevel, numControlPointsYThisLevel, numControlPointsZThisLevel);
  controlPointsTemp1 = Fvector3dalloc3d(numControlPointsXThisLevel, numControlPointsYThisLevel, numControlPointsZThisLevel);
  controlPointsTemp2 = Fvector3dalloc3d(numControlPointsXThisLevel, numControlPointsYThisLevel, numControlPointsZThisLevel);
  
  float maxDisplacementAllowedX, maxDisplacementAllowedY;
  if (distBetweenControlPointsX<14)
    maxDisplacementAllowedX = 1.05*distBetweenControlPointsX * pow(0.95, numLevels+1-levelIndex);
  else
    maxDisplacementAllowedX = 0.9*distBetweenControlPointsX * pow(0.90, numLevels-levelIndex);
	
  if (distBetweenControlPointsY<14)
    maxDisplacementAllowedY = 1.05*distBetweenControlPointsY * pow(0.95, numLevels+1-levelIndex);
  else
    maxDisplacementAllowedY = 0.9*distBetweenControlPointsY * pow(0.90, numLevels-levelIndex);
	
  float maxDisplacementAllowedZ;
  if (distBetweenControlPointsZ<4)
    maxDisplacementAllowedZ = 0.5*distBetweenControlPointsZ * pow(0.95, numLevels+1-levelIndex);
  else 
    maxDisplacementAllowedZ = 0.65*distBetweenControlPointsZ * pow(0.90, numLevels-levelIndex);
	
  // upsample and interpolate
  for (k=0;k<MIN(numControlPointsZPreviousLevel*2, numControlPointsZThisLevel-1);k++)
   for (i=0;i<MIN(numControlPointsXPreviousLevel*2, numControlPointsXThisLevel-1);i++)
    for (j=0;j<MIN(numControlPointsYPreviousLevel*2, numControlPointsYThisLevel-1);j++)
	  {
	  realCoordinate.x = i*distBetweenControlPointsX;
	  realCoordinate.y = j*distBetweenControlPointsY;
	  realCoordinate.z = k*distBetweenControlPointsZ;
	
      //printf("(i,j) = (%d, %d)\n",i,j);
	  if (i%2==0 & j%2==0 & k%2==0) // all even numbers
	    {
		controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[k/2][i/2][j/2].x;
		controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[k/2][i/2][j/2].y;
		controlPointsThisLevelBackup[k][i][j].z = controlPointsPreviousLevel[k/2][i/2][j/2].z;
		}
	  else if (i%2==0 & j%2!=0 & k%2==0) // i even, j odd, k even
	    {
		jfloor = (int)floor((float)j/2.0);
		jceil = (int)ceil((float)j/2.0);
		if (jceil>(numControlPointsYPreviousLevel-1))
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[k/2][i/2][jfloor].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[k/2][i/2][jfloor].y;
		  controlPointsThisLevelBackup[k][i][j].z = controlPointsPreviousLevel[k/2][i/2][jfloor].z;
		  }
		else
		  {
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (jOffsetFloor*controlPointsPreviousLevel[k/2][i/2][jfloor].x + jOffsetCeil*controlPointsPreviousLevel[k/2][i/2][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (jOffsetFloor*controlPointsPreviousLevel[k/2][i/2][jfloor].y + jOffsetCeil*controlPointsPreviousLevel[k/2][i/2][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (jOffsetFloor*controlPointsPreviousLevel[k/2][i/2][jfloor].z + jOffsetCeil*controlPointsPreviousLevel[k/2][i/2][jceil].z);
		  }
		}
	  else if (i%2!=0 & j%2==0 & k%2==0) // i odd, j even, k even
	    {
		ifloor = (int)floor((float)i/2.0);
		iceil = (int)ceil((float)i/2.0);
		if (iceil>(numControlPointsXPreviousLevel-1))
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[k/2][ifloor][j/2].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[k/2][ifloor][j/2].y;
		  controlPointsThisLevelBackup[k][i][j].z = controlPointsPreviousLevel[k/2][ifloor][j/2].z;
		  }
		else
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][j/2].x + iOffsetCeil*controlPointsPreviousLevel[k/2][iceil][j/2].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][j/2].y + iOffsetCeil*controlPointsPreviousLevel[k/2][iceil][j/2].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][j/2].z + iOffsetCeil*controlPointsPreviousLevel[k/2][iceil][j/2].z);
		  }
		} //else if (i%2!=0 & j%2==0) // i odd, j even
	  else if (i%2!=0 & j%2!=0 & k%2==0) // both i and j are odd, k even
	    {
		ifloor = (int)floor((float)i/2.0);
		iceil = (int)ceil((float)i/2.0);
		jfloor = (int)floor((float)j/2.0);
		jceil = (int)ceil((float)j/2.0);
		if (iceil>(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1))
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[k/2][ifloor][jfloor].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[k/2][ifloor][jfloor].y;
		  }
		else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1))
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].x + iOffsetCeil*controlPointsPreviousLevel[k/2][iceil][jfloor].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].y + iOffsetCeil*controlPointsPreviousLevel[k/2][iceil][jfloor].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].z + iOffsetCeil*controlPointsPreviousLevel[k/2][iceil][jfloor].z);
		  }
		else if (iceil>(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1))
		  {
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (jOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].x + jOffsetCeil*controlPointsPreviousLevel[k/2][ifloor][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (jOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].y + jOffsetCeil*controlPointsPreviousLevel[k/2][ifloor][jceil].y);
		  }
		else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1))
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].x + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[k/2][ifloor][jceil].x + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[k/2][iceil][jfloor].x + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[k/2][iceil][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].y + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[k/2][ifloor][jceil].y + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[k/2][iceil][jfloor].y + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[k/2][iceil][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[k/2][ifloor][jfloor].z + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[k/2][ifloor][jceil].z + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[k/2][iceil][jfloor].z + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[k/2][iceil][jceil].z);
		  }
		} // else if (i%2!=0 & j%2!=0) // i and j both odd, k even
		
	  if (i%2==0 & j%2==0 & k%2!=0) // i and j both even numbers, k odd
	    {
		kceil = (int)ceil((float)k/2.0);
		kfloor = (int)floor((float)k/2.0);
		
		if (kceil>(numControlPointsZPreviousLevel-1))
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[kfloor][i/2][j/2].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[kfloor][i/2][j/2].y;
		  controlPointsThisLevelBackup[k][i][j].z = controlPointsPreviousLevel[kfloor][i/2][j/2].z;
		  }
		else
		  {
		  kOffsetFloor = ceil((float)k/2.0) - (float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][j/2].x + kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][j/2].x;
		  controlPointsThisLevelBackup[k][i][j].y = kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][j/2].y + kOffsetFloor*controlPointsPreviousLevel[kceil][i/2][j/2].y;
		  controlPointsThisLevelBackup[k][i][j].z = kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][j/2].z + kOffsetFloor*controlPointsPreviousLevel[kceil][i/2][j/2].z;
		  }
		}
	  else if (i%2==0 & j%2!=0 & k%2!=0) // i even, j odd, k odd
	    {
		jfloor = (int)floor((float)j/2.0);
		jceil = (int)ceil((float)j/2.0);
		kfloor = (int)floor((float)k/2.0);
		kceil = (int)ceil((float)k/2.0);
		if ( jceil>(numControlPointsYPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1) )
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[kfloor][i/2][jfloor].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[kfloor][i/2][jfloor].y;
		  }
		else if ( jceil<=(numControlPointsYPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1) )
		  {
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (jOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].x + jOffsetCeil*controlPointsPreviousLevel[kfloor][i/2][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (jOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].y + jOffsetCeil*controlPointsPreviousLevel[kfloor][i/2][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (jOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].z + jOffsetCeil*controlPointsPreviousLevel[kfloor][i/2][jceil].z);
		  }
		else if ( jceil>(numControlPointsYPreviousLevel-1) & kceil<=(numControlPointsZPreviousLevel-1) )
		  {
		  kOffsetFloor = ceil((float)k/2.0) - (float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].x + kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jfloor].x);
		  controlPointsThisLevelBackup[k][i][j].y = (kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].y + kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jfloor].y);
		  controlPointsThisLevelBackup[k][i][j].z = (kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].z + kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jfloor].z);
		  }
		else if  ( jceil<=(numControlPointsYPreviousLevel-1) & kceil<=(numControlPointsZPreviousLevel-1) )
		  {
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  kOffsetFloor = ceil((float)k/2.0) - (float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].x + jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jfloor].x + jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jceil].x + jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].y + jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jfloor].y + jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jceil].y + jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jfloor].z + jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jfloor].z + jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][i/2][jceil].z + jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][i/2][jceil].z);
		  }
		}
	  else if (i%2!=0 & j%2==0 & k%2!=0) // i odd, j even, k odd
	    {
		ifloor = (int)floor((float)i/2.0);
		iceil = (int)ceil((float)i/2.0);
		kfloor = (int)floor((float)k/2.0);
		kceil = (int)ceil((float)k/2.0);
		if ( iceil>(numControlPointsXPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1) )
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[kfloor][ifloor][j/2].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[kfloor][ifloor][j/2].y;
		  }
		else if ( iceil<=(numControlPointsXPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1) )
		  {
		  iOffsetFloor = ceil((float)i/2.0)-(float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].x + iOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][j/2].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].y + iOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][j/2].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].z + iOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][j/2].z);
		  }
		else if ( iceil>(numControlPointsXPreviousLevel-1) & kceil<=(numControlPointsZPreviousLevel-1) )
		  {
		  kOffsetFloor = ceil((float)k/2.0) - (float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].x + kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][j/2].x);
		  controlPointsThisLevelBackup[k][i][j].y = (kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].y + kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][j/2].y);
		  controlPointsThisLevelBackup[k][i][j].z = (kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].z + kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][j/2].z);
		  }
		else if  ( iceil<=(numControlPointsXPreviousLevel-1) & kceil<=(numControlPointsZPreviousLevel-1) )
		  {
		  iOffsetFloor = ceil((float)i/2.0)-(float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  kOffsetFloor = ceil((float)k/2.0) - (float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].x + iOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][j/2].x + iOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][j/2].x + iOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][j/2].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].y + iOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][j/2].y + iOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][j/2].y + iOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][j/2].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][j/2].z + iOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][j/2].z + iOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][j/2].z + iOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][j/2].z);
		  }
		}
	  else if (i%2!=0 & j%2!=0 & k%2!=0) // all odd
	    {
		ifloor = (int)floor((float)i/2.0);
		iceil = (int)ceil((float)i/2.0);
		jfloor = (int)floor((float)j/2.0);
		jceil = (int)ceil((float)j/2.0);
		kfloor = (int)floor((float)k/2.0);
		kceil = (int)ceil((float)k/2.0);
		if (iceil>(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1))
		  {
		  controlPointsThisLevelBackup[k][i][j].x = controlPointsPreviousLevel[kfloor][ifloor][jfloor].x;
		  controlPointsThisLevelBackup[k][i][j].y = controlPointsPreviousLevel[kfloor][ifloor][jfloor].y;
		  controlPointsThisLevelBackup[k][i][j].z = controlPointsPreviousLevel[kfloor][ifloor][jfloor].z;
		  }
		else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1))
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + iOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][jfloor].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + iOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][jfloor].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + iOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][jfloor].z);
		  }
		else if (iceil>(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1) & kceil>(numControlPointsZPreviousLevel-1))
		  {
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (jOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + jOffsetCeil*controlPointsPreviousLevel[kfloor][ifloor][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (jOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + jOffsetCeil*controlPointsPreviousLevel[kfloor][ifloor][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (jOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + jOffsetCeil*controlPointsPreviousLevel[kfloor][ifloor][jceil].z);
		  }
		else if (iceil>(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1) & kceil<=(numControlPointsZPreviousLevel-1))
		  {
		  kOffsetFloor = ceil((float)k/2.0)-(float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].x);
		  controlPointsThisLevelBackup[k][i][j].y = (kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].y);
		  controlPointsThisLevelBackup[k][i][j].z = (kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].z);
		  }
		else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1) & k>(numControlPointsZPreviousLevel-1))
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  jOffsetFloor = ceil((float)j/2.0)-(float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[kfloor][ifloor][jceil].x + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].x + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[kfloor][ifloor][jceil].y + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].y + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*jOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + iOffsetFloor*jOffsetCeil*controlPointsPreviousLevel[kfloor][ifloor][jceil].z + iOffsetCeil*jOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].z + iOffsetCeil*jOffsetCeil*controlPointsPreviousLevel[kfloor][iceil][jceil].z);
		  }
		else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil>(numControlPointsYPreviousLevel-1) & k<=(numControlPointsZPreviousLevel-1))
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  kOffsetFloor = ceil((float)k/2.0)-(float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + iOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].x + iOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].x + iOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jfloor].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + iOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].y + iOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].y + iOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jfloor].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + iOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].z + iOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].z + iOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jfloor].z);
		  }
		else if (iceil>(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1) & k<=(numControlPointsZPreviousLevel-1))
		  {
		  jOffsetFloor = ceil((float)j/2.0) - (float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  kOffsetFloor = ceil((float)k/2.0)-(float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].x + jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jceil].x + jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].y + jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jceil].y + jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].z + jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jceil].z + jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jceil].z);
		  }
		else if (iceil<=(numControlPointsXPreviousLevel-1) & jceil<=(numControlPointsYPreviousLevel-1) & k<=(numControlPointsZPreviousLevel-1))
		  {
		  iOffsetFloor = ceil((float)i/2.0) - (float)i/2.0;
		  iOffsetCeil = (float)i/2.0 - floor((float)i/2.0);
		  jOffsetFloor = ceil((float)j/2.0) - (float)j/2.0;
		  jOffsetCeil = (float)j/2.0 - floor((float)j/2.0);
		  kOffsetFloor = ceil((float)k/2.0)-(float)k/2.0;
		  kOffsetCeil = (float)k/2.0 - floor((float)k/2.0);
		  controlPointsThisLevelBackup[k][i][j].x = (iOffsetFloor*jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].x + iOffsetFloor*jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].x + iOffsetFloor*jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jceil].x + iOffsetFloor*jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jceil].x + iOffsetCeil*jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].x + iOffsetCeil*jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jfloor].x + iOffsetCeil*jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jceil].x + iOffsetCeil*jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jceil].x);
		  controlPointsThisLevelBackup[k][i][j].y = (iOffsetFloor*jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].y + iOffsetFloor*jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].y + iOffsetFloor*jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jceil].y + iOffsetFloor*jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jceil].y + iOffsetCeil*jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].y + iOffsetCeil*jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jfloor].y + iOffsetCeil*jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jceil].y + iOffsetCeil*jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jceil].y);
		  controlPointsThisLevelBackup[k][i][j].z = (iOffsetFloor*jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jfloor].z + iOffsetFloor*jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jfloor].z + iOffsetFloor*jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][ifloor][jceil].z + iOffsetFloor*jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][ifloor][jceil].z + iOffsetCeil*jOffsetFloor*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jfloor].z + iOffsetCeil*jOffsetFloor*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jfloor].z + iOffsetCeil*jOffsetCeil*kOffsetFloor*controlPointsPreviousLevel[kfloor][iceil][jceil].z + iOffsetCeil*jOffsetCeil*kOffsetCeil*controlPointsPreviousLevel[kceil][iceil][jceil].z);
		  }
		} // else if (i%2!=0 & j%2!=0) // both odd
		
		controlPointsThisLevelBackup[k][i][j].x *= 2;
		controlPointsThisLevelBackup[k][i][j].y *= 2;
		controlPointsThisLevelBackup[k][i][j].z *= 2;
	  }// for each control point this level

  // step 1: copy controlPointsBackup to controlPoints
  for (k=0;k<numControlPointsZThisLevel;k++)
   for (i=0;i<numControlPointsXThisLevel;i++)
	for (j=0;j<numControlPointsYThisLevel;j++)
	  {
	  controlPointsTemp1[k][i][j].x = controlPointsThisLevelBackup[k][i][j].x;
	  controlPointsTemp1[k][i][j].y = controlPointsThisLevelBackup[k][i][j].y;
	  controlPointsTemp1[k][i][j].z = controlPointsThisLevelBackup[k][i][j].z;
	  controlPointsTemp2[k][i][j].x = controlPointsThisLevelBackup[k][i][j].x;
	  controlPointsTemp2[k][i][j].y = controlPointsThisLevelBackup[k][i][j].y;
	  controlPointsTemp2[k][i][j].z = controlPointsThisLevelBackup[k][i][j].z;
	  }

	// step 3: if displacement at a control point exceeds maximum allowed, distribute it along the gradient direction.
	for (k=0;k<numControlPointsZThisLevel;k++)
	 for (i=0;i<numControlPointsXThisLevel;i++)
      for (j=0;j<numControlPointsYThisLevel;j++)
	    {
		// x directions, -
		if (fabs(controlPointsTemp1[k][i][j].x)>maxDisplacementAllowedX)
		  {
		  offset = controlPointsTemp1[k][i][j].x - sign(controlPointsTemp1[k][i][j].x)*maxDisplacementAllowedX;
		  controlPointsTemp1[k][i][j].x = sign(controlPointsTemp1[k][i][j].x)*maxDisplacementAllowedX;
		  
		  ii = i-1;
		  while (ii>=0 && ii<numControlPointsXThisLevel)
		     {
			 controlPointsTemp1[k][ii][j].x += offset*pow(0.85, fabs((float)(ii-i)));
			 if (fabs(controlPointsTemp1[k][ii][j].x)>maxDisplacementAllowedX)
			   {
			   offset = controlPointsTemp1[k][ii][j].x - sign(controlPointsTemp1[k][i][j].x)*maxDisplacementAllowedX;
			   controlPointsTemp1[k][ii][j].x = sign(controlPointsTemp1[k][i][j].x)*maxDisplacementAllowedX;
			   ii += -1;
			   } // if(ii,j) 
			 else
			   break;
			 } // while
		  }// if (i,j)
		  
		// y direction, -
		if (fabs(controlPointsTemp1[k][i][j].y)>maxDisplacementAllowedX)
		  {
		  offset = controlPointsTemp1[k][i][j].y - sign(controlPointsTemp1[k][i][j].y)* maxDisplacementAllowedY;
		  controlPointsTemp1[k][i][j].y = sign(controlPointsTemp1[k][i][j].y)*maxDisplacementAllowedY;
		  
		  jj = j-1;
		  while (jj>=0 && jj<numControlPointsYThisLevel)
		     {
			 controlPointsTemp1[k][i][jj].y += offset*pow(0.85, fabs((float)(jj-j)));
			 if (fabs(controlPointsTemp1[k][i][jj].y)>maxDisplacementAllowedY)
			   {
			   offset = controlPointsTemp1[k][i][jj].y - sign(controlPointsTemp1[k][i][j].y)*maxDisplacementAllowedY;
			   controlPointsTemp1[k][i][jj].y = sign(controlPointsTemp1[k][i][j].y)*maxDisplacementAllowedY;
			   jj += -1;
			   } // if(i,jj) 
			 else
			   break;
			 } // while
		  }// if (i,j)
		
		// z direction, -
		if (fabs(controlPointsTemp1[k][i][j].z)>maxDisplacementAllowedZ)
		  {
		  offset = controlPointsTemp1[k][i][j].z - sign(controlPointsTemp1[k][i][j].z)* maxDisplacementAllowedZ;
		  controlPointsTemp1[k][i][j].z = sign(controlPointsTemp1[k][i][j].z)*maxDisplacementAllowedZ;
		  
		  kk = k-1;
		  while (kk>=0 && kk<numControlPointsZThisLevel)
		     {
			 controlPointsTemp1[kk][i][j].z += offset*pow(0.85, fabs((float)(kk-k)));
			 if (fabs(controlPointsTemp1[kk][i][j].z)>maxDisplacementAllowedZ)
			   {
			   offset = controlPointsTemp1[kk][i][j].z - sign(controlPointsTemp1[k][i][j].z)*maxDisplacementAllowedZ;
			   controlPointsTemp1[kk][i][j].z = sign(controlPointsTemp1[k][i][j].z)*maxDisplacementAllowedZ;
			   kk += -1;
			   } // if(i,jj) 
			 else
			   break;
			 } // while
		  }// if (i,j)
		} // for for
	
	for (k=(numControlPointsZThisLevel-1);k>=0;k--)
	 for (i=(numControlPointsXThisLevel-1);i>=0;i--)
      for (j=(numControlPointsYThisLevel-1);j>=0;j--)
	    {
		// x directions, +
		if (fabs(controlPointsTemp2[k][i][j].x)>maxDisplacementAllowedX)
		  {
		  offset = controlPointsTemp2[k][i][j].x - sign(controlPointsTemp2[k][i][j].x)*maxDisplacementAllowedX;
		  controlPointsTemp2[k][i][j].x = sign(controlPointsTemp2[k][i][j].x)*maxDisplacementAllowedX;
		  
		  ii = i+1;
		  while (ii>=0 && ii<numControlPointsXThisLevel)
		     {
			 controlPointsTemp2[k][ii][j].x += offset*pow(0.85, fabs((float)(ii-i)));
			 if (fabs(controlPointsTemp2[k][ii][j].x)>maxDisplacementAllowedX)
			   {
			   offset = controlPointsTemp2[k][ii][j].x - sign(controlPointsTemp2[k][i][j].x)*maxDisplacementAllowedX;
			   controlPointsTemp2[k][ii][j].x = sign(controlPointsTemp2[k][i][j].x)*maxDisplacementAllowedX;
			   ii += 1;
			   } // if(ii,j) 
			 else
			   break;
			 } // while
		  }// if (i,j)
		  
		// y direction, +
		if (fabs(controlPointsTemp2[k][i][j].y)>maxDisplacementAllowedY)
		  {
		  offset = controlPointsTemp2[k][i][j].y - sign(controlPointsTemp2[k][i][j].y)*maxDisplacementAllowedY;
		  controlPointsTemp2[k][i][j].y = sign(controlPointsTemp2[k][i][j].y)*maxDisplacementAllowedY;
		  
		  jj = j+1;
		  while (jj>=0 && jj<numControlPointsYThisLevel)
		     {
			 controlPointsTemp2[k][i][jj].y += offset*pow(0.85, fabs((float)(jj-j)));
			 if (fabs(controlPointsTemp2[k][i][jj].y)>maxDisplacementAllowedY)
			   {
			   offset = controlPointsTemp2[k][i][jj].y - sign(controlPointsTemp2[k][i][j].y)*maxDisplacementAllowedY;
			   controlPointsTemp2[k][i][jj].y = sign(controlPointsTemp2[k][i][j].y)*maxDisplacementAllowedY;
			   jj += 1;
			   } // if(i,jj) 
			 else
			   break;
			 } // while
		  }// if (i,j)
		  
		// z direction, +
		if (fabs(controlPointsTemp2[k][i][j].z)>maxDisplacementAllowedZ)
		  {
		  offset = controlPointsTemp2[k][i][j].z - sign(controlPointsTemp2[k][i][j].z)*maxDisplacementAllowedZ;
		  controlPointsTemp2[k][i][j].z = sign(controlPointsTemp2[k][i][j].z)*maxDisplacementAllowedZ;
		  
		  kk = k+1;
		  while (kk>=0 && kk<numControlPointsZThisLevel)
		     {
			 controlPointsTemp2[kk][i][j].z += offset*pow(0.85, fabs((float)(kk-k)));
			 if (fabs(controlPointsTemp2[kk][i][j].z)>maxDisplacementAllowedZ)
			   {
			   offset = controlPointsTemp2[kk][i][j].z - sign(controlPointsTemp2[k][i][j].z)*maxDisplacementAllowedZ;
			   controlPointsTemp2[kk][i][j].z = sign(controlPointsTemp2[k][i][j].z)*maxDisplacementAllowedZ;
			   kk += 1;
			   } // if(i,jj) 
			 else
			   break;
			 } // while
		  }// if (i,j)
		} // for for
	
	for (k=0;k<numControlPointsZThisLevel;k++)
	 for (i=0;i<numControlPointsXThisLevel;i++)
      for (j=0;j<numControlPointsYThisLevel;j++)
	    {
		// x
		if (fabs(controlPointsTemp1[k][i][j].x)>fabs(controlPointsTemp2[k][i][j].x))
		  controlPointsThisLevel[k][i][j].x = 0.7 * controlPointsTemp1[k][i][j].x + 0.3* controlPointsTemp2[k][i][j].x;
		else
		  controlPointsThisLevel[k][i][j].x = 0.3 * controlPointsTemp1[k][i][j].x + 0.7* controlPointsTemp2[k][i][j].x;
		  
		//y
		if (fabs(controlPointsTemp1[k][i][j].y)>fabs(controlPointsTemp2[k][i][j].y))
		  controlPointsThisLevel[k][i][j].y = 0.7 * controlPointsTemp1[k][i][j].y + 0.3* controlPointsTemp2[k][i][j].y;
		else
		  controlPointsThisLevel[k][i][j].y = 0.3 * controlPointsTemp1[k][i][j].y + 0.7* controlPointsTemp2[k][i][j].y;
		  
		//z
		if (fabs(controlPointsTemp1[k][i][j].z)>fabs(controlPointsTemp2[k][i][j].z))
		  controlPointsThisLevel[k][i][j].z = 0.7 * controlPointsTemp1[k][i][j].z + 0.3* controlPointsTemp2[k][i][j].z;
		else
		  controlPointsThisLevel[k][i][j].z = 0.3 * controlPointsTemp1[k][i][j].z + 0.7* controlPointsTemp2[k][i][j].z;  
		}
	
	Fvector3dfree3d(controlPointsThisLevelBackup, numControlPointsZThisLevel, numControlPointsXThisLevel);
	Fvector3dfree3d(controlPointsTemp1, numControlPointsZThisLevel, numControlPointsXThisLevel);
	Fvector3dfree3d(controlPointsTemp2, numControlPointsZThisLevel, numControlPointsXThisLevel);
}

// ---------------------------------------------------------------------------
void UpsampleDeformationField(Fvector3d ***dfFieldB2APreviousLevel, Fvector3d ***dfFieldB2AThisLevel, Ivector3d imageSizePreviousLevel, Ivector3d imageSizeThisLevel, Fvector3d resolutionRatio, int smoothOrNot, int levelIndex, int numLevels, float foregroundRatio)
{
    int i,j,k;
	int xFloor,yFloor,zFloor;
	int xCeil,yCeil,zCeil;
	float upsampleFactor; 
	
	
	if (foregroundRatio>0.35)		
		upsampleFactor=1.9;
	else
		{
		if ( (MAX(imageSizeThisLevel.x*resolutionRatio.x, imageSizeThisLevel.y*resolutionRatio.y)>150*pow(0.5, levelIndex-1)) && (imageSizeThisLevel.z*resolutionRatio.z>150*pow(0.5,levelIndex-1)) )
			upsampleFactor=1.6;
		else
			upsampleFactor=1.2;
		}
	printf("upsample factor=%f\n", upsampleFactor);
	
		
	#pragma omp paralle for private(i,j,k,xFloor,yFloor,zFloor,xCeil,yCeil,zCeil) num_threads(50)
	for (k=0;k<imageSizeThisLevel.z;k++)
	  for (i=0;i<imageSizeThisLevel.x;i++)
	    for (j=0;j<imageSizeThisLevel.y;j++)
		  {
		  xFloor = (int)MIN((float)(imageSizePreviousLevel.x-1),floor((float)i/2.0));
		  yFloor = (int)MIN((float)(imageSizePreviousLevel.y-1),floor((float)j/2.0));
		  zFloor = (int)MIN((float)(imageSizePreviousLevel.z-1),floor((float)k/2.0));
		  xCeil = (int)MIN((float)(imageSizePreviousLevel.x-1),floor(0.5+(float)i/2.0));
		  yCeil = (int)MIN((float)(imageSizePreviousLevel.y-1),floor(0.5+(float)j/2.0));
		  zCeil = (int)MIN((float)(imageSizePreviousLevel.z-1),floor(0.5+(float)k/2.0));
		  
		  dfFieldB2AThisLevel[k][i][j].x = 0.125* (
		    dfFieldB2APreviousLevel[zFloor][xFloor][yFloor].x +
			dfFieldB2APreviousLevel[zFloor][xFloor][yCeil].x +
			dfFieldB2APreviousLevel[zFloor][xCeil][yFloor].x +
			dfFieldB2APreviousLevel[zFloor][xCeil][yCeil].x +
			dfFieldB2APreviousLevel[zCeil][xFloor][yFloor].x +
			dfFieldB2APreviousLevel[zCeil][xFloor][yCeil].x +
			dfFieldB2APreviousLevel[zCeil][xCeil][yFloor].x +
			dfFieldB2APreviousLevel[zCeil][xCeil][yCeil].x
			) * upsampleFactor;
		
		  dfFieldB2AThisLevel[k][i][j].y = 0.125* (
		    dfFieldB2APreviousLevel[zFloor][xFloor][yFloor].y +
			dfFieldB2APreviousLevel[zFloor][xFloor][yCeil].y +
			dfFieldB2APreviousLevel[zFloor][xCeil][yFloor].y +
			dfFieldB2APreviousLevel[zFloor][xCeil][yCeil].y +
			dfFieldB2APreviousLevel[zCeil][xFloor][yFloor].y +
			dfFieldB2APreviousLevel[zCeil][xFloor][yCeil].y +
			dfFieldB2APreviousLevel[zCeil][xCeil][yFloor].y +
			dfFieldB2APreviousLevel[zCeil][xCeil][yCeil].y
			) * upsampleFactor;
			
		  dfFieldB2AThisLevel[k][i][j].z = 0.125* (
		    dfFieldB2APreviousLevel[zFloor][xFloor][yFloor].z +
			dfFieldB2APreviousLevel[zFloor][xFloor][yCeil].z +
			dfFieldB2APreviousLevel[zFloor][xCeil][yFloor].z +
			dfFieldB2APreviousLevel[zFloor][xCeil][yCeil].z +
			dfFieldB2APreviousLevel[zCeil][xFloor][yFloor].z +
			dfFieldB2APreviousLevel[zCeil][xFloor][yCeil].z +
			dfFieldB2APreviousLevel[zCeil][xCeil][yFloor].z +
			dfFieldB2APreviousLevel[zCeil][xCeil][yCeil].z
			) * upsampleFactor;
		  
		  
		  }
	
	if (smoothOrNot==YYES)
		{
		Fvector3d ***defTemp;
		defTemp = Fvector3dalloc3d(imageSizeThisLevel.x, imageSizeThisLevel.y, imageSizeThisLevel.z);
		
		smoothDeformationField(dfFieldB2AThisLevel, defTemp, imageSizeThisLevel, levelIndex, numLevels);
		
		for (k=0;k<imageSizeThisLevel.z;k++)
		  for (i=0;i<imageSizeThisLevel.x;i++)
			for (j=0;j<imageSizeThisLevel.y;j++)
			  {
			  dfFieldB2AThisLevel[k][i][j].x = defTemp[k][i][j].x;
			  dfFieldB2AThisLevel[k][i][j].y = defTemp[k][i][j].y;
			  dfFieldB2AThisLevel[k][i][j].z = defTemp[k][i][j].z;
			  }
			  
		Fvector3dfree3d(defTemp, imageSizeThisLevel.z, imageSizeThisLevel.x);
		}
}

// ---------------------------------------------------------------------------
void UpsampleConfidenceMap(float ***MSMapPreviousLevel, float ***MSMapThisLevel, Ivector3d imageSizePreviousLevel, Ivector3d imageSizeThisLevel)
{
    int i,j,k;
	int xFloor,yFloor,zFloor;
	int xCeil,yCeil,zCeil;
		
	#pragma omp paralle for private(i,j,k,xFloor,yFloor,zFloor,xCeil,yCeil,zCeil) num_threads(50)
	for (k=0;k<imageSizeThisLevel.z;k++)
	  for (i=0;i<imageSizeThisLevel.x;i++)
	    for (j=0;j<imageSizeThisLevel.y;j++)
		  {
		  xFloor = (int)MIN((float)(imageSizePreviousLevel.x-1),floor((float)i/2.0));
		  yFloor = (int)MIN((float)(imageSizePreviousLevel.y-1),floor((float)j/2.0));
		  zFloor = (int)MIN((float)(imageSizePreviousLevel.z-1),floor((float)k/2.0));
		  xCeil = (int)MIN((float)(imageSizePreviousLevel.x-1),floor(0.5+(float)i/2.0));
		  yCeil = (int)MIN((float)(imageSizePreviousLevel.y-1),floor(0.5+(float)j/2.0));
		  zCeil = (int)MIN((float)(imageSizePreviousLevel.z-1),floor(0.5+(float)k/2.0));

		  MSMapThisLevel[k][i][j] = 0.125* (
		    MSMapPreviousLevel[zFloor][xFloor][yFloor] +
			MSMapPreviousLevel[zFloor][xFloor][yCeil] +
			MSMapPreviousLevel[zFloor][xCeil][yFloor] +
			MSMapPreviousLevel[zFloor][xCeil][yCeil] +
			MSMapPreviousLevel[zCeil][xFloor][yFloor] +
			MSMapPreviousLevel[zCeil][xFloor][yCeil] +
			MSMapPreviousLevel[zCeil][xCeil][yFloor] +
			MSMapPreviousLevel[zCeil][xCeil][yCeil]
			);
		  }
}

// ---------------------------------------------------------------------------
int GenerateMask(unsigned char ***imageA, unsigned char ***imageB, Ivector3d imageSize, unsigned char ***maskForImages, unsigned char ***maskForControlPoints, int levelIndex, int numLevels, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int foregroundThre)
{
   int i,j,k, x,y,z;
   int ii,jj,kk;
   int counter=0;
   int dilateRadiusX = (int)(2.3*distBetweenControlPointsX*pow(1.3, numLevels-levelIndex)+0.5);
   int dilateRadiusY = (int)(2.3*distBetweenControlPointsY*pow(1.3, numLevels-levelIndex)+0.5);
   int dilateRadiusZ = (int)(2.3*distBetweenControlPointsZ*pow(1.3, numLevels-levelIndex)+0.5);
   float dilationFactorX, dilationFactorY, dilationFactorZ;
   int s,t,r;
   int sNbrMin, tNbrMin, rNbrMin;
   int sNbrMax, tNbrMax, rNbrMax;
   int sFillMin, tFillMin, rFillMin;
   int sFillMax, tFillMax, rFillMax;
   int chk1, chk2, chkX, chkY, chkZ;
   int deltaX1, deltaY1, deltaZ1;
   int deltaX2, deltaY2, deltaZ2;
   
   printf("dilation radius = (%d,%d,%d)\n", dilateRadiusX, dilateRadiusY, dilateRadiusZ);
   
   
   
   #pragma omp parallel for private(i,j,k) reduction(+:counter) num_threads(100)
   for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
     for (j=0;j<imageSize.y;j++)
	   {  
		maskForImages[k][i][j]=0;
		maskForControlPoints[k][i][j]=0;
		
		if (imageB[k][i][j]>foregroundThre || imageA[k][i][j]>foregroundThre)
		  {
		  maskForImages[k][i][j]=255;
		  maskForControlPoints[k][i][j]=255;
		  counter++;
		  }		
	   } // for each point

  
  #pragma omp parallel for shared(maskForImages, maskForControlPoints, imageSize, dilateRadiusX, dilateRadiusY, dilateRadiusZ) private(i,j,k,chk1,chk2,rNbrMin,rNbrMax,sNbrMin,sNbrMax,tNbrMin,tNbrMax,r,s,t,rFillMin,rFillMax,sFillMin,sFillMax,tFillMin,tFillMax) num_threads(100)
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
     for (j=0;j<imageSize.y;j++)
	   {
	   // detect if [x,y,z] is an edge point
	   
	   if (maskForImages[k][i][j]==0)
	   {
	   
	   chk1=0;
	   chk2=0;
	   
	   rNbrMin=MAX(0,k-1);
	   rNbrMax=MIN(imageSize.z-1,k+1);
	   sNbrMin=MAX(0,i-1);
	   sNbrMax=MIN(imageSize.x-1,i+1);
	   tNbrMin=MAX(0,j-1);
	   tNbrMax=MIN(imageSize.y-1,j+1);
	   for (r=rNbrMin;r<=rNbrMax;r++)
		for (s=sNbrMin;s<=sNbrMax;s++)
	     for (t=tNbrMin;t<=tNbrMax;t++)
		   {
		   if (maskForImages[r][s][t]>0)   chk1=1;
		   if (maskForImages[r][s][t]==0)  chk2=1;
		   }
	   
	   if ( (chk1==1)&&(chk2==1) )
	     {
		 rFillMin = MAX(0, k-dilateRadiusZ);
		 rFillMax = MIN(imageSize.z-1, k+dilateRadiusZ);
		 sFillMin = MAX(0, i-dilateRadiusX);
		 sFillMax = MIN(imageSize.x-1, i+dilateRadiusX);
		 tFillMin = MAX(0, j-dilateRadiusY);
		 tFillMax = MIN(imageSize.y-1, j+dilateRadiusY);
		 for (r=rFillMin;r<=rFillMax;r++)
		  for (s=sFillMin;s<=sFillMax;s++)
		    for (t=tFillMin;t<=tFillMax;t++)
		      maskForControlPoints[r][s][t]=255;			   
		 }
	   }
	   }
	   
	   
   printf("done!\n");
   
   printf("total number of foreground voxels = %d\n", counter);
   return counter;
}

// ---------------------------------------------------------------------------
int DilateImageMaskIntoControlPointMask(unsigned char ***maskForImages, unsigned char ***maskForControlPoints, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int levelIndex, int numLevels)
{
  int ii,jj,kk;
  int count=0;
  int chk1, chk2;
  int r,s,t;
  int dilateRadiusX = (int)(2.3*distBetweenControlPointsX*pow(1.3, numLevels-levelIndex)+0.5)+1;
  int dilateRadiusY = (int)(2.3*distBetweenControlPointsY*pow(1.3, numLevels-levelIndex)+0.5)+1;
  int dilateRadiusZ = (int)(2.3*distBetweenControlPointsZ*pow(1.3, numLevels-levelIndex)+0.5)+1;
  int smin, smax, rmin, rmax, tmin, tmax;  
  
  int maxCoorInRows[imageSize.z][imageSize.x];  // to record the maximum coordinate of the point in the foreground
  int minCoorInRows[imageSize.z][imageSize.x];
  int maxCoorInColumns[imageSize.z][imageSize.y];
  int minCoorInColumns[imageSize.z][imageSize.y];
  int maxCoorVertical[imageSize.x][imageSize.y];
  int minCoorVertical[imageSize.x][imageSize.y];  
  
  unsigned char ***maskForImagesWithHolesFilled;
  maskForImagesWithHolesFilled = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);
  
  // --------------------------
  //   fill in all holes , begin (Yangming added on Sept 14, 2010)
  // --------------------------
  
  // initialize
  for (kk=0;kk<imageSize.z;kk++)
    {
    for (ii=0;ii<imageSize.x;ii++)
	  {
	  maxCoorInRows[kk][ii]=0;
	  minCoorInRows[kk][ii]=imageSize.y-1;
	  }
	for (jj=0;jj<imageSize.y;jj++)
      {
	  maxCoorInColumns[kk][jj]=0;
	  minCoorInColumns[kk][jj]=imageSize.x-1;
	  }
	}
  for (ii=0;ii<imageSize.x;ii++)
    for (jj=0;jj<imageSize.y;jj++)
      {
	  maxCoorVertical[ii][jj]=0;
	  minCoorVertical[ii][jj]=imageSize.z-1;
	  }	
  
  
  for (kk=0;kk<imageSize.z;kk++)
    for (ii=0;ii<imageSize.x;ii++)
	  for (jj=0;jj<imageSize.y;jj++)
	    { //for...for...for...
		maskForImagesWithHolesFilled[kk][ii][jj]=maskForImages[kk][ii][jj];
		if (maskForImages[kk][ii][jj]>0)
		   {
		   count++;
		   if (jj>maxCoorInRows[kk][ii])      maxCoorInRows[kk][ii]=jj;
		   if (jj<minCoorInRows[kk][ii])      minCoorInRows[kk][ii]=jj;
		   if (ii>maxCoorInColumns[kk][jj])   maxCoorInColumns[kk][jj]=ii;
		   if (ii<minCoorInColumns[kk][jj])   minCoorInColumns[kk][jj]=ii;
		   if (kk>maxCoorVertical[ii][jj])    maxCoorVertical[ii][jj]=kk;
		   if (kk<minCoorVertical[ii][jj])    minCoorVertical[ii][jj]=kk;
		   }
		}
		
  for (kk=0;kk<imageSize.z;kk++)
    {
    for (ii=0;ii<imageSize.x;ii++)
	  for (jj=minCoorInRows[kk][ii];jj<maxCoorInRows[kk][ii];jj++)
	    maskForImagesWithHolesFilled[kk][ii][jj]=255;
	
	for (jj=0;jj<imageSize.y;jj++)
	  for (ii=minCoorInColumns[kk][jj];ii<maxCoorInColumns[kk][jj];ii++)
	    maskForImagesWithHolesFilled[kk][ii][jj]=255;
	}
	
  for (ii=0;ii<imageSize.x;ii++)
    for (jj=0;jj<imageSize.y;jj++)
      for (kk=minCoorVertical[ii][jj];kk<maxCoorVertical[ii][jj];kk++)
        maskForImagesWithHolesFilled[kk][ii][jj]=255;	  
  // --------------------------
  //   fill in all holes , end (Yangming added on Sept 14, 2010)
  // --------------------------
  
  
  
  //  -------------------------
  //  dilate
  // --------------------------
  printf("dilation radius = (%d,%d,%d)\n", dilateRadiusX, dilateRadiusY, dilateRadiusZ);
  
  
  for (kk=0;kk<imageSize.z;kk++)
    for (ii=0;ii<imageSize.x;ii++)
     for (jj=0;jj<imageSize.y;jj++)
	   {	   
	   maskForControlPoints[kk][ii][jj]=0;  // initialize
	   
	   if (maskForImagesWithHolesFilled[kk][ii][jj]>0)
	     maskForControlPoints[kk][ii][jj]=255;	     
	   else // detect if [x,y,z] is an edge point
	    {
	   chk1=0;
	   chk2=0;
	   
	   rmin=MAX(kk-1,0);
	   rmax=MIN(kk+1,imageSize.z-1);
	   smin=MAX(ii-1,0);
	   smax=MIN(ii+1,imageSize.x-1);
	   tmin=MAX(jj-1,0);
	   tmax=MIN(jj+1,imageSize.y-1);
	   for (r=rmin;r<=rmax;r++)
		for (s=smin;s<=smax;s++)
	     for (t=tmin;t<=tmax;t++)
		   {
		   if (maskForImagesWithHolesFilled[r][s][t]>0)   chk1=1;
		   if (maskForImagesWithHolesFilled[r][s][t]==0)  chk2=1;
		   }
	   
	   if ( (chk1==1)&&(chk2==1) )
	     {
		 rmin = MAX(0, kk-dilateRadiusZ);
		 rmax = MIN(imageSize.z-1, kk+dilateRadiusZ);
		 smin = MAX(0, ii-dilateRadiusX);
		 smax = MIN(imageSize.x-1, ii+dilateRadiusX);
		 tmin = MAX(0, jj-dilateRadiusY);
		 tmax = MIN(imageSize.y-1, jj+dilateRadiusY);
		 for (r=rmin;r<=rmax;r++)
		  for (s=smin;s<=smax;s++)
		    for (t=tmin;t<=tmax;t++)
		      maskForControlPoints[r][s][t]=255;			   
		 }
	   
	    }
	   }	
  
  printf("\ntotal number of foreground voxels = %d\n", count);
  
  UCfree3d(maskForImagesWithHolesFilled, imageSize.z, imageSize.x);
  
  return count;	   
}

// ---------------------------------------------------------------------------
void saveDisplacementAtControlPoints3D(Fvector3d ***controlPoints, int numControlPointsXThisLevel, int numControlPointsYThisLevel, int numControlPointsZThisLevel)
{
   int i, j, k;
   char fileNameX[1024], fileNameY[1024], fileNameZ[1024];
   FILE *fp1, *fp2, *fp3;
   
   sprintf(fileNameX, "%s", "dvXMat.txt");
   sprintf(fileNameY, "%s", "dvYMat.txt");
   sprintf(fileNameZ, "%s", "dvZMat.txt");
   
   printf("Saving displacement at control points...\n");
   fp1 = fopen(fileNameX, "w");
   fp2 = fopen(fileNameY, "w");
   fp3 = fopen(fileNameZ, "w");
   for (k=0;k<numControlPointsZThisLevel;k++)
    for (i=0;i<numControlPointsXThisLevel;i++)
     {
	 printf("i=%d\n", i);
     for (j=0;j<numControlPointsYThisLevel;j++)
	   {
	    fprintf(fp1, "%f ", controlPoints[k][i][j].y);
		fprintf(fp2, "%f ", controlPoints[k][i][j].x);
		fprintf(fp3, "%f ", controlPoints[k][i][j].z);
	   }
	 }
  
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}

// ---------------------------------------------------------------------------
void SmoothUCImage3D(unsigned char ***image, unsigned char ***imageSmoothed, Ivector3d imageSize)
{
  int i,j,k;
  float ***smoothMat;
  Ivector3d inputSize, smoothMatSize;
  float ***imageFloat, ***imageSmoothedFloat;
  int gaussianRadiusXY = 2;
  int gaussianRadiusZ = 1;
  float gaussianKernelXY = 1.4;
  float gaussianKernelZ = 0.6;
  
  if (imageSize.x>100 && imageSize.x<=200)
    {
	gaussianRadiusXY = 3;
	gaussianKernelXY = 2;
	}
  else if (imageSize.x>200)
    {
	gaussianRadiusXY = 4;
	gaussianKernelXY = 2.5;
	}
  
  inputSize.x = imageSize.x+2*gaussianRadiusXY;
  inputSize.y = imageSize.y+2*gaussianRadiusXY;
  inputSize.z = imageSize.z+2*gaussianRadiusZ;
  smoothMatSize.x = 2*gaussianRadiusXY+1;
  smoothMatSize.y = 2*gaussianRadiusXY+1;
  smoothMatSize.z = 2*gaussianRadiusZ+1;
  
  smoothMat = Falloc3d(smoothMatSize.x, smoothMatSize.y, smoothMatSize.z);
  gauss3D(smoothMatSize.x, gaussianKernelXY, smoothMatSize.y, gaussianKernelXY, smoothMatSize.z, gaussianKernelZ, smoothMat);
  
  imageFloat = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  imageSmoothedFloat = Falloc3d(inputSize.x, inputSize.y, inputSize.z);
  
  // display the smoothing matrix
  printf("smoothing matrix:");
  for (k=0;k<smoothMatSize.z;k++)
   {
   printf("k=%d\n",k);
   for (i=0;i<smoothMatSize.x;i++)
    {
	printf("\n");
	for (j=0;j<smoothMatSize.y;j++)
	  printf("%f ", smoothMat[k][i][j]);
	}
   }
  printf("\n\n");
	
	
  for (k=0;k<inputSize.z;k++)
   for (i=0;i<inputSize.x;i++)
    for (j=0;j<inputSize.y;j++)
	  {
	  if ( (i<gaussianRadiusXY) | (j<gaussianRadiusXY) | (k<gaussianRadiusZ) | (i>=(imageSize.x+gaussianRadiusXY)) | (j>=(imageSize.y+gaussianRadiusXY)) | (k>=(imageSize.z+gaussianRadiusZ)) )
	    {
		imageFloat[k][i][j] = 0.0;
		imageSmoothedFloat[k][i][j] = 0.0;
		}
	  else
	    {
		imageFloat[k][i][j] = (float)image[k-gaussianRadiusZ][i-gaussianRadiusXY][j-gaussianRadiusXY];
		imageSmoothedFloat[k][i][j] = 0.0;
	    }
	  }
	  
  Convolution3DFloat(imageFloat, inputSize, smoothMat, smoothMatSize, imageSmoothedFloat);
  
  for (k=0;k<imageSize.z;k++)
   for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
	  imageSmoothed[k][i][j] = (int)(imageSmoothedFloat[k+gaussianRadiusZ][i+gaussianRadiusXY][j+gaussianRadiusXY]+0.5);
	  
  // release momery
  Ffree3d(imageSmoothedFloat, inputSize.z, inputSize.x);
  Ffree3d(imageFloat, inputSize.z, inputSize.x);
  Ffree3d(smoothMat, smoothMatSize.z, smoothMatSize.x);
}

// ---------------------------------------------------------------------------
bool ReadFeaturesFromFeatureList(const char* featureImageListFile, unsigned char ****featureMap, const Image::Region& region)
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
  string featureDir   = os::path::dirname(featureImageListFile);
  char   filename[1024];

  for (featureIndex = 0; featureIndex < featureNum; featureIndex++) {
    if (fscanf(fp, "%s", filename) != 1) {
      fprintf(stderr, "Failed to read filename of feature number %d from file %s!\n", featureIndex + 1, featureImageListFile);
      break;
    }
    string filepath = os::path::join(os::path::abspath(featureDir.c_str()), filename);
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
    if (featureimage->region != region) {
        Image* tmp = ResizeImage(featureimage, region, 0, true);
        if (tmp == NULL) {
            fprintf(stderr, "Failed to pad feature image %d read from file %s!\n", featureIndex + 1, filepath.c_str());
            break;
        }
        delete featureimage;
        featureimage = tmp;
    }
    for (int k = 0; k < region.nz; k++) {
      for (int i = 0; i < region.nx; i++) {
        for (int j = 0; j < region.ny; j++) {
          featureMap[featureIndex][k][i][j] = featureimage->img.uc[k][i][j];
        }
      }
    }
    delete featureimage;
    featureimage = NULL;
  }

  if (featureimage) delete featureimage;
  fclose(fp);

  return (featureIndex == featureNum); // all feature images read successfully
}

// ---------------------------------------------------------------------------
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPoint(unsigned char ****featureMapA, unsigned char ****featureMapB, float ***confidenceMap, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defField, int numFeatures, float ****featureMapA2BUpdated, float ****featureMapA2BOld, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, float increment, char orientation, int NthLevel, int updateFeaturesOrNot)
{
  int i,j,k, featureIndex;
  int xcoor, ycoor, zcoor;
  float deltaEnergySquared;
  float energyChange=0;
  
  for (k=-affectingRadiusZAhead;k<affectingRadiusZBehind;k++)
   for (i=-affectingRadiusXAhead;i<affectingRadiusXBehind;i++)
	for (j=-affectingRadiusYAhead;j<affectingRadiusYBehind;j++)
	  {
		xcoor = realCoorThisControlPoint.x + j;
		ycoor = realCoorThisControlPoint.y + i;
		zcoor = realCoorThisControlPoint.z + k;
		
		if (xcoor>=0 && xcoor<=(imageSize.x-1) && ycoor>=0 && ycoor<=(imageSize.y-1) && zcoor>=0 && zcoor<imageSize.z)
		  {
		  deltaEnergySquared = 0.0;
		  for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
			{
			featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor] = FFDLocalEffect(featureMapA[featureIndex], xcoor, ycoor, zcoor, controlPoints, defField, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPointIndex1, controlPointIndex2, controlPointIndex3, increment, orientation, NNO, featureIndex, numFeatures);
			deltaEnergySquared += (featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor]+featureMapA2BOld[featureIndex][zcoor][xcoor][ycoor]-2*(float)featureMapB[featureIndex][zcoor][xcoor][ycoor])*(featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor]-featureMapA2BOld[featureIndex][zcoor][xcoor][ycoor]);
			
			if (updateFeaturesOrNot==YYES)
			  featureMapA2BOld[featureIndex][zcoor][xcoor][ycoor] = featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor];
			}
		  energyChange += deltaEnergySquared * confidenceMap[zcoor][xcoor][ycoor];  // confidence map involved   
		  }
	  } // calculate energy change brought by deltaX
	
  energyChange /= numFeatures;
  
  return energyChange;
}

// ---------------------------------------------------------------------------
float UpdateFeaturesAndCalculateEnergyChangeInTheBlockInducedByThisControlPointTemp(unsigned char ****featureMapA, unsigned char ****featureMapB, float ***confidenceMap, Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d ***controlPoints, Fvector3d ***defField, int numFeatures, float ****featureMapA2BUpdated, float ****featureMapA2BOld, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind, int controlPointIndex1, int controlPointIndex2, int controlPointIndex3, float increment, char orientation, int NthLevel)
{
  int i,j,k, featureIndex;
  int xcoor, ycoor, zcoor;
  float deltaEnergySquared;
  float energyChange=0;
  
  // to speed up, sample the voxels in the neighborhood of control point (controlPointIndex1, controlPointIndex2)
  
  int startPointX = -(int)((float)affectingRadiusXAhead*pow(0.7, NthLevel));
  int endPointX = (int)((float)affectingRadiusXBehind*pow(0.7, NthLevel));
  int startPointY = -(int)((float)affectingRadiusYAhead*pow(0.7, NthLevel));
  int endPointY = (int)((float)affectingRadiusYBehind*pow(0.7, NthLevel));
  int startPointZ = -(int)((float)affectingRadiusZAhead*pow(0.8, NthLevel));
  int endPointZ = (int)((float)affectingRadiusZBehind*pow(0.8, NthLevel));
  int intervalX = (int)ceil((float)(endPointX-startPointX)/5.0);
  int intervalY = (int)ceil((float)(endPointY-startPointY)/5.0);
  int intervalZ = (int)ceil((float)(endPointZ-startPointZ)/2.0);
  int numPointsSampled=0;
  
  
  for (k=startPointZ;k<=endPointZ;k=k+intervalZ)
   for (i=startPointX;i<=endPointX;i=i+intervalX)
	for (j=startPointY;j<=endPointY;j=j+intervalY)
	  {
		xcoor = realCoorThisControlPoint.x + j;
		ycoor = realCoorThisControlPoint.y + i;
		zcoor = realCoorThisControlPoint.z + k;
		
		if (xcoor>=1 && xcoor<(imageSize.x-1) && ycoor>=1 && ycoor<(imageSize.y-1) && zcoor>=0 && zcoor<imageSize.z)
		  {
		  deltaEnergySquared = 0.0;
		  for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
			{
			featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor] = FFDLocalEffect(featureMapA[featureIndex], xcoor, ycoor, zcoor, controlPoints, defField, imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, controlPointIndex1, controlPointIndex2, controlPointIndex3, increment, orientation, YYES, featureIndex, numFeatures);
			deltaEnergySquared += (featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor]+featureMapA2BOld[featureIndex][zcoor][xcoor][ycoor]-2*(float)featureMapB[featureIndex][zcoor][xcoor][ycoor])*(featureMapA2BUpdated[featureIndex][zcoor][xcoor][ycoor]-featureMapA2BOld[featureIndex][zcoor][xcoor][ycoor]);
			}
		  energyChange += deltaEnergySquared * confidenceMap[zcoor][xcoor][ycoor];  // confidence map involved
		  numPointsSampled++;
		  }
	  } // calculate energy change brought by deltaX
	  
  float samplingRate = pow(1.15,NthLevel)*(float)distBetweenControlPointsX*(float)distBetweenControlPointsY*distBetweenControlPointsZ*64.0/(float)numPointsSampled;	  
  energyChange *= (samplingRate/(float)numFeatures);

  return energyChange;
}

// ---------------------------------------------------------------------------
void CopyFeatureMapFloatAtTheBlockInducedByThisControlPoint(float ****featureMapOrig, float ****featureMapCopied, Ivector3d imageSize, int numFeatures, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind)
{
  int i,j,k, featureIndex;
  int xcoor, ycoor, zcoor;
  
  for (k=-affectingRadiusZAhead;k<affectingRadiusZBehind;k++)
   for (i=-affectingRadiusXAhead;i<affectingRadiusXBehind;i++)
	for (j=-affectingRadiusYAhead;j<affectingRadiusYBehind;j++)
	  {
		xcoor = realCoorThisControlPoint.x + j;
		ycoor = realCoorThisControlPoint.y + i;
		zcoor = realCoorThisControlPoint.z + k;
		
		if (xcoor>=0 && xcoor<imageSize.x && ycoor>=0 && ycoor<imageSize.y && zcoor>=0 && zcoor<imageSize.z)
		  {
			for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
			  featureMapCopied[featureIndex][zcoor][xcoor][ycoor] = featureMapOrig[featureIndex][zcoor][xcoor][ycoor];
		  }
	  }
  
}

// ---------------------------------------------------------------------------
void CopyDispAtTheBlockInducedByThisControlPoint(Fvector3d ***defFieldOrig, Fvector3d ***defFieldCopied, Ivector3d imageSize, Ivector3d realCoorThisControlPoint, int affectingRadiusXAhead, int affectingRadiusXBehind, int affectingRadiusYAhead, int affectingRadiusYBehind, int affectingRadiusZAhead, int affectingRadiusZBehind)
{
  int i,j,k;
  int xcoor, ycoor, zcoor;
  
  for (k=-affectingRadiusZAhead;k<affectingRadiusZBehind;k++)
   for (i=-affectingRadiusXAhead;i<affectingRadiusYBehind;i++)
	for (j=-affectingRadiusYAhead;j<affectingRadiusYBehind;j++)
	  {
		xcoor = realCoorThisControlPoint.x + j;
		ycoor = realCoorThisControlPoint.y + i;
		zcoor = realCoorThisControlPoint.z + k;
		
		if (xcoor>=0 && xcoor<imageSize.x && ycoor>=0 && ycoor<imageSize.y && zcoor>=0 && zcoor<imageSize.z)
		    {
			defFieldCopied[zcoor][xcoor][ycoor].x = defFieldOrig[zcoor][xcoor][ycoor].x;
			defFieldCopied[zcoor][xcoor][ycoor].y = defFieldOrig[zcoor][xcoor][ycoor].y;
			defFieldCopied[zcoor][xcoor][ycoor].z = defFieldOrig[zcoor][xcoor][ycoor].z;
			}
	  } 
}

// ---------------------------------------------------------------------------
void initializeConfidenceMapWithOnes(float ***confidenceMap, Ivector3d imageSize)
{
  int i,j,k;
  
  #pragma omp parallel for shared(confidenceMap) private(i,j,k) num_threads(100)
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++)
	    confidenceMap[k][i][j] = 1.0;
}


// ---------------------------------------------------------------------------
void UpdateConfidenceMapByInputMask(float ***confidenceMap, unsigned char ***maskForImage, Ivector3d imageSize)
{
  int i,j,k;
  
  #pragma omp parallel for shared(confidenceMap) private(i,j,k) num_threads(100)
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++){
	    if (maskForImage[k][i][j]==0) confidenceMap[k][i][j] = 0.0; }
}


// ---------------------------------------------------------------------------
void calculateConfidenceMap(float ****featureMapA2BIterative, unsigned char ****featureMapB, unsigned char ***mask, Ivector3d imageSize, int numFeatures, int levelIndex, int distBetweenControlPointsX, int distBetweenControlPointsY, int SimilarityMeasure, float ***confidenceMap)
{
  // *************************
  //  REMEMBER: Confidence Map is Always Defined on the space of Image B
  // *************************
  int i,j,k;
  int xB,yB,zB;
  int xA2B,yA2B,zA2B;
  int r,s,t;
  int ss,tt;
  int o,p,q;
  int featureIndex;
  int counter;
  float sum_CN, sum_PN;
  float mean_CN, mean_PN;
  float confidenceTimesSim, maxConfidenceTimesSim;
  float confidence;
  float dist, sim;
  float maxConfidenceInImage=0.0;
  float maxSim;
  float simThre;
  int searchRadiusXY;  
  int radiusXY_CN;   // Core Neighborhood
  int radiusXY_TN;   // Transitional Neighborhood
  int radiusXY_PN;  // Peripheral Neighborhood
  
  if (levelIndex==3)  // coarsest level (resolution)
    {
	radiusXY_CN = 1;
	radiusXY_TN = 2;
	radiusXY_PN = 4;
	searchRadiusXY = (int)(0.5+MAX(distBetweenControlPointsX,distBetweenControlPointsY)*0.4);
	}
  else if (levelIndex==2)  // middel level (resolution)
    {
	radiusXY_CN = 1;
	radiusXY_TN = 3;
	radiusXY_PN = 5;
	searchRadiusXY = (int)(0.5+MAX(distBetweenControlPointsX,distBetweenControlPointsY)*0.6);
	}
  else if (levelIndex==1)  // finest level (resolution)
    {
	radiusXY_CN = 2;
	radiusXY_TN = 4;
	radiusXY_PN = 6;
	searchRadiusXY = (int)(0.5+MAX(distBetweenControlPointsX,distBetweenControlPointsY)*0.8);
	}
  
  int matRadiusXY = searchRadiusXY+radiusXY_PN;
  Matrix *matSim;
  CreateMatrix(&matSim, 2*matRadiusXY+1, 2*matRadiusXY+1);
  
  
  if (SimilarityMeasure==SSD)
  {
  #pragma omp parallel for private(xB,yB,zB,confidence,maxSim,s,t,xA2B,yA2B,zA2B,dist,sim,maxConfidenceTimesSim,simThre,sum_CN,mean_CN,sum_PN,mean_PN,counter,p,q,confidenceTimesSim) num_threads(100)
  for (zB=0;zB<imageSize.z;zB++)
   { //for (zB=0;zB<imageSize.z;zB++)
   printf("\tslice %d\n", zB);
   for (xB=0;xB<imageSize.x;xB++)
    for (yB=0;yB<imageSize.y;yB++)
	  {
	  confidence = 1.0;
	  confidenceMap[zB][xB][yB]=1.0;  // initialize
	  
	  if (mask[zB][xB][yB]>0)
	    {
		// calculate similarity between point (xB,yB,zB) on image B and every point (xA2B,yA2B,zA2B) within the (searchRadiuxXY+radiuXY_PN) on image A2BIterative
		maxSim=0.0;
		for (s=-matRadiusXY;s<=matRadiusXY;s++)
		  for (t=-matRadiusXY;t<=matRadiusXY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB;
			
			ss = s+matRadiusXY;
			tt = t+matRadiusXY;
			matSim->data[ss][tt] = 0.0;  // initialize
			
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) )
			  {
			  dist = calculateEuclideanDistanceBetweenTwoFeatureVectors(featureMapA2BIterative, featureMapB, xA2B, yA2B, zA2B, xB, yB, zB, numFeatures);
			  sim = 2.0/(2.0+dist);  // similarity between (xB,yB) and (xA2B,yA2B)                    // Note: Yangming changed 1.0/(1.0+dist) into 5.0/(5.0+dist) on 11/27/2011
			  matSim->data[ss][tt] = sim;
			  
			  if ( (sim>maxSim)&&(abs(s)<=searchRadiusXY)&&(abs(t)<=searchRadiusXY) )   maxSim=sim;
			  }
			}  // for...s...for...t...
			
		// find the best match (highest similarity*confidence) in image A2BIterative for point (xB,yB) in image B 
		maxConfidenceTimesSim = 0.0;
		simThre = 0.99*maxSim;
		for (s=-searchRadiusXY;s<=searchRadiusXY;s++)
		  for (t=-searchRadiusXY;t<=searchRadiusXY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB;
			
			ss = s+matRadiusXY;
			tt = t+matRadiusXY;
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) && (matSim->data[ss][tt]>simThre) )
			  {
			  // need to calculate confidence metric between (xA2B,yA2B) in image A2BIterative and (xB,yB) in image B
			  
			  // step 1: calculate mean_CN
			  sum_CN=0.0;
			  mean_CN=0.0;
			  counter = 0;
			  for (p=-radiusXY_CN;p<=radiusXY_CN;p++)
			    for (q=-radiusXY_CN;q<=radiusXY_CN;q++)
				  {
				  if ( (xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_CN += matSim->data[ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_CN = sum_CN/counter;
				
			  // step 2: calculate mean_PN
			  sum_PN = 0.0;
			  mean_PN = 0.0;
			  counter = 0;
			  for (p=-radiusXY_PN;p<=radiusXY_PN;p++)
			    for (q=-radiusXY_PN;q<=radiusXY_PN;q++)
				  {
				  if ( (abs(p)>radiusXY_TN)&(abs(q)>radiusXY_TN)&(xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_PN += matSim->data[ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_PN = sum_PN/counter;
				
			  // sim and confidence			  
			  if (mean_PN!=0)
			    //confidence = mean_CN/mean_PN;
			    confidence = MAX(0.0, 5*log(mean_CN/mean_PN));
			    
				
			  // test if sim*confidence reaches a high peak between (xB,yB) in image B and (xA2B,yA2B)=(xB+s,yB+t) in image A
			  confidenceTimesSim = confidence*matSim->data[ss][tt];
			  if (confidenceTimesSim>maxConfidenceTimesSim)
				{
			    maxConfidenceTimesSim = confidenceTimesSim;				
				confidenceMap[zB][xB][yB] = confidence;   
				}
			  if (confidence > maxConfidenceInImage)
			    maxConfidenceInImage=confidence;
			  } // if (xA2B,yA2B,zA2B) is within the image A2BIterative
			} // for...s...for...t
		} // if (mask[zB][xB][yB]>0)
	  } // for...xB...for...yB...for...zB...
   } // for (zB=0;zB<imageSize.z;zB++)
  } // if (SimilarityMeasure==SSD)

  if (SimilarityMeasure==CC)
  {
  #pragma omp parallel for private(xB,yB,zB,confidence,maxSim,s,t,xA2B,yA2B,zA2B,dist,sim,maxConfidenceTimesSim,simThre,sum_CN,mean_CN,sum_PN,mean_PN,counter,p,q,confidenceTimesSim) num_threads(100)
  for (zB=0;zB<imageSize.z;zB++)
   { //for (zB=0;zB<imageSize.z;zB++)
   printf("\tslice %d\n", zB);
   for (xB=0;xB<imageSize.x;xB++)
    for (yB=0;yB<imageSize.y;yB++)
	  {
	  confidence = 1.0;
	  confidenceMap[zB][xB][yB]=1.0;  // initialize
	  
	  if (mask[zB][xB][yB]>0)
	    {
		// calculate similarity between point (xB,yB,zB) on image B and every point (xA2B,yA2B,zA2B) within the (searchRadiuxXY+radiuXY_PN) on image A2BIterative
		maxSim=0.0;
		for (s=-matRadiusXY;s<=matRadiusXY;s++)
		  for (t=-matRadiusXY;t<=matRadiusXY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB;
			
			ss = s+matRadiusXY;
			tt = t+matRadiusXY;
			matSim->data[ss][tt] = 0.0;  // initialize
			
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) )
			  {
			  sim = calculateCorrelationCoefficientBetweenTwoFeatureVectors(featureMapA2BIterative, featureMapB, xA2B, yA2B, zA2B, xB, yB, zB, numFeatures);
			  matSim->data[ss][tt] = sim;
			  if (sim>maxSim)   maxSim=sim;
			  }
			}  // for...s...for...t...
			
		// find the best match (highest similarity*confidence) in image A2BIterative for point (xB,yB) in image B 
		maxConfidenceTimesSim = 0.0;
		//maxSim = matSim->data[matRadiusXY][matRadiusXY];
		simThre = 0.99*maxSim;
		for (s=-searchRadiusXY;s<=searchRadiusXY;s++)
		  for (t=-searchRadiusXY;t<=searchRadiusXY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB;
			
			ss = s+matRadiusXY;
			tt = t+matRadiusXY;
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) && (matSim->data[ss][tt]>simThre) )
			  {
			  // need to calculate confidence metric between (xA2B,yA2B) in image A2BIterative and (xB,yB) in image B
			  
			  // step 1: calculate mean_CN
			  sum_CN=0.0;
			  mean_CN=0.0;
			  counter = 0;
			  for (p=-radiusXY_CN;p<=radiusXY_CN;p++)
			    for (q=-radiusXY_CN;q<=radiusXY_CN;q++)
				  {
				  if ( (xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_CN += matSim->data[ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_CN = sum_CN/counter;
				
			  // step 2: calculate mean_PN
			  sum_PN = 0.0;
			  mean_PN = 0.0;
			  counter = 0;
			  for (p=-radiusXY_PN;p<=radiusXY_PN;p++)
			    for (q=-radiusXY_PN;q<=radiusXY_PN;q++)
				  {
				  if ( (abs(p)>radiusXY_TN)&(abs(q)>radiusXY_TN)&(xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_PN += matSim->data[ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_PN = sum_PN/counter;
				
			  // sim and confidence
			  
			  if (mean_PN!=0)
			    //confidence = mean_CN/mean_PN;
			    confidence = MAX(0.0, 5*log(mean_CN/mean_PN));
				
			    
			  // test if sim*confidence reaches a high peak between (xB,yB) in image B and (xA2B,yA2B)=(xB+s,yB+t) in image A
			  confidenceTimesSim = confidence*matSim->data[ss][tt];
			  if (confidenceTimesSim>maxConfidenceTimesSim)
			    {
			    maxConfidenceTimesSim = confidenceTimesSim;				
				confidenceMap[zB][xB][yB] = confidence; 
				}
			  if (confidence > maxConfidenceInImage)
			    maxConfidenceInImage=confidence;
			  } // if (xA2B,yA2B,zA2B) is within the image A2BIterative
			} // for...s...for...t
		} // if (mask[zB][xB][yB]>0)
	  } // for...xB...for...yB...for...zB...
   } // for (zB=0;zB<imageSize.z;zB++)
  } // if (SimilarityMeasure==SSD)
  printf("\nMax confidence in image (before smoothing) = %f\n\n", maxConfidenceInImage);
  FreeMatrix(matSim);



  // smooth the confidence map slice by slice
  printf("smooth mutual-saliency map...\n");
  float **smoothMat;
  Ivector2d smoothMatSize, sliceSize;
  sliceSize.x = imageSize.x;
  sliceSize.y = imageSize.y;
  smoothMatSize.x = static_cast<int>(ceil(7.0*pow(0.65,(levelIndex-1))));
  smoothMatSize.y = static_cast<int>(ceil(7.0*pow(0.65,(levelIndex-1))));
  if (smoothMatSize.x%2==0) smoothMatSize.x++;
  if (smoothMatSize.y%2==0) smoothMatSize.y++;
  
  smoothMat = Falloc2d(smoothMatSize.x, smoothMatSize.y);
  gauss2D(smoothMatSize.x, 3.5*pow(0.65, (levelIndex-1)), smoothMatSize.y, 2.0*pow(0.65, (levelIndex-1)), 0, smoothMat);
  float ***confidenceMapTemp;
  confidenceMapTemp = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
  
  maxConfidenceInImage=0.0;
  for (k=0;k<imageSize.z;k++)
   {
   Convolution2DFloat(confidenceMap[k], sliceSize, smoothMat, smoothMatSize, confidenceMapTemp[k]);
   for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
	  {
	  if (mask[k][i][j]==0)
	    confidenceMap[k][i][j]=1.0;
	  else
	    confidenceMap[k][i][j]=confidenceMapTemp[k][i][j];
		
	  if (confidenceMap[k][i][j]>maxConfidenceInImage)
		maxConfidenceInImage=confidenceMap[k][i][j];
	  }
   }
  printf("\nMax confidence in image (after smoothing) = %f\n\n", maxConfidenceInImage);
  
   
  //Ffree3d(smoothMat, smoothMatSize.z, smoothMatSize.x);
  Ffree2d(smoothMat, smoothMatSize.x);
  Ffree3d(confidenceMapTemp, imageSize.z, imageSize.x);  
}

// ---------------------------------------------------------------------------
void calculateConfidenceMap3D(float ****featureMapA2BIterative, unsigned char ****featureMapB, unsigned char ***mask, Ivector3d imageSize, int numFeatures, int levelIndex, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int SimilarityMeasure, float ***confidenceMap)
{
  // *************************
  //  REMEMBER: Confidence Map is Always Defined on the space of Image B
  // *************************
  int i,j,k;
  int xB,yB,zB;
  int xA2B,yA2B,zA2B;
  int r,s,t;
  int rr,ss,tt;
  int o,p,q;
  int featureIndex;
  int counter;
  float sum_CN, sum_PN;
  float mean_CN, mean_PN;
  float confidenceTimesSim, maxConfidenceTimesSim;
  float confidence;
  float dist, sim;
  float maxConfidenceInImage=0.0;
  float maxSim;
  float simThre;
  int searchRadiusX, searchRadiusY, searchRadiusZ;  
  int radiusXY_CN;   // Core Neighborhood
  int radiusXY_TN;   // Transitional Neighborhood
  int radiusXY_PN;  // Peripheral Neighborhood
  
  if (levelIndex==3)  // coarsest level (resolution)
    {
	radiusXY_CN = 1;
	radiusXY_TN = 2;
	radiusXY_PN = 4;
	searchRadiusX = (int)(0.5+distBetweenControlPointsX*0.4);
	searchRadiusY = (int)(0.5+distBetweenControlPointsY*0.4);
	searchRadiusZ = (int)(0.5+distBetweenControlPointsZ*0.4);
	}
  else if (levelIndex==2)  // middel level (resolution)
    {
	radiusXY_CN = 1;
	radiusXY_TN = 3;
	radiusXY_PN = 5;
	searchRadiusX = (int)(0.5+distBetweenControlPointsX*0.6);
	searchRadiusY = (int)(0.5+distBetweenControlPointsY*0.6);
	searchRadiusZ = (int)(0.5+distBetweenControlPointsZ*0.4);
	}
  else if (levelIndex==1)  // finest level (resolution)
    {
	radiusXY_CN = 2;
	radiusXY_TN = 4;
	radiusXY_PN = 6;
	searchRadiusX = (int)(0.5+distBetweenControlPointsX*0.8);
	searchRadiusY = (int)(0.5+distBetweenControlPointsY*0.8);
	searchRadiusZ = (int)(0.5+distBetweenControlPointsZ*0.6);
	}
  
  int matRadiusX = searchRadiusX+radiusXY_PN;
  int matRadiusY = searchRadiusY+radiusXY_PN;
  int matRadiusZ = searchRadiusZ;
  float ***matSim;
  matSim=Falloc3d(2*matRadiusX+1, 2*matRadiusY+1, 2*matRadiusZ+1);
  
  
  if (SimilarityMeasure==SSD)
  {
  #pragma omp parallel for private(xB,yB,zB,confidence,maxSim,s,t,xA2B,yA2B,zA2B,dist,sim,maxConfidenceTimesSim,simThre,sum_CN,mean_CN,sum_PN,mean_PN,counter,p,q,confidenceTimesSim) num_threads(100)
  for (zB=0;zB<imageSize.z;zB++)
   { //for (zB=0;zB<imageSize.z;zB++)
   printf("\tslice %d\n", zB);
   for (xB=0;xB<imageSize.x;xB++)
    for (yB=0;yB<imageSize.y;yB++)
	  {
	  confidence = 1.0;
	  confidenceMap[zB][xB][yB]=1.0;  // initialize
	  
	  if (mask[zB][xB][yB]>0)
	    {
		// calculate similarity between point (xB,yB,zB) on image B and every point (xA2B,yA2B,zA2B) within the (searchRadiuxXY+radiuXY_PN) on image A2BIterative
		maxSim=0.0;
		for (r=-matRadiusZ; r<=matRadiusZ; r++)
		for (s=-matRadiusX;s<=matRadiusX;s++)
		  for (t=-matRadiusY;t<=matRadiusY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB+r;
			
			rr = r+matRadiusZ;
			ss = s+matRadiusX;
			tt = t+matRadiusY;
			matSim[rr][ss][tt] = 0.0;  // initialize
			
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) )
			  {
			  dist = calculateEuclideanDistanceBetweenTwoFeatureVectors(featureMapA2BIterative, featureMapB, xA2B, yA2B, zA2B, xB, yB, zB, numFeatures);
			  sim = 10.0/(10.0+dist);  // similarity between (xB,yB) and (xA2B,yA2B)                    // Note: Yangming changed 1.0/(1.0+dist) into 5.0/(5.0+dist) on 11/27/2011
			  matSim[rr][ss][tt] = sim;
			  
			  if (sim>maxSim)   maxSim=sim;
			  }
			}  // for...s...for...t...
			
		// find the best match (highest similarity*confidence) in image A2BIterative for point (xB,yB) in image B 
		maxConfidenceTimesSim = 0.0;
		simThre = 0.99*maxSim;
		for (r=-searchRadiusZ; r<=searchRadiusZ; r++)
		for (s=-searchRadiusX;s<=searchRadiusX;s++)
		  for (t=-searchRadiusY;t<=searchRadiusY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB+r;
			
			rr = r+matRadiusZ;
			ss = s+matRadiusX;
			tt = t+matRadiusY;
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) && (matSim[rr][ss][tt]>simThre) )
			  {
			  // need to calculate confidence metric between (xA2B,yA2B) in image A2BIterative and (xB,yB) in image B
			  
			  // step 1: calculate mean_CN
			  sum_CN=0.0;
			  mean_CN=0.0;
			  counter = 0;
			  for (p=-radiusXY_CN;p<=radiusXY_CN;p++)
			    for (q=-radiusXY_CN;q<=radiusXY_CN;q++)
				  {
				  if ( (xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_CN += matSim[rr][ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_CN = sum_CN/counter;
				
			  // step 2: calculate mean_PN
			  sum_PN = 0.0;
			  mean_PN = 0.0;
			  counter = 0;
			  for (p=-radiusXY_PN;p<=radiusXY_PN;p++)
			    for (q=-radiusXY_PN;q<=radiusXY_PN;q++)
				  {
				  if ( (abs(p)>radiusXY_TN)&(abs(q)>radiusXY_TN)&(xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_PN += matSim[rr][ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_PN = sum_PN/counter;
				
			  // sim and confidence
			  
			  if (mean_PN!=0)
			    confidence = mean_CN/mean_PN;
				
				
			  // test if sim*confidence reaches a high peak between (xB,yB) in image B and (xA2B,yA2B)=(xB+s,yB+t) in image A
			  confidenceTimesSim = confidence*matSim[rr][ss][tt];
			  if (confidenceTimesSim>maxConfidenceTimesSim)
			    {
			    maxConfidenceTimesSim = confidenceTimesSim;				
				confidenceMap[zB][xB][yB] = confidence;   
				}
			  if (confidence > maxConfidenceInImage)
			    maxConfidenceInImage=confidence;
			  } // if (xA2B,yA2B,zA2B) is within the image A2BIterative
			} // for...s...for...t
		} // if (mask[zB][xB][yB]>0)
	  } // for...xB...for...yB...for...zB...
   } // for (zB=0;zB<imageSize.z;zB++)
  } // if (SimilarityMeasure==SSD)

  if (SimilarityMeasure==CC)
  {
  #pragma omp parallel for private(xB,yB,zB,confidence,maxSim,s,t,xA2B,yA2B,zA2B,dist,sim,maxConfidenceTimesSim,simThre,sum_CN,mean_CN,sum_PN,mean_PN,counter,p,q,confidenceTimesSim) num_threads(100)
  for (zB=0;zB<imageSize.z;zB++)
   { //for (zB=0;zB<imageSize.z;zB++)
   printf("\tslice %d\n", zB);
   for (xB=0;xB<imageSize.x;xB++)
    for (yB=0;yB<imageSize.y;yB++)
	  {
	  confidence = 1.0;
	  confidenceMap[zB][xB][yB]=1.0;  // initialize
	  
	  if (mask[zB][xB][yB]>0)
	    {
		// calculate similarity between point (xB,yB,zB) on image B and every point (xA2B,yA2B,zA2B) within the (searchRadiuxXY+radiuXY_PN) on image A2BIterative
		maxSim=0.0;
		for (r=-matRadiusZ; r<=matRadiusZ; r++)
		 for (s=-matRadiusX; s<=matRadiusX;s++)
		  for (t=-matRadiusY; t<=matRadiusY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB+r;
			
			rr = r+matRadiusZ;
			ss = s+matRadiusX;
			tt = t+matRadiusY;
			matSim[rr][ss][tt] = 0.0;  // initialize
			
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) )
			  {
			  sim = calculateCorrelationCoefficientBetweenTwoFeatureVectors(featureMapA2BIterative, featureMapB, xA2B, yA2B, zA2B, xB, yB, zB, numFeatures);
			  matSim[rr][ss][tt] = sim;
			  if (sim>maxSim)   maxSim=sim;
			  }
			}  // for...s...for...t...
			
		// find the best match (highest similarity*confidence) in image A2BIterative for point (xB,yB) in image B 
		maxConfidenceTimesSim = 0.0;
		//maxSim = matSim->data[matRadiusXY][matRadiusXY];
		simThre = 0.99*maxSim;
		for (r=-searchRadiusZ; r<=searchRadiusZ; r++)
		for (s=-searchRadiusX;s<=searchRadiusX;s++)
		  for (t=-searchRadiusY;t<=searchRadiusY;t++)
		    {
			xA2B = xB+s;
			yA2B = yB+t;
			zA2B = zB+r;
			
			rr = r+matRadiusZ;
			ss = s+matRadiusX;
			tt = t+matRadiusY;
			if ( (xA2B>=0)&(xA2B<imageSize.x)&(yA2B>=0)&(yA2B<imageSize.y)&(zA2B>=0)&(zA2B<imageSize.z) && (matSim[rr][ss][tt]>simThre) )
			  {
			  // need to calculate confidence metric between (xA2B,yA2B) in image A2BIterative and (xB,yB) in image B
			  
			  // step 1: calculate mean_CN
			  sum_CN=0.0;
			  mean_CN=0.0;
			  counter = 0;
			  for (p=-radiusXY_CN;p<=radiusXY_CN;p++)
			    for (q=-radiusXY_CN;q<=radiusXY_CN;q++)
				  {
				  if ( (xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_CN += matSim[rr][ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_CN = sum_CN/counter;
				
			  // step 2: calculate mean_PN
			  sum_PN = 0.0;
			  mean_PN = 0.0;
			  counter = 0;
			  for (p=-radiusXY_PN;p<=radiusXY_PN;p++)
			    for (q=-radiusXY_PN;q<=radiusXY_PN;q++)
				  {
				  if ( (abs(p)>radiusXY_TN)&(abs(q)>radiusXY_TN)&(xA2B+p>=0)&(xA2B+p<imageSize.x)&(yA2B+q>=0)&(yA2B+q<imageSize.y) )
				    {
					sum_PN += matSim[rr][ss+p][tt+q];
					counter++;
					}
				  }
			  if (counter!=0)
			    mean_PN = sum_PN/counter;
				
			  // sim and confidence
			  
			  if (mean_PN!=0)
			    confidence = mean_CN/mean_PN;
				
			  // test if sim*confidence reaches a high peak between (xB,yB) in image B and (xA2B,yA2B)=(xB+s,yB+t) in image A
			  confidenceTimesSim = confidence*matSim[rr][ss][tt];
			  if (confidenceTimesSim>maxConfidenceTimesSim)
			    {
				maxConfidenceTimesSim = confidenceTimesSim;				
				confidenceMap[zB][xB][yB] = confidence; 
				}
			  if (confidence > maxConfidenceInImage)
			    maxConfidenceInImage=confidence;
			  } // if (xA2B,yA2B,zA2B) is within the image A2BIterative
			} // for...s...for...t
		} // if (mask[zB][xB][yB]>0)
	  } // for...xB...for...yB...for...zB...
   } // for (zB=0;zB<imageSize.z;zB++)
  } // if (SimilarityMeasure==SSD)
  printf("\nMax confidence in image (before smoothing) = %f\n\n", maxConfidenceInImage);
  //FreeMatrix(matSim);
  Ffree3d(matSim, 2*matRadiusZ+1, 2*matRadiusX+1);

  // smooth the confidence map slice by slice
  printf("smooth mutual-saliency map...\n");
  float **smoothMat;
  Ivector2d smoothMatSize, sliceSize;
  sliceSize.x = imageSize.x;
  sliceSize.y = imageSize.y;
  smoothMatSize.x = (int)ceil(5*pow(0.65,(levelIndex-1)));
  smoothMatSize.y = (int)ceil(5*pow(0.65,(levelIndex-1)));
  smoothMat = Falloc2d(smoothMatSize.x, smoothMatSize.y);
  gauss2D(smoothMatSize.x, 2.5*pow(0.65, (levelIndex-1)), smoothMatSize.y, 2.0*pow(0.65, (levelIndex-1)), 0, smoothMat);
  float ***confidenceMapTemp;
  confidenceMapTemp = Falloc3d(imageSize.x, imageSize.y, imageSize.z);

  maxConfidenceInImage=0.0;
  for (k=0;k<imageSize.z;k++)
   {
   Convolution2DFloat(confidenceMap[k], sliceSize, smoothMat, smoothMatSize, confidenceMapTemp[k]);
   for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
	  {
	  if (mask[k][i][j]==0)
	    confidenceMap[k][i][j]=1.0;
	  else
	    confidenceMap[k][i][j]=confidenceMapTemp[k][i][j];
		
	  if (confidenceMap[k][i][j]>maxConfidenceInImage)
		maxConfidenceInImage=confidenceMap[k][i][j];
	  }
   }
  printf("\nMax confidence in image (after smoothing) = %f\n\n", maxConfidenceInImage);

  Ffree2d(smoothMat, smoothMatSize.x);
  Ffree3d(confidenceMapTemp, imageSize.z, imageSize.x);
}

// ---------------------------------------------------------------------------
// calculate distance (difference) between feature vectors
float calculateEuclideanDistanceBetweenTwoFeatureVectors(float ****featureMapA2B, unsigned char ****featureMapB, int xA, int yA, int zA, int xB, int yB, int zB, int numFeatures)
{
   float dist=0.0;
   int featureIndex;
   
   #pragma omp parallel for private(featureIndex) reduction(+: dist) num_threads(numFeatures)
   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
     {
      dist += pow((featureMapA2B[featureIndex][zA][xA][yA]-(float)featureMapB[featureIndex][zB][xB][yB]), 2.0);
	 }
   
   dist = sqrt(dist/(float)numFeatures);
   
   return dist;
}

// ---------------------------------------------------------------------------
// calculate correlation coefficient between feature vectors
float calculateCorrelationCoefficientBetweenTwoFeatureVectors(float ****featureMapA2B, unsigned char ****featureMapB, int xA, int yA, int zA, int xB, int yB, int zB, int numFeatures)
{
   float cc=0.0;
   int featureIndex;
   float sumA=0.0;
   float sumB=0.0;
   float sumAB=0.0;
   float sumAsq=0.0;
   float sumBsq=0.0;
   float aveA=0.0;
   float aveB=0.0;
   float covAB=0.0;
   float varA=0.0;
   float varB=0.0;
   
   #pragma omp parallel for private(featureIndex) reduction(+: dist) num_threads(numFeatures)
   for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
     {
     sumA += featureMapA2B[featureIndex][zA][xA][yA];
	 sumB += (float)featureMapB[featureIndex][zA][xA][yA];
	 sumAB += featureMapA2B[featureIndex][zA][xA][yA]*(float)featureMapB[featureIndex][zA][xA][yA];
	 sumAsq += featureMapA2B[featureIndex][zA][xA][yA]*featureMapA2B[featureIndex][zA][xA][yA];
	 sumBsq += (float)featureMapB[featureIndex][zA][xA][yA]*(float)featureMapB[featureIndex][zA][xA][yA];
	 }
   
   aveA = sumA/(float)numFeatures;
   aveB = sumB/(float)numFeatures;
   
   covAB = sumAB-aveB*sumA-aveA*sumB+(float)numFeatures*aveA*aveB;
   varA  = sumAsq-2.0*aveA*sumA+(float)numFeatures*aveA*aveA;
   varB  = sumBsq-2.0*aveB*sumB+(float)numFeatures*aveB*aveB;
   
   
   cc = covAB/sqrt(varA*varB);
   
   return cc;
}

// ---------------------------------------------------------------------------
void ApplyDeformationFieldOnImage(unsigned char ***image, Fvector3d ***dfFieldB2A, Ivector3d imageSize, unsigned char ***imageA2B)
{
   int i,j,k;
   float dx,dy,dz;
   float Fintensity;
   
   printf("Warping image by tentative deformations...\n\n");
   
   #pragma omp parallel for private(i,j,k,dx,dy,dz,Fintensity) num_threads(100)
   for (k=0;k<imageSize.z;k++)
     {
     for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 dx = dfFieldB2A[k][i][j].x;
		 dy = dfFieldB2A[k][i][j].y;
		 dz = dfFieldB2A[k][i][j].z;
		 
		 Fintensity = trilinearInterpolation(image, (float)i+dx, (float)j+dy, (float)k+dz, imageSize);	
		 imageA2B[k][i][j] = (int)Fintensity;
		 }
	 }
}

// ---------------------------------------------------------------------------
void ApplyDeformationFieldOnFeatures(unsigned char ****featureMapA, unsigned char ****featureMapB, int numFeatures, Ivector3d imageSize, Fvector3d ***dfFieldB2A, float ****featureMapA2BInitialized)
{
   int i,j,k;
   int featureIndex;
   float dx,dy,dz;
   float newFeature;
   
   int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
   float x_src, y_src, z_src;
   float w1, w2, w3, w4, w5, w6, w7, w8;
   float distXtoXFloor,distXtoXCeil,distYtoYFloor,distYtoYCeil,distZtoZFloor,distZtoZCeil;
  
  
   printf("Warping image by tentative deformations...\n\n");
   #pragma omp parallel for private(i,j,k,dx,dy,dz,x_src,y_src,z_src,xFloor,xCeil,yFloor,yCeil,zFloor,zCeil,distXtoXFloor,distXtoXCeil,distYtoYFloor,distYtoYCeil,distZtoZFloor,distZtoZCeil,w1,w2,w3,w4,w5,w6,w7,w8,featureIndex,newFeature)  num_threads(20)
   for (k=0;k<imageSize.z;k++)
     { //for (k=0;k<imageSize.z;k++)
	 for (i=0;i<imageSize.x;i++)
	   for (j=0;j<imageSize.y;j++)
	     {
		 dx = dfFieldB2A[k][i][j].x;
		 dy = dfFieldB2A[k][i][j].y;
		 dz = dfFieldB2A[k][i][j].z;
		 
		 
		 x_src = (float)i+dx;
		 y_src = (float)j+dy;
		 z_src = (float)k+dz;
	
		 if (x_src<0)  x_src=0.0;
		 if (y_src<0)  y_src=0.0;
		 if (z_src<0)  z_src=0.0;
		 if (x_src>imageSize.x-1) x_src=(float)(imageSize.x-1);
		 if (y_src>imageSize.y-1) y_src=(float)(imageSize.y-1);
		 if (z_src>imageSize.z-1) z_src=(float)(imageSize.z-1);
	
		 xFloor = (int)floor(x_src);
		 xCeil  = (int)ceil(x_src);
		 yFloor = (int)floor(y_src);
		 yCeil  = (int)ceil(y_src);
		 zFloor = (int)floor(z_src);
		 zCeil  = (int)ceil(z_src);
	
		 if ( (xFloor==xCeil)&(xFloor!=(imageSize.x-1)) )  xCeil+=1;
		 if ( (xFloor==xCeil)&(xFloor==(imageSize.x-1)) )  xFloor-=1;
		 if ( (yFloor==yCeil)&(yFloor!=(imageSize.y-1)) )  yCeil+=1;
		 if ( (yFloor==yCeil)&(yFloor==(imageSize.y-1)) )  yFloor-=1;
		 if ( (zFloor==zCeil)&(zFloor!=(imageSize.z-1)) )  zCeil+=1;
		 if ( (zFloor==zCeil)&(zFloor==(imageSize.z-1)) )  zFloor-=1;
	
		 distXtoXFloor = x_src-xFloor;
		 distXtoXCeil  = xCeil-x_src;
		 distYtoYFloor = y_src-yFloor;
		 distYtoYCeil  = yCeil-y_src;
		 distZtoZFloor = z_src-zFloor;
		 distZtoZCeil  = zCeil-z_src;
	
		 w1 = distXtoXFloor*distYtoYFloor*distZtoZFloor;
		 w2 = distXtoXFloor*distYtoYFloor*distZtoZCeil;
		 w3 = distXtoXFloor*distYtoYCeil*distZtoZFloor;
		 w4 = distXtoXFloor*distYtoYCeil*distZtoZCeil;
		 w5 = distXtoXCeil*distYtoYFloor*distZtoZFloor;
		 w6 = distXtoXCeil*distYtoYFloor*distZtoZCeil;
		 w7 = distXtoXCeil*distYtoYCeil*distZtoZFloor;
		 w8 = distXtoXCeil*distYtoYCeil*distZtoZCeil;
		 
		 for (featureIndex=0;featureIndex<numFeatures;featureIndex++)
		   {
		   //newFeature = trilinearInterpolation(featureMapA[featureIndex], (float)i+dx, (float)j+dy, (float)k+dz, imageSize);	
		   newFeature = (float)featureMapA[featureIndex][zCeil][xCeil][yCeil]   * w1 +
						(float)featureMapA[featureIndex][zFloor][xCeil][yCeil]  * w2 +
						(float)featureMapA[featureIndex][zCeil][xCeil][yFloor]  * w3 +
						(float)featureMapA[featureIndex][zFloor][xCeil][yFloor] * w4 +
						(float)featureMapA[featureIndex][zCeil][xFloor][yCeil]  * w5 +
						(float)featureMapA[featureIndex][zFloor][xFloor][yCeil] * w6 +
						(float)featureMapA[featureIndex][zCeil][xFloor][yFloor] * w7 +
						(float)featureMapA[featureIndex][zFloor][xFloor][yFloor]* w8 ;
		 
		   //if ( (newFeature<=0.0)&&(i+dx<0 || i+dx>=(imageSize.x-1) || j+dy<0 || j+dy>=(imageSize.y-1) || k+dz<0 || k+dz>=(imageSize.z-1)) )
		      //newFeature = MAX(0.0, (float)featureMapA[featureIndex][k][i][j]);
		  
		   featureMapA2BInitialized[featureIndex][k][i][j] = newFeature;
		   //if (featureMapA[featureIndex][z][x][y]!=0)
		   //printf("\tat (%d,%d,%d), #%d: %d->%2.2f, def = (%2.2f, %2.2f, %2.2f)\n", i,j,k,featureIndex, featureMapA[featureIndex][k][i][j], newFeature, dx, dy, dz);
		   }
		 }
	 } //for (k=0;k<imageSize.z;k++)
}



// ---------------------------------------------------------------------------
void IncorporateInitialDeformation(Fvector3d*** def, Fvector3d*** init_def, Ivector3d defsize, int levelIndex, int numImgLevels)
{
	int i,j,k;
	float lambda=(numImgLevels-levelIndex)*0.1;
	if (lambda<0.0) lambda=0.0;
	
	for (k=0;k<defsize.z;k++)
	  for (i=0;i<defsize.x;i++)
	    for (j=0;j<defsize.y;j++)
		 {
		  def[k][i][j].x = lambda*def[k][i][j].x + (1-lambda)*init_def[k][i][j].x;
		  def[k][i][j].y = lambda*def[k][i][j].y + (1-lambda)*init_def[k][i][j].y;
		  def[k][i][j].z = lambda*def[k][i][j].z + (1-lambda)*init_def[k][i][j].z;
		 }
}



// ---------------------------------------------------------------------------
float discreteOptimization(unsigned char ****SF, unsigned char ****TF, Fvector3d ***defField, int numSamples, Ivector3d imageSize, Fvector3d resolutionRatio, int levelIndex, int method, int distMethod, float threshold, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int numFeatures, int distBasedWeight, float regWeight,  float ***confidenceMap, unsigned char ***mask, Fvector3d ***controlPointDisp, int AdditionOrComposition, int chk, float *weightsX, float *weightsY, float *weightsZ, float label_factor, int iter,  int MaxNumIterInThisResolution, int indexGridLevel, int numGridLevels, int fastApproximationOrNot, int SimilarityMeasure, float initEnergyPrevIter)
{
	//FUNCTION INPUTS:
	//	SF: feature map of the source image
	//	TF: feature map of the target image
	//	detField: deformation field from the previous iterationf
	//  numnSamples: is the number of the samples that will be sampled along the principal diagonals and
	//				 the diagonals (18-neighborhood is considered in the 3D space)
	//	imageSize: is a vector that contains the size of the image
	//  resolutionRatio: is a vector that contains the ratio in resolution (voxel size)
	//	method: is the parameter that will define the size of the patch that will be used to calculate the singleton
	//			potential 
	//	distMethod: is the parameter that defines what pair-wise distance is going to be used.
	//	distBetweenControlPointsX: is the distance between the nodes of the deformation grid along the x axis
	//	distBetweenControlPointsY: is the distance between the nodes of the deformation grid along the y axis
	//  distBetweenControlPointsZ: is the distance between the nodes of the deformation grid along the z axis
	//	numFeatures: is the number of features used in order to calculate the feature based cost.
	//  distBasedWeight: if equal to 1, the contribution of each voxel to the control point of the deformation grid
	//					 will be a function of its distance to it. if equal to zero only the weights based on the
	//					 saliency map will be taken into account.
	//	regWeight: is the weight for the regularization term of the MRF energy (the weight with whom all edge potentials
	//			   are multiplied.
	//  threshold: is the theshold that is used for the pair-wise distance.
	//	confidenceMAP: is the map that gives the weights per point based on the computation of the saliency map
	//	mask: masks out the background.
	//
	//FUNCTION OUTPUT:
	//	controlPointDisp: contains the optimal displacement that should be applied to every node.

	//declaration of variables
	float *labels; //matrix that contains the labels that should be applied to the nodes
	int numPoints; //number of the nodes of the grids
	float *lcosts; //array of the size numLabels*numPoints containing the label costs (singleton potential)
	int *pairs;	   //array that contains the nodes' indices for each MRF edge.
	float *dist;   //the distance function used for defining the MRF pairwise potential
	float *wcosts; //array containing the weights for each edge (or else how important the regularization
				   //will be.
	int *optimalLabels; //array that contains the optimal labels per node
	//number of control points in each direction
	int numControlPointsX;
	int numControlPointsY;
	int numControlPointsZ;
	int numPairs;
	int maxIters = 10; //maximum number of iterations for the discrete optimization algorithm
	int i; //index
	int numLabels;
	//end of declaration, time for some action
	
	float initEnergyThisIter;
	float regWeightThisIter;
	
	// by default, search in 18 directions in the discretized space (+x, -x, +y, -y, +z, -z, (+x)(+y), (+x)(-y), (-x)(+y), (-x)(-y), (+x)(+z), (+x)(-z), (-x)(+z), (-x)(-z), (+y)(+z), (+y)(-z), (-y)(+z), (-y)(-z) )
	// in fast approximation mode, search only in 6 directions (+x, -x, +y, -y, +z, -z)
	if (fastApproximationOrNot==NNO)	
		numLabels = 18*numSamples + 1;
	else	
		numLabels = 6*numSamples + 1;
		
	numControlPointsX = (int)ceil((float)imageSize.x/(float)distBetweenControlPointsX);
	numControlPointsY = (int)ceil((float)imageSize.y/(float)distBetweenControlPointsY);
	numControlPointsZ = (int)ceil((float)(imageSize.z-1)/(float)distBetweenControlPointsZ);
	numPairs = ((numControlPointsX-1)*numControlPointsY + (numControlPointsY-1)*numControlPointsX)*numControlPointsZ
		+ numControlPointsX*numControlPointsY*(numControlPointsZ-1);
	numPoints = numControlPointsX*numControlPointsY*numControlPointsZ;

	//attribute memory	
	optimalLabels = (int *)malloc(numPoints*sizeof(int));
	labels = (float *)malloc(3*numLabels*sizeof(float));
	pairs = (int *)malloc(2*numPairs*sizeof(int));
	wcosts = (float *)malloc(numPairs*sizeof(float));
	dist = (float *)malloc(numLabels*numLabels*sizeof(float));
	lcosts = (float *)calloc(numPoints*numLabels,sizeof(float));

	//computations
	printf("computing labels (discretized displacements)...\n");
	computeLabels(numSamples, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, resolutionRatio, levelIndex, method, labels, label_factor, iter, MaxNumIterInThisResolution, indexGridLevel, numGridLevels, fastApproximationOrNot);
	printf("computing pairs...\n");
	computePairs(numControlPointsX, numControlPointsY, numControlPointsZ, pairs);	
	printf("computing wcosts...\n");
	//computeWcosts(numPairs, regWeight, wcosts);
	printf("computing dists...\n");
	computeDist(labels, numLabels, resolutionRatio, distMethod, threshold, dist, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ);	
	printf("computing unary... (this takes some time, be patient~ )\n");
	computeUnary(SF, TF, labels, numControlPointsX, numControlPointsY, numControlPointsZ, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, method, distBasedWeight, imageSize, resolutionRatio, confidenceMap, defField, numLabels, numFeatures, SimilarityMeasure, chk, mask, lcosts, AdditionOrComposition, weightsX, weightsY, weightsZ);
	
	regWeightThisIter = regWeight * pow(0.25, 3-levelIndex) * pow(0.75,  iter) * pow(0.75, indexGridLevel-1);
	computeWcosts_Ou(numPairs, numPoints, numLabels, regWeightThisIter, lcosts, dist, wcosts);

	//time to optimize
	printf("\nPreparing Optimization...");
	CV_Fast_PD pd( numPoints, numLabels, lcosts,
	               numPairs, pairs, dist, maxIters,
				   wcosts );
	printf("\nOptimizing...\n");
	initEnergyThisIter = pd.run(initEnergyPrevIter);
	printf("free memory...\n");
	free(pairs);
	free(wcosts);
	free(dist);

	// need to calculate the MRF
	printf("done!\n");
	//save computed labels to nodes
	for (i=0; i<numPoints; i++)
	{
		//to every node the optimal label is attributed
		optimalLabels[i] = pd._pinfo[i].label;
	}
	//assign optimal deformations to the graph nodes
	labels2deformations(optimalLabels, labels, numControlPointsX, numControlPointsY, 
						 numControlPointsZ, controlPointDisp, numLabels);
	
	//release memory
	free(optimalLabels); 
	free(labels); 
	free(lcosts);

	return initEnergyThisIter;
}
