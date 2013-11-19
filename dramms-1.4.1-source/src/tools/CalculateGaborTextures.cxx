/**
 * @file  CalculateGaborTextures.cxx
 * @brief Calculate Gabor textures of an image.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
//#include <omp.h>

#include <sys/times.h>
#include <unistd.h>

#include <time.h>

#include <common/general.h>
#include <common/cres.h>
#include <common/mvcd.h>
#include <common/matrix.h>  
#include <common/imageio.h>

#include <dramms/basis.h> // exename(), print_contact()

// acceptable in .cxx file
using namespace std;
using namespace basis;
using namespace dramms;



#define  ORIENT_SCALE  24
#define  NumOfDimension 2
#define  AlongXAxis 1
#define  AlongYAxis 2
#define  AlongZAxis 3


int Orig_X, Orig_Y, Orig_Z;
int SubSampled_X,SubSampled_Y,SubSampled_Z;
float res_X,res_Y,res_Z;
int SmoothLongAxis,SmoothShortAxis;

int side, ScaleNUM,OrientationNUM,DCOrNot;
float Ul,Uh;
int StartScale,EndScale,StartOrient,EndOrient;

int Sub_X,Sub_Y,Sub_Z;

/* Guassian variation */
float CannyGuassianSigma ;
int UsePolarCoordinate ;

Fvector2d CenterOfProbe;
double Radius;

double *WeightingOnFrequency ;

unsigned char**** featureImage;
float ****featureImageFloat;
float ***minImaginaryHori, ***maxImaginaryHori, ***minRealHori, ***maxRealHori, ***minImaginaryVert, ***maxImaginaryVert, ***minRealVert, ***maxRealVert;

char ResultDirectory[1024];

int SmoothOrNot;
int DebugOrNot;
int HoriOrVert;
int DetectOrNot;
int regionGrowingOrNot;
int fillHoleOrNot;
int saveMaskOrNot;
int inputMaskOrNot;



// ===========================================================================
// help
// ===========================================================================

// ---------------------------------------------------------------------------
void print_help()
{ 
    string exec_name = exename();
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input_image> <output_GaborFeaturePrefix>" << endl;
    cout << endl;
    cout << "Description:" << endl;
    cout << "  This program calculates Gabor textures of an input image and outputs" << endl;
    cout << "  a number of multi-scale and multi-orientation feature images." << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
    cout << "  <input_image>                 File to calculate Gabor features of." << endl;
    cout << "  <output_GaborFeaturePrefix>   Prefix for multiple output feature images." << endl;
    cout << "                                The output Jacobian images have datatype DT_FLOAT." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
    cout << "  -s <int>                     Number of Gabor scales. (default: 2)" << endl;
    cout << "  -o <int>                     Number of Gabor orientations. (default: 6)" <<endl;
    cout << "  -u <float>,<float>           Frequency range of Gabor filter bank. (default: [0.2, 0.8])" <<endl;
    cout << "  -i <int>                     Foreground threshold for generating mask for feature calculation." << endl;
    cout << "                               (default: 12, on a 0-255 intensity scale)" << endl;
    cout << "  -L                           Take log(JacobianDet). All non-positive Jacobian determinant" << endl;
    cout << "                               values will be mapped to -100.0 in log(JacDet) calculation." << endl;
    cout << "                               (default: off)" << endl;
    cout << "  -m <int>,<int>               The size of smoothing filter (LongAxis, ShortAxis). (default: 3, 2)" << endl;
    cout << "  -w <int>                     Half size of Gabor window. (default: 15)" << endl;
    cout << "  -b <int>,<int>               Subsample ratio (xsubsample, zsubsample). (default: 1, 1)" << endl;
    cout << "  -t <int>                     Method to convert Gabor responses to attribute image values (0-255)" << endl;
    cout << "                                 1: translate response values" << endl;
    cout << "                                 2: linear scaling of response values" << endl;
    cout << "                               (default: 1)" << endl;
    cout << "  -x <int>                     Save given number of features of Gabor scales" << endl;
    cout << "                               (e.g., -x1 means to save features only at finest Gabor scale," << endl;
    cout << "                               -x2 will save features from the 2 most finest scales...)." << endl;
    cout << "                               (default: save all scales)." << endl;
    cout << "  -S                           Smooth original image. (default: no)" << endl;
    cout << "  -F                           Use DC part. (default: no)" << endl;
    cout << "  -e <float>                   Guassian variation. (default: 1.0)" << endl;
    cout << "  -I <char>                    Input mask image. Only Gabor features inside the mask are extracted." << endl;
    cout << "  -P                           Use polar coordinates. (default: no)" << endl;
    cout << "  -D                           Display intermediate result. (default: off)" << endl;
    cout << "  -R                           Region growing to remove regions containing voxels less than certain amount." << endl;
    cout << "                               (default: off)" << endl;
    cout << "  -V                           Fill in any holes when generating mask. (default: off)" << endl;
    cout << "  -M                           Save foreground mask. (default: off)" << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << exec_name <<" DF.def JacobianMap.img -s0 -S0 " << endl;
    cout << endl;
    print_contact();
}



// sub-functions
void generateMask(unsigned char***image, Ivector3d imageSize, unsigned char***mask, int foregroundThre, int fillHoleOrNot, int regionGrowingOrNot);
void GetGaborFeatureOfOneImage(double **Features_Real, double **Features_Imginary, unsigned char **UC_img, int image_X, int image_Y,  int side, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag) ;
void GaborFunction(Matrix *Gr, Matrix *Gi, int s, int n, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag) ;
void WriteGaborFiltersToFile(double *GaborFilterBank, int side, int ScaleNUM, int OrientationNUM) ;
void CalculateAndSaveGaborFilters_ForVisualEffect(int side, int ScaleNUM, int OrientationNUM, int flag, double Ul, double Uh) ;
void GetGaborFeatureOfOneImage_ThroughPolarCoordinate(double **Features_Real, double **Features_Imginary, unsigned char **UC_img, int image_X, int image_Y,  int side, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag) ;
void CalculateGaborRotationInvariants_ToProbeCenter(double **Features_Real, double **Features_Imginary, int image_X, int image_Y,  int ScaleNUM, int OrientationNUM ) ;

/* edge detection */
void canny_core(float s, int cols, int rows, int FeaNUM, double **Features_Real,double **Features_Imginary,unsigned char *derivative_mag,unsigned char *magnitude,unsigned char *orientation) ;
void DetectEdgeFromGaborFeatureVectors(unsigned char **Img_edge, int image_X, int image_Y, int FeaNUM, double **Features_Real, double **Features_Imginary, float CannyGuassianSigma) ;
void WriteImg2D(const char* filename, unsigned char **data, int image_X, int image_Y) ;
void Write3DImage(const char* filename,unsigned char*** Image,int x_size,int y_size,int z_size);

void SmoothingFilterOnTheOrientedEllipseDomain( int SmoothLongAxis, int SmoothShortAxis, unsigned char **UC_img, int image_X, int image_Y ) ;
void SmoothingImgByRemovingHighFrequencyFourierParemters(unsigned char **UC_img, int image_X, int image_Y) ;
void Put2DGaborResponseInto3DArray(double **Features_Real, double **Features_Imginary, unsigned char **mask, int image_X, int image_Y, int StartScale,int EndScale,int StartOrient,int EndOrient,int ScaleNUM, int OrientationNUM,int SliceIndex, int method);


/* polar */
void PolarizeImg( unsigned char **UC_img, int image_X, int image_Y ) ;

void DetectCircle(unsigned char** i_image,Fvector2d* center,double* radius);

void CalculateGaborFeatureFromOneSlice(unsigned char** i_image,unsigned char **mask, int image_X,int image_Y, float res_X, float res_Y, int RotateInvariantOrNot,int SliceIndex, int method);

double GaborMotherFrequencyFunction(double u,double v,double Uvar,double Vvar,double Shift);
void GaborFrequencyFunction(Matrix *Gr, Matrix *Gi, int Scale,int Orientation, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag);
void GetGaborFeatureOfOneImageNew(double **Features_Real, double **Features_Imginary, unsigned char **UC_img, unsigned char **mask, int image_X, int image_Y, float res_X, float res_Y, double Ul,double Uh, int StartScale,int EndScale,int StartOrient,int EndOrient,int ScaleNUM, int OrientationNUM, int flag);

void GaborSpatialFunction(Matrix *Gr, Matrix *Gi, int Scale,int Orientation, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag);
void GaborMotherSpatialFunction(double x,double y,double Xvar,double Yvar,double Shift,double* real,double* imag);

void ReadImg2D(const char* filename,unsigned char **data,int image_X,int image_Y);
void Construct3dImg(unsigned char*** F_img3d,int Scale,int Orientation,int image_X,int image_Y,int SliceNum,int RealOrImag);

double **DDalloc2d(int i_size,int j_size);
float GetGreatestCommonDivisor(float aa, float bb);
void zoomInImage(unsigned char **UC_Img, unsigned char **zoomedUC_Img, int origSizeX, int origSizeY, int newSizeX, int newSizeY, int zoomFactor_X, int zoomFactor_Y);


void linearScalingGaborResponseIntoAttributeValues(int scaleIndex, int orientationIndex, int ScaleNUM, int OrientationNUM, int SliceNUM, int sizeX, int sizeY, int sizeZ, unsigned char ***mask);
void RegionGrowingInBinaryImage(unsigned char ***mask, Ivector3d imageSize, int RegionPointNumThre);
void NeighboringSimilarPointsSearch( long int *PointNum, int ***Detected_Region, int order, unsigned char ***src, unsigned char ***status, unsigned char ***VoxelInStack, Ivector3d *Stack, int i, int j, int k, int Greythreshold, int x_size, int y_size, int z_size ) ;


// ===========================================================================
// IO
// ===========================================================================

// ---------------------------------------------------------------------------
bool WriteUCImage(const char* filename, const Image* image, unsigned char*** img)
{
    Image tmp;
    tmp.hdr      = image->hdr;
    tmp.compress = image->compress;
    tmp.filefmt  = image->filefmt;
    tmp.imgfmt   = image->imgfmt;
    tmp.img.uc   = img;
    tmp.ownsimg  = false;
    tmp.region   = image->region;
    return WriteNiftiImage(filename, &tmp);
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


   int    num ;
   FILE *fp;
   unsigned char ***UC_img3d, ***Fea_img3d,***Img_edge, ***mask3d;
   unsigned char **UC_img2d, **mask2d;
   int method;
   int foregroundThre;


    /* Gabor filters */
   int    i, j, k,temp ;
   int s,n;
   char filename[1024];
   char inputMaskName[1024];
   struct tms startTime,endTime;
   double durTime;
   Ivector3d imageSize3d;
   Ivector2d imageSize2d;

   SmoothLongAxis=3; SmoothShortAxis=2;
   side=15;
   ScaleNUM=2;
   OrientationNUM=6;
   Ul = 0.2 ;       
   Uh = 0.8 ;        
   Sub_X=1;
   Sub_Y=1;
   Sub_Z=1;
   method=1; // 1 for translate, 2 for linear scaling, of gabor responses into attribute values
   foregroundThre=12;
   EndScale=9999;
   DCOrNot=NNO;
   SmoothOrNot=NNO;
   DebugOrNot=NNO;
   HoriOrVert=AlongZAxis;
   fillHoleOrNot = NNO;
   saveMaskOrNot = NNO;
   regionGrowingOrNot = NNO;
 
   CannyGuassianSigma = 1.0 ;      
   UsePolarCoordinate = NNO ;
   inputMaskOrNot=NNO;
   DetectOrNot=YYES;

   // parse arguments
   int c = -1;
   while ((c=getopt(argc, argv, "m:w:s:o:u:b:Fe:t:x:i:I:PDSRVM")) != -1) {
      switch(c)
        {
        case 'm':
          sscanf(optarg, "%d,%d", &SmoothLongAxis, &SmoothShortAxis) ;
          break;
          
        case 'w':
          sscanf(optarg,"%d",&side);
          break;

        case 's':
          sscanf(optarg,"%d",&ScaleNUM);
          break;

        case 'o':
          sscanf(optarg,"%d",&OrientationNUM);
          break;

        case 'u':
          sscanf(optarg,"%f,%f",&Ul,&Uh);
          break;
        
        case 'b':
          sscanf(optarg,"%d,%d",&Sub_X,&Sub_Z);
          break;

        case 't':
          sscanf(optarg,"%d",&method);
          break;
          
        case 'x':
          sscanf(optarg, "%d", &EndScale);
          break;
         
        case 'i':
          sscanf(optarg,"%d",&foregroundThre);
          break;
          
        case 'I':
          sscanf(optarg,"%s",inputMaskName);
          if (strcmp(inputMaskName,"NULL")!=0)
            inputMaskOrNot=YYES;
          break;
          
        case 'S':
          SmoothOrNot=YYES;
          break;
          
        case 'F':
          DCOrNot=YYES;
          break;          

        case 'e':
          CannyGuassianSigma=atof(optarg);
          break;
                  
        case 'P':
          UsePolarCoordinate = YYES ;
          break;

        case 'D':
          DebugOrNot=YYES;
          break;
          
        case 'R':
          regionGrowingOrNot=YYES;
          break;
          
        case 'V':
          fillHoleOrNot=YYES;
          break;
         
        case 'M':
          saveMaskOrNot=YYES;
          break;
          
        default:
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

   sprintf(ResultDirectory,"%s",argv[1]);

   StartScale=0;
   EndScale=MIN(ScaleNUM-1, EndScale-1);
   StartOrient=0;
   EndOrient=OrientationNUM-1;
   
   times(&startTime);

      
   // open image
   printf("\nReading input image...\n");
   Image* image = ReadNiftiImage(argv[0], DT_UNSIGNED_CHAR);
   if (image == NULL) {
        cerr << "Failed to read image from file " << argv[0] << endl;
        exit(1);
   }
   UC_img3d = image->img.uc;

   Orig_X = image->region.nx;
   Orig_Y = image->region.ny;
   Orig_Z = image->region.nz;
   imageSize3d.x = Orig_X;
   imageSize3d.y = Orig_Y;
   imageSize3d.z = Orig_Z;

   res_X  = image->hdr.pixdim[1];
   res_Y  = image->hdr.pixdim[2];
   res_Z  = image->hdr.pixdim[3];

   float minRes;
   minRes = MIN(res_X,res_Y);
   minRes = MIN(res_Z,minRes);
   float resRatio_X = res_X/minRes;
   float resRatio_Y = res_Y/minRes;
   float resRatio_Z = res_Z/minRes;

   printf("image size: (%d,%d)\n", Orig_X, Orig_Y) ;
   printf("resolution: (%f %f %f)\n",res_X, res_Y,res_Z);
   printf("resolution: min=%f, ratio = (%f,%f,%f)\n", minRes, resRatio_X, resRatio_Y, resRatio_Z);
   printf("smooth axis: %dx%d \n",SmoothLongAxis, SmoothShortAxis);
   printf("side=%d\n", side) ;
   printf("ScaleNum=%d OrientationNum=%d\n",ScaleNUM,OrientationNUM);
   printf("method to convert Gabor responses into attributes = %d (1 for translate; 2 for linear scaling)\n", method);
   printf("foreground threshod for generating mask = %d\n", foregroundThre);
   printf("when generating mask, fill hole or not = %d (0 for No and 1 for Yes)\n", fillHoleOrNot);
   printf("when generating mask, region growing or not = %d (0 for No and 1 for Yes)\n", regionGrowingOrNot);
   printf("Ul=%f Uh=%f\n",Ul,Uh);
   printf("DCOrNot=%d\n", DCOrNot) ;
   printf("SmoothOrNot=%d\n",SmoothOrNot);
   printf("Sub_X=%d Sub_Z=%d\n",Sub_X,Sub_Z);   
   printf("CannyGuassianSigma=%f\n", CannyGuassianSigma) ;
   printf("ResultDirectory=%s\n",ResultDirectory);
   printf("DebugOrNot=%d\n",DebugOrNot);
   printf("slicesnum=%d\n", Orig_Z);

   featureImage=(unsigned char****)malloc(sizeof(unsigned char***)*ScaleNUM*OrientationNUM*2);
   featureImageFloat=(float****)malloc(sizeof(float***)*ScaleNUM*OrientationNUM*2);
   #pragma omp parallel for private(s,n) num_threads(5)
   for (s=0;s<ScaleNUM;s++)
     for (n=0;n<OrientationNUM;n++)
       {
        featureImage[(s*OrientationNUM+n)*2]=UCalloc3d(Orig_X,Orig_Y,Orig_Z);   
        featureImage[(s*OrientationNUM+n)*2+1]=UCalloc3d(Orig_X,Orig_Y,Orig_Z);   
        featureImageFloat[(s*OrientationNUM+n)*2]=Falloc3d(Orig_X,Orig_Y,Orig_Z);  
        featureImageFloat[(s*OrientationNUM+n)*2+1]=Falloc3d(Orig_X,Orig_Y,Orig_Z);   
       }
 
    Image* mask = NULL;
    if (inputMaskOrNot) {
        mask = ReadNiftiImage(inputMaskName);
        if (mask == NULL) {
            cerr << "Failed to read mask image from file " << inputMaskName << endl;
            // TODO free already allocated memory...
            exit(1);
        }
        mask3d = mask->img.uc;
    } else {
        mask = new Image(image->region.nx, image->region.ny, image->region.nz,
                         image->hdr.pixdim[1], image->hdr.pixdim[2], image->hdr.pixdim[3],
                         DT_UNSIGNED_CHAR, 1, Image::FORMAT_DRAMMS);
        if (mask == NULL) {
            fprintf(stderr, "Failed to allocate memory!\n");
            // TODO free already allocated memory...
            exit(1);
        }
        printf("\nGenerate mask for feature extraction...\n");
        printf("----------------------------\n");
        mask3d = mask->img.uc;
        generateMask(UC_img3d, imageSize3d, mask3d, foregroundThre, fillHoleOrNot, regionGrowingOrNot);
        printf("----------------------------\n");
    }

   SubSampled_X=Orig_X/Sub_X;
   SubSampled_Y=Orig_Y/Sub_X;
   SubSampled_Z=Orig_Z/Sub_Z;
      
   printf("\nSubsampled Volumn Image:(%d %d %d)\n",SubSampled_X,SubSampled_Y,SubSampled_Z);
   
   Fea_img3d=UCalloc3d(SubSampled_X,SubSampled_Y,SubSampled_Z);      
   
   //Detect Probe Center   
   UC_img2d = UCalloc2d(SubSampled_X,SubSampled_Y);
   mask2d = UCalloc2d(SubSampled_X,SubSampled_Y);
   k=Orig_Z/2;
   for (i=0;i<SubSampled_X;i++)
     for (j=0;j<SubSampled_Y;j++)
       {
       UC_img2d[i][j]=UC_img3d[k][i*Sub_X][j*Sub_X];
       mask2d[i][j] = mask3d[k][i*Sub_X][j*Sub_X];
       }
   
   if (DetectOrNot)
     DetectCircle(UC_img2d,&CenterOfProbe,&Radius);

   
   if (method==2)  
     {
     minImaginaryHori=(float***)malloc(sizeof(float**)*SubSampled_Z);
     maxImaginaryHori=(float***)malloc(sizeof(float**)*SubSampled_Z);
     minImaginaryVert=(float***)malloc(sizeof(float**)*SubSampled_Y);
     maxImaginaryVert=(float***)malloc(sizeof(float**)*SubSampled_Y);
     minRealHori=(float***)malloc(sizeof(float**)*SubSampled_Z);
     maxRealHori=(float***)malloc(sizeof(float**)*SubSampled_Z);
     minRealVert=(float***)malloc(sizeof(float**)*SubSampled_Y);
     maxRealVert=(float***)malloc(sizeof(float**)*SubSampled_Y);
     for (i=0;i<SubSampled_Z;i++)
       {
       minImaginaryHori[i]=Falloc2d(ScaleNUM,OrientationNUM);
       maxImaginaryHori[i]=Falloc2d(ScaleNUM,OrientationNUM);
       minRealHori[i]=Falloc2d(ScaleNUM,OrientationNUM);
       maxRealHori[i]=Falloc2d(ScaleNUM,OrientationNUM);
       }
     for (i=0;i<SubSampled_Y;i++)
       {
       minImaginaryVert[i]=Falloc2d(ScaleNUM,OrientationNUM);
       maxImaginaryVert[i]=Falloc2d(ScaleNUM,OrientationNUM);
       minRealVert[i]=Falloc2d(ScaleNUM,OrientationNUM);
       maxRealVert[i]=Falloc2d(ScaleNUM,OrientationNUM);
       }
     }
     
     
     
     
   //Calculate Gabor feature along Z-axis
   
   HoriOrVert=AlongZAxis;
   
   printf("\nCalculating Horizontal Slice (along z axis)...\n\n");
   printf("StartScale=%d EndScale=%d StartOrient=%d EndOrient=%d\n",StartScale,EndScale,StartOrient,EndOrient);

   //Calculating Gabor features
   for (k=0;k<SubSampled_Z;k++)   
     {
       for (i=0;i<SubSampled_X;i++)
         for (j=0;j<SubSampled_Y;j++) 
           {         
           UC_img2d[i][j]=UC_img3d[k*Sub_Z][i*Sub_X][j*Sub_X];
           mask2d[i][j]=mask3d[k*Sub_Z][i*Sub_X][j*Sub_X];
           }

       //printf("\n\nCalculating Gabor Features of %d Slice...\n",k*Sub_Z);
       
       //For prostate,we need rotation invariant
       //CalculateGaborFeatureFromOneSlice(UC_img2d,SubSampled_X,SubSampled_Y,YYES,k);

       //For brain, we needn't rotation invariant
       CalculateGaborFeatureFromOneSlice(UC_img2d,mask2d,SubSampled_X,SubSampled_Y,res_X, res_Y, NNO,k, method);
     }
   
   
   //Save result         
   for (s=StartScale;s<=EndScale;s++)
     for (n=StartOrient;n<=EndOrient;n++)
       {
         if (method==2) { // linear scaling of gabor response into attribute values
            linearScalingGaborResponseIntoAttributeValues(s,n,ScaleNUM, OrientationNUM, SubSampled_Z, Orig_X, Orig_Y, Orig_Z, mask3d);
         }

         sprintf(filename,"%s_3dHori_F_imag.%d_%d",ResultDirectory,s,n);         
         printf("Save feature image to %s...\n",filename);
         WriteUCImage(filename, image, featureImage[(s*OrientationNUM+n)*2]);

         sprintf(filename,"%s_3dHori_F_real.%d_%d",ResultDirectory,s,n);
         printf("Save feature image to %s...\n",filename);
         WriteUCImage(filename, image, featureImage[(s*OrientationNUM+n)*2 + 1]);
       }

   UCfree2d(UC_img2d,SubSampled_X);
   UCfree2d(mask2d, SubSampled_X);


   if (SubSampled_Z>1)  // extracting vertical features is meaningful only when there is more than 1 slice in z direction
   {
   
   // changed on 10/27/2011, calculate gabor features along X or Y axis, instead of always along Y axis in the original version. Whether along X or Y axis depends on which axis has less slices.
   if (SubSampled_Y<=SubSampled_X)  // Calculate Gabor features along Y axis if there are less slices along Y axis
   { // if calculating features along Y axis
   HoriOrVert = AlongYAxis;

   UC_img2d = UCalloc2d(SubSampled_Z,SubSampled_X);
   mask2d = UCalloc2d(SubSampled_Z,SubSampled_X);
   
   printf("\nCalculating Vertical Slices (along Y axis)...\n\n");
   printf("StartScale=%d EndScale=%d StartOrient=%d EndOrient=%d\n",StartScale,EndScale,StartOrient,EndOrient);

   
   for (j=0;j<SubSampled_Y;j++)
     {
       for (k=0;k<SubSampled_Z;k++)
         for (i=0;i<SubSampled_X;i++)
           {
           UC_img2d[k][i]=UC_img3d[k*Sub_Z][i*Sub_X][j*Sub_X];
           mask2d[k][i] = mask3d[k*Sub_Z][i*Sub_X][j*Sub_X];
           }
              
       //printf("\n\nCalculating Gabor Features of %d Slice...\n",j*Sub_X);
       CalculateGaborFeatureFromOneSlice(UC_img2d,mask2d,SubSampled_Z,SubSampled_X, res_Z, res_X, NNO,j, method);
       
     }   
   
   
   for (s=StartScale;s<=EndScale;s++)
     for (n=StartOrient;n<=EndOrient;n++)
       {
         if (method==2) // linear scaling of gabor response into attribute values
            linearScalingGaborResponseIntoAttributeValues(s,n,ScaleNUM, OrientationNUM, SubSampled_Y, Orig_X, Orig_Y, Orig_Z, mask3d);

         // imaginary component
         sprintf(filename,"%s_3dVert_F_imag.%d_%d",ResultDirectory,s,n);         
         printf("Save feature image to %s...\n",filename);
         WriteUCImage(filename, image, featureImage[(s*OrientationNUM+n)*2]);

         // real component
         sprintf(filename,"%s_3dVert_F_real.%d_%d",ResultDirectory,s,n);
         printf("Save feature image to %s...\n",filename);
         WriteUCImage(filename, image, featureImage[(s*OrientationNUM+n)*2 + 1]);
       }
    UCfree2d(UC_img2d,SubSampled_Z);
    UCfree2d(mask2d, SubSampled_Z);
    } // if calculating features along Y axis
   else  // Calculate Gabor features along X axis if there are more slices along X axis
   { // if calculating features along X axis
   HoriOrVert = AlongXAxis;

   UC_img2d = UCalloc2d(SubSampled_Z,SubSampled_Y);
   mask2d = UCalloc2d(SubSampled_Z,SubSampled_Y);
   
   printf("\nCalculating Vertical Slices (along x axis)...\n\n");
   printf("StartScale=%d EndScale=%d StartOrient=%d EndOrient=%d\n",StartScale,EndScale,StartOrient,EndOrient);

   for (i=0;i<SubSampled_X;i++)
     {
       for (k=0;k<SubSampled_Z;k++)
         for (j=0;j<SubSampled_Y;j++)
           {
           UC_img2d[k][j]=UC_img3d[k*Sub_Z][i*Sub_Y][j*Sub_Y];
           mask2d[k][j] = mask3d[k*Sub_Z][i*Sub_Y][j*Sub_Y];
           }
              
       //printf("\n\nCalculating Gabor Features of %d Slice...\n",i*Sub_Y);
       CalculateGaborFeatureFromOneSlice(UC_img2d,mask2d,SubSampled_Z,SubSampled_Y, res_Z, res_Y, NNO,i, method);
       
     }   
   
   
   for (s=StartScale;s<=EndScale;s++)
     for (n=StartOrient;n<=EndOrient;n++)
       {
         if (method==2) // linear scaling of gabor response into attribute values
            linearScalingGaborResponseIntoAttributeValues(s,n,ScaleNUM, OrientationNUM, SubSampled_X, Orig_X, Orig_Y, Orig_Z, mask3d);
            
         // imaginary component
         sprintf(filename,"%s_3dVert_F_imag.%d_%d",ResultDirectory,s,n);         
         printf("Save feature image to %s...\n",filename);
         WriteUCImage(filename, image, featureImage[(s*OrientationNUM+n)*2]);

         // real component
         sprintf(filename,"%s_3dVert_F_real.%d_%d",ResultDirectory,s,n);
         printf("Save feature image to %s...\n",filename);
         WriteUCImage(filename, image, featureImage[(s*OrientationNUM+n)*2 + 1]);
       }
    UCfree2d(UC_img2d,SubSampled_Z);
    UCfree2d(mask2d, SubSampled_Z);
    } // if calculating features along X axis
    
    } // if (SubSampled_Z>1)
    
    
   
   if (saveMaskOrNot)
    {
    sprintf(filename,"%s_mask",ResultDirectory);
    printf("save mask to %s\n", filename);
    WriteUCImage(filename, image, mask3d);
    }

   UCfree3d(Fea_img3d,SubSampled_Z,SubSampled_X);

   for (s=0;s<ScaleNUM;s++)
     for (n=0;n<OrientationNUM;n++)
       {
     UCfree3d(featureImage[(s*OrientationNUM+n)*2],Orig_Z,Orig_X);
     UCfree3d(featureImage[(s*OrientationNUM+n)*2+1],Orig_Z,Orig_X);
     Ffree3d(featureImageFloat[(s*OrientationNUM+n)*2],Orig_Z,Orig_X);
     Ffree3d(featureImageFloat[(s*OrientationNUM+n)*2+1],Orig_Z,Orig_X);
       }
   
   free(featureImage);
   free(featureImageFloat);
   if (image) delete image;
   if (mask) delete mask;

   times(&endTime);
   durTime=((double)endTime.tms_utime-(double)startTime.tms_utime);
   printf("\n\nDuration=%.2f seconds\n",durTime/100.);

   return 0;   
}

void DetectCircle(unsigned char** i_image,Fvector2d* center,double* radius)
{
  int i,j,i1,i2;  
  unsigned char** subimage,**edgeimage;
  int subX,subY;
  Fvector2d edgepoint[10000];
  int edgepointnum;
  int binthre,edgethre;
  double centerX,centerY;
  double pre,var;
  Matrix* temp;
  char filename[500];
      
  printf("Detect circle center...\n");

  //Get subimage
  subX=static_cast<int>(SubSampled_X*.75);
  subY=static_cast<int>(SubSampled_Y*.5);
  subimage=UCalloc2d(subX,subY);

  if (DebugOrNot)
    printf("SubX=%d SubY=%d\n",subX,subY);
      
  for (i=0;i<subX;i++)
    for (j=0;j<subY;j++)
      subimage[i][j]=i_image[i+SubSampled_X/2-subX/2][j];
    
  //Binarize subimage
  binthre=0;
  
  for (i=0;i<subX;i++)
    for (j=0;j<subY;j++)
      {
        if (subimage[i][j]>binthre)
          subimage[i][j]=255;
        else
          subimage[i][j]=0;
      }
  
  //Detect edge point
  edgethre=100;
  
  edgeimage=UCalloc2d(subX,subY);
  for (i=1;i<subX-1;i++)    
    for(j=1;j<subY-1;j++)              
      if (abs(subimage[i][j-1]-subimage[i][j+1])>edgethre||abs(subimage[i-1][j]-subimage[i+1][j])>edgethre)
        edgeimage[i][j]=255;
      else
        edgeimage[i][j]=0;      
  
  //Detect coordinate of probe center
  
  //Get rid of noise edge point
  edgepointnum=0;  
  pre=INFINITE;
  
  for (j=0;j<subY;j++)
    {
      i1=subX/2;

      while ((i1>=0)&&(edgeimage[i1][j]!=255))
        i1--;

      i2=subX/2+1;
      while ((i2<subX)&&(edgeimage[i2][j]!=255))
        i2++;

      if (i1>0&&i2<subX)
        {
          if (i2-i1<pre*1.2)
            {
              edgepoint[edgepointnum].x=i1;
              edgepoint[edgepointnum].y=j;
              edgepoint[edgepointnum+1].x=i2;
              edgepoint[edgepointnum+1].y=j;
              
              edgepointnum+=2;

              pre=i2-i1;
            }
          else
            break;
        }                        
    }
  
  //Calculate CenterX
  CreateMatrix(&temp,edgepointnum,1);
  for (i=0;i<edgepointnum;i++)    
    temp->data[i][0]=edgepoint[i].x;
  
  Mat_Mean(&centerX,temp);
  
  j=(int)edgepoint[edgepointnum-1].y;
    
  //Calculate CenterY
  pre=INFINITE;
  while (j>0)
    {     
      for (i=0;i<edgepointnum;i++)        
        temp->data[i][0]=sqrt(pow(static_cast<double>(centerX - edgepoint[i].x), 2.0) + pow(static_cast<double>(j - edgepoint[i].y), 2.0));
      
      Mat_Variance(&var, temp);

      if (var<pre)        
        {
          centerY=j;
          pre=var;
        }        
      j--;
    }

  if (pre>1)
    {
      printf("Var=%f Something strange happened in \"DetectCircle\"!\n",pre);
      return;
    }
          
    
  for (i=0;i<edgepointnum;i++)        
    temp->data[i][0]=sqrt(pow((centerX-edgepoint[i].x),2.0)+pow((centerY-edgepoint[i].y),2.0));

  Mat_Mean(radius,temp);

  center->x=centerX+SubSampled_X/2-subX/2;
  center->y=centerY;
    
  if (DebugOrNot)
    {
      printf("Probe:X(%.3f) Y:(%.3f) R:(%.3f)\n",center->x,center->y,*radius);

      for (i=-1;i<=1;i++)
        for (j=-1;j<=1;j++)      
          edgeimage[(int)centerX+i][(int)centerY+j]=255;              
      
      sprintf(filename,"%s_ProbeCenter.raw",ResultDirectory);
      WriteImg2D(filename, edgeimage, subX, subY) ;
    }
    
  UCfree2d(subimage,subX);
  UCfree2d(edgeimage,subX);
  FreeMatrix(temp);
         
}

void CalculateGaborFeatureFromOneSlice(unsigned char** i_image,unsigned char **mask, int image_X,int image_Y, float res_X, float res_Y, int RotateInvariantOrNot,int SliceIndex, int method)
{
  char filename[500];
  int FeatureNum;
  double **Features_Real, **Features_Imginary ;
  int i;
  
  unsigned char** Img_edge;

  if ( (method!=1)&&(method!=2) )
    method=1; // default
    
  
  if (SmoothOrNot)
    {
      SmoothingImgByRemovingHighFrequencyFourierParemters(i_image, image_X, image_Y) ;
      
      if (DebugOrNot)
        {
          sprintf(filename,"%s_SmoothedByFourier.raw",ResultDirectory);
          WriteImg2D(filename, i_image, image_X, image_Y);
        }
      
      if (RotateInvariantOrNot)
        {
          SmoothingFilterOnTheOrientedEllipseDomain( SmoothLongAxis, SmoothShortAxis, i_image, image_X, image_Y );
          
          if (DebugOrNot)
            {
              sprintf(filename,"%s_SmoothedImg.raw",ResultDirectory);
              WriteImg2D(filename, i_image, image_X, image_Y);
            }
        }
    }
  
  // allocate the memory for texture feature
  FeatureNum = ScaleNUM*OrientationNUM;
  Features_Real     = DDalloc2d(image_X*image_Y, FeatureNum) ;
  Features_Imginary = DDalloc2d(image_X*image_Y, FeatureNum) ;
  WeightingOnFrequency = Dalloc1d(ScaleNUM) ;
  
  // setting weightings
  for(i=0; i<ScaleNUM; i++)
    {
      //WeightingOnFrequency[i] = exp(-pow((double)i/(double)ScaleNUM,2.0)/(2.0*4.0));      
      WeightingOnFrequency[i]=1.0;
    }
  GetGaborFeatureOfOneImageNew(Features_Real, Features_Imginary,i_image,mask,image_X,image_Y,res_X,res_Y, Ul,Uh,StartScale,EndScale,StartOrient,EndOrient,ScaleNUM, OrientationNUM,DCOrNot);
  
  if (RotateInvariantOrNot)
    CalculateGaborRotationInvariants_ToProbeCenter(Features_Real, Features_Imginary, image_X,image_Y,ScaleNUM,OrientationNUM);

  Put2DGaborResponseInto3DArray(Features_Real, Features_Imginary, mask, image_X, image_Y, StartScale,EndScale,StartOrient,EndOrient,ScaleNUM,OrientationNUM,SliceIndex, method);

  /*
  Img_edge=UCalloc2d(image_X,image_Y);
  DetectEdgeFromGaborFeatureVectors(Img_edge, image_X, image_Y, FeatureNum, Features_Real,Features_Imginary,CannyGuassianSigma);

  WriteImg2D("Edge.raw", Img_edge, image_X, image_Y);
  UCfree2d(Img_edge,image_X);
  */
  Dfree2d(Features_Real,image_X*image_Y);
  Dfree2d(Features_Imginary,image_X*image_Y);
  free(WeightingOnFrequency);
  
}

/* subprograms ... */
float gaussian(float x, float s)
{
  return(exp((-x*x)/(2*s*s)));
}

double hypotenuse(double x,double y)
{
    if (x==0.0 && y==0.0) return(0.0);
    else return(hypot(x,y));
}

void SwapMatrix(Matrix* A)
{
  Matrix* Temp;
  int h,w,xs,ys;

  CreateMatrix(&Temp,A->height,A->width);

  xs=A->height;
  ys=A->width;

  for (h=0;h<xs/2;h++)
    for (w=0;w<ys/2;w++)
      {
        Temp->data[h][w]=A->data[h+xs/2][w+ys/2];
        Temp->data[h+xs/2][w+ys/2]=A->data[h][w];
      }
  
  for (h=0;h<xs/2;h++)
    for (w=ys/2;w<ys;w++)
      {
        Temp->data[h][w]=A->data[h+xs/2][w-ys/2];
        Temp->data[h+xs/2][w-ys/2]=A->data[h][w];
      }

  for (h=0;h<xs;h++)
    for (w=0;w<ys;w++)
      A->data[h][w]=Temp->data[h][w];
  
  FreeMatrix(Temp);
}


void SmoothingImgByRemovingHighFrequencyFourierParemters(unsigned char **UC_img, int image_X, int image_Y)
{
  int h, w, xs, ys;
  int hei, wid;
  Matrix *IMG, *IMG_imag,*F_real, *F_imag;
  float weight ;
  
  printf("Smoothing image using Fourier transformation ... \n\n") ;
  
  hei = image_X ;
  wid = image_Y ;
  
  xs = (int) pow(2.0, ceil(log2((double) image_X)));
  ys = (int) pow(2.0, ceil(log2((double) image_Y)));
  
  if (DebugOrNot)
    printf("(xs,ys)=(%d,%d)\n", xs, ys) ;
  
  CreateMatrix(&IMG, xs, ys);
  for(h=0; h<hei; h++) 
    for(w=0; w<wid; w++)
      IMG->data[h][w] = (double)UC_img[h][w] ;   


  CreateMatrix(&F_real, xs, ys);
  CreateMatrix(&F_imag, xs, ys);
  CreateMatrix(&IMG_imag, xs, ys);
  
  Mat_FFT2(F_real, F_imag, IMG, IMG_imag);

  //Swap matrix, shift the DC conefficient to the center of the matrix
  SwapMatrix(F_real);
  SwapMatrix(F_imag);
  
  /* directly setting high frequency to zeros */
  for(h=0;h<xs;h++)
    for(w=0;w<ys;w++) 
      {
        //weight = exp(-(double)(h*h+w*w)/(2.*(double)(xs*xs+ys*ys)/1.0) ) ; /* 4.0 */

        weight = exp(-((h-xs/2)*(h-xs/2)+(w-ys/2)*(w-ys/2))/(2.*(xs/2*xs/2+ys/2*ys/2)/4.0));
        F_real->data[h][w]*= weight;
        F_imag->data[h][w]*= weight;
      }

  //Swap matrix back
  SwapMatrix(F_real);
  SwapMatrix(F_imag);

  Mat_IFFT2(IMG, IMG_imag, F_real, F_imag);
      
  for(h=0; h<hei; h++)
    for(w=0; w<wid; w++)
      UC_img[h][w] = (int)(IMG_imag->data[h][w]);  

  /* save smoothing result */
  for(h=0; h<hei; h++)
    for(w=0; w<wid; w++)
      UC_img[h][w] = (int)(IMG->data[h][w]) ;
  
  /* free */
  FreeMatrix(F_real);
  FreeMatrix(F_imag);
  FreeMatrix(IMG);
  FreeMatrix(IMG_imag);
}


// ShortAxis is along the direction from the probe center, LongAxis is orthogonal to it.
/* Oct 10, 2001 */
#define  MaxIntensity 256 /* 255+1 */
void SmoothingFilterOnTheOrientedEllipseDomain( int SmoothLongAxis, int SmoothShortAxis, unsigned char **UC_img, int image_X, int image_Y )
{
  int i, j, d0, d1, c ;
  unsigned char **Smth_Img ;
  Fvector2d  *AxisOfEllipseWindow, CurrentP ;
  float mag ;
  int Count[MaxIntensity], TotalPixels, HalfNum, max, BestValue, average, LargestIntensity ;
  float IncreaseRate ;

  printf("Smoothing image using oriented Median filter ... \n\n") ;

  /* apply for space */
  AxisOfEllipseWindow = Fvector2dalloc1d(NumOfDimension); /* 2D */
  Smth_Img = UCalloc2d(image_X,image_Y);

  /* process */
  for(i=0; i<image_X; i++) 
    for(j=0; j<image_Y; j++)
      {
        /*printf("i=%d j=%d\n", i,j) ;*/
        if( i!=CenterOfProbe.x || j!=CenterOfProbe.y )
          {
            /* get the axis of the filter window */
            /* first */
            AxisOfEllipseWindow[0].x = i-CenterOfProbe.x ;
            AxisOfEllipseWindow[0].y = j-CenterOfProbe.y ;
            /* normalization ... */
            mag = sqrt(AxisOfEllipseWindow[0].x*AxisOfEllipseWindow[0].x + AxisOfEllipseWindow[0].y*AxisOfEllipseWindow[0].y) ;
            AxisOfEllipseWindow[0].x /= mag ;
            AxisOfEllipseWindow[0].y /= mag ; 
            /* second */
            AxisOfEllipseWindow[1].x =  AxisOfEllipseWindow[0].y ;
            AxisOfEllipseWindow[1].y = -AxisOfEllipseWindow[0].x ;
            
            IncreaseRate = 1.0+mag/image_X ;
            
            /* reset */
            for(c=0; c<MaxIntensity; c++)    Count[c] = 0 ;
            TotalPixels = 0;
            
            for(d0=-SmoothShortAxis; d0<=SmoothShortAxis; d0++) 
              for(d1=static_cast<int>(-SmoothLongAxis*IncreaseRate);  d1<=static_cast<int>(SmoothLongAxis*IncreaseRate);  d1++) 
                {
                  CurrentP.x = d0*AxisOfEllipseWindow[0].x + d1*AxisOfEllipseWindow[1].x + i + 0.5 ;
                  CurrentP.y = d0*AxisOfEllipseWindow[0].y + d1*AxisOfEllipseWindow[1].y + j + 0.5;
                  
                  if( CurrentP.x<image_X && CurrentP.y<image_Y &&  CurrentP.x>=0 && CurrentP.y>=0 )
                    {
                      TotalPixels ++ ;
                      Count[ UC_img[(int)CurrentP.x][(int)CurrentP.y] ] ++ ;
                    }
                }
            
            /* get the median one */
            /*printf("TotalPixels=%d\n", TotalPixels) ;*/
            HalfNum = 0 ;
            max = 0 ;   BestValue=-1 ;
            average = 0;
            LargestIntensity = 0 ;
            for(c=0; c<MaxIntensity; c++)    
              {
                /* statistically majority greyvalue */
                if( Count[c]>max ) { max = Count[c] ; BestValue=c; }
                
                /* average greyvalue */
                average += Count[c]*c ;
                
                /* largest intensity */
                if( Count[c]!=0 && c>LargestIntensity)  LargestIntensity = c ;  

                /* middle position of histogram  */
                HalfNum += Count[c] ;
                if( HalfNum>=TotalPixels/2 )  
                  {
                    Smth_Img[i][j] = c ;
                    c = MaxIntensity ; /* for setting STOP sign*/
                  }
              }
            /*Smth_Img[i][j] = BestValue ;*/
            /*Smth_Img[i][j] = LargestIntensity ;*/
            
            /*Smth_Img[i][j] = average/TotalPixels ;*/
            /*printf("average=%d TotalPixels=%d average/TotalPixels=%d\n", average, TotalPixels, average/TotalPixels) ;*/
          }
      }
  
  /* back tt */
  for(i=0; i<image_X; i++) 
    for(j=0; j<image_Y; j++)
      UC_img[i][j] = Smth_Img[i][j] ;

  /* free */
  UCfree2d(Smth_Img, image_X) ;
  free(AxisOfEllipseWindow) ;
}

/* extract texture feature */
void GetGaborFeatureOfOneImage(double **Features_Real, double **Features_Imginary, unsigned char **UC_img, int image_X, int image_Y,int side,double Ul,double Uh, int ScaleNUM, int OrientationNUM, int flag)
{
  int h, w, xs, ys, hei, wid, s, n, base;
  Matrix *IMG, *IMG_imag, *Gr, *Gi, *Tmp_1, *Tmp_2, *F_1, *F_2, *G_real, *G_imag, *F_real, *F_imag, *F;
  double m, v;
  
  hei = image_X ;
  wid = image_Y ;
  
  xs = (int) pow(2.0, ceil(log2((double) image_X)));
  ys = (int) pow(2.0, ceil(log2((double) image_Y)));
  
  CreateMatrix(&IMG, xs, ys);
  for(h=0; h<hei; h++) 
    for(w=0; w<wid; w++)
      IMG->data[h][w] = (double)UC_img[h][w] ;
  
  CreateMatrix(&F_real, xs, ys);
  CreateMatrix(&F_imag, xs, ys);
  CreateMatrix(&IMG_imag, xs, ys);
  
  Mat_FFT2(F_real, F_imag, IMG, IMG_imag);
  
  
  /* ----------- compute the Gabor filtered output ------------- */  
  CreateMatrix(&Gr, 2*side+1, 2*side+1);
  CreateMatrix(&Gi, 2*side+1, 2*side+1);
  CreateMatrix(&Tmp_1, xs, ys);
  CreateMatrix(&Tmp_2, xs, ys);
  CreateMatrix(&F_1, xs, ys);
  CreateMatrix(&F_2, xs, ys);
  CreateMatrix(&G_real, xs, ys);
  CreateMatrix(&G_imag, xs, ys);
  CreateMatrix(&F, hei, wid);
  
  base = ScaleNUM*OrientationNUM;
  
  for(s=0; s<ScaleNUM; s++) 
    {
      printf("s=%d\n", s) ;      
      for(n=0; n<OrientationNUM; n++) 
        {
          printf("n=%d\n", n);
          
          //GaborFunction(Gr, Gi, s+1, n+1, Ul, Uh, ScaleNUM, OrientationNUM,flag); 
          GaborFunction(Gr, Gi, s+1, n+1, Ul, Uh, ScaleNUM, OrientationNUM,flag); 

          Mat_Copy(F_1, Gr, 0, 0, 0, 0, 2*side, 2*side);
          Mat_Copy(F_2, Gi, 0, 0, 0, 0, 2*side, 2*side);


          for (h=0;h<5;h++)
            for (w=0;w<5;w++)              
              printf("(%d %d):%f\n",h,w,F_1->data[h][w]);          

          Mat_FFT2(G_real, G_imag, F_1, F_2);
          
          for (h=0;h<5;h++)
            for (w=0;w<5;w++)              
              printf("Real(%d %d):%f\n",h,w,G_real->data[h][w]);

          printf("\n\n");
          for (h=0;h<5;h++)
            for (w=0;w<5;w++)              
              printf("Imag(%d %d):%f\n",h,w,G_real->data[h][w]);          
          
          Mat_Product(Tmp_1, G_real, F_real);
          Mat_Product(Tmp_2, G_imag, F_imag);
          Mat_Substract(IMG, Tmp_1, Tmp_2);
          
          Mat_Product(Tmp_1, G_real, F_imag);
          Mat_Product(Tmp_2, G_imag, F_real);
          Mat_Sum(IMG_imag, Tmp_1, Tmp_2);
          
          Mat_IFFT2(Tmp_1, Tmp_2, IMG, IMG_imag);
          Mat_Shift(IMG, Tmp_1, side);
          Mat_Shift(IMG_imag, Tmp_2, side);

          /* save features, for returning back */
          for(h=0;h<hei;h++)
            for(w=0;w<wid;w++)
              {
                Features_Real[h*wid+w][s*OrientationNUM+n] = IMG->data[h][w]*WeightingOnFrequency[s] ;
                Features_Imginary[h*wid+w][s*OrientationNUM+n] = IMG_imag->data[h][w]*WeightingOnFrequency[s] ; 
              }          
        }
    }
  
  
  /* free */
  FreeMatrix(Gr);
  FreeMatrix(Gi);
  FreeMatrix(Tmp_1);
  FreeMatrix(Tmp_2);
  FreeMatrix(F_1);
  FreeMatrix(F_2);
  FreeMatrix(G_real);
  FreeMatrix(G_imag);
  FreeMatrix(F_real);
  FreeMatrix(F_imag);
  FreeMatrix(IMG);
  FreeMatrix(IMG_imag);
  FreeMatrix(F);
}

void GetGaborFeatureOfOneImageNew(double **Features_Real, double **Features_Imginary, unsigned char **UC_img, unsigned char **mask, int image_X, int image_Y, float res_X, float res_Y, double Ul,double Uh, int StartScale,int EndScale,int StartOrient,int EndOrient,int ScaleNUM, int OrientationNUM, int flag)
{
  Matrix *IMG_real,*IMG_imag,*F_real,*F_imag,*G_real,*G_imag,*Out_real,*Out_imag;
  int xs,ys;
  int h,w,hei,wid;
  int hh,ww;
  int s,n;
  float minResolution = MIN(res_X, res_Y);   
  float resolutionRatio_X = floor(res_X/minResolution*2.0+0.5)/2.0; 
  float resolutionRatio_Y = floor(res_Y/minResolution*2.0+0.5)/2.0; 
  
  
  if (resolutionRatio_X == resolutionRatio_Y)   // isotropic pixel
  { //if (resolutionRatio_X == resolutionRatio_Y)  
  hei=image_X;
  wid=image_Y;

  xs = (int) pow(2.0, ceil(log2((double) image_X)));
  ys = (int) pow(2.0, ceil(log2((double) image_Y)));
  
  
  CreateMatrix(&IMG_real, xs, ys);
  CreateMatrix(&IMG_imag, xs, ys);
  CreateMatrix(&F_real, xs, ys);
  CreateMatrix(&F_imag, xs, ys);
  CreateMatrix(&G_real, xs, ys);
  CreateMatrix(&G_imag, xs, ys);
  CreateMatrix(&Out_real, xs, ys);
  CreateMatrix(&Out_imag, xs, ys);
  
  for(h=0; h<hei; h++) 
    for(w=0; w<wid; w++)
     {
      IMG_real->data[h][w] = (double)UC_img[h][w] ;
      IMG_imag->data[h][w] = 0.0;
     }

  
  Mat_FFT2(F_real, F_imag, IMG_real, IMG_imag);

  SwapMatrix(F_real);
  SwapMatrix(F_imag);
  
  
  for(s=StartScale; s<=EndScale; s++)
    {
      //printf("\nScale:%d\n", s) ;      
      for(n=StartOrient; n<=EndOrient; n++)
        {
          //printf("Orientation=%d\n", n);

          GaborFrequencyFunction(G_real,G_imag,s,n,Ul,Uh,ScaleNUM,OrientationNUM,DCOrNot);

          Mat_Product(Out_real, G_real, F_real);
          Mat_Product(Out_imag, G_real, F_imag);

          SwapMatrix(Out_real);
          SwapMatrix(Out_imag);
          
          Mat_IFFT2(IMG_real, IMG_imag,Out_real, Out_imag);

          #pragma omp parallel for private(h,w,hh,ww) num_threads(5)
          /* save features, for returning back */
          for(h=0;h<hei;h++)
            for(w=0;w<wid;w++)
              {
              Features_Real[h*wid+w][s*OrientationNUM+n] = IMG_real->data[h][w]*WeightingOnFrequency[s] ;
              Features_Imginary[h*wid+w][s*OrientationNUM+n] = IMG_imag->data[h][w]*WeightingOnFrequency[s] ; 
              }
        }//for
    } //for
  }  //if (resolutionRatio_X == resolutionRatio_Y)  
  else  // anisotropic pixel
  {// anisotropic pixel
  float gcd = GetGreatestCommonDivisor(resolutionRatio_X, resolutionRatio_Y);
  int zoomFactor_X = (int)(resolutionRatio_X/gcd+0.49);
  int zoomFactor_Y = (int)(resolutionRatio_Y/gcd+0.49);
  
  hei=image_X*zoomFactor_X;
  wid=image_Y*zoomFactor_Y;

  xs = (int) pow(2.0, ceil(log2((double) hei)));
  ys = (int) pow(2.0, ceil(log2((double) wid)));
  
  CreateMatrix(&IMG_real, xs, ys);
  CreateMatrix(&IMG_imag, xs, ys);
  CreateMatrix(&F_real, xs, ys);
  CreateMatrix(&F_imag, xs, ys);
  CreateMatrix(&G_real, xs, ys);
  CreateMatrix(&G_imag, xs, ys);
  CreateMatrix(&Out_real, xs, ys);
  CreateMatrix(&Out_imag, xs, ys);
  
  unsigned char **zoomedUC_img;
  zoomedUC_img = UCalloc2d(hei,wid);
  zoomInImage(UC_img, zoomedUC_img, image_X, image_Y, hei, wid, zoomFactor_X, zoomFactor_Y);
  for(h=0; h<hei; h++) 
    for(w=0; w<wid; w++)
     {
      IMG_real->data[h][w] = (double)zoomedUC_img[h][w] ;
      IMG_imag->data[h][w] = 0.0;
     }
  
  Mat_FFT2(F_real, F_imag, IMG_real, IMG_imag);

  SwapMatrix(F_real);
  SwapMatrix(F_imag);
  
  for(s=StartScale; s<=EndScale; s++)
    {
      //printf("\nScale:%d\n", s) ;      
      for(n=StartOrient; n<=EndOrient; n++)
        {
          //printf("Orientation=%d\n", n);
                  
          GaborFrequencyFunction(G_real,G_imag,s,n,Ul,Uh,ScaleNUM,OrientationNUM,DCOrNot);

          Mat_Product(Out_real, G_real, F_real);
          Mat_Product(Out_imag, G_real, F_imag);

          SwapMatrix(Out_real);
          SwapMatrix(Out_imag);
          
          Mat_IFFT2(IMG_real, IMG_imag,Out_real, Out_imag);

          #pragma omp parallel for private(h,w,hh,ww) num_threads(5)
          /* save features, for returning back */
          for(h=0;h<image_X;h++)
            for(w=0;w<image_Y;w++)
              {
              //printf("thread %d \t",omp_get_thread_num());
              
              hh = (h+1)*zoomFactor_X-1;
              ww = (w+1)*zoomFactor_Y-1;
                
              Features_Real[h*image_Y+w][s*OrientationNUM+n] = IMG_real->data[hh][ww]*WeightingOnFrequency[s] ;
              Features_Imginary[h*image_Y+w][s*OrientationNUM+n] = IMG_imag->data[hh][ww]*WeightingOnFrequency[s] ; 
              }
        }//for
    } //for
  }// anisotropic pixel
  
  FreeMatrix(IMG_real);
  FreeMatrix(IMG_imag);
  FreeMatrix(F_real);
  FreeMatrix(F_imag);
  FreeMatrix(G_real);
  FreeMatrix(G_imag);
  FreeMatrix(Out_real);
  FreeMatrix(Out_imag);
  
}

//This function generate Gabor function in frequency domain
void GaborFrequencyFunction(Matrix *Gr, Matrix *Gi, int Scale,int Orientation, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag)
{
  double Uvar,Vvar,a;
  double tmp;
  double u,v,i1,j1;
  int i,j,height,width;  
  double dilation,theta;
  

  a=pow(Uh/Ul,1.0/(double)(ScaleNUM-1));
  dilation=pow(a,(double)Scale);
  
  
  Uvar=(a-1)*Uh/((a+1)*sqrt(2.0*log(2.0)));
  tmp=2.0*log(2.0)*Uvar/Uh;
  Vvar = tan(PI/(2.0*(OrientationNUM-1)))*(Uh-2.0*log(Uvar*Uvar/Uh))/sqrt(2.0*log(2.0)-tmp*tmp);
      
  height=Gr->height;
  width=Gr->width;

  
  theta=E_PI/(double)OrientationNUM*(double)Orientation;

  //printf("dilation=%f,theta=%f Uvar=%f Vvar=%f Uh=%f\n",dilation,theta,Uvar,Vvar,Uh);

  for (i=-height/2;i<height/2;i++)
    for (j=-width/2;j<width/2;j++)
      {
        i1=2*E_PI/(double)height*(double)i;
        j1=2*E_PI/(double)width*(double)j;
        
        u=(i1*cos(theta)+j1*sin(theta))*dilation;
        v=(-i1*sin(theta)+j1*cos(theta))*dilation;

        Gr->data[i+(int)(height/2)][j+(int)(width/2)]=GaborMotherFrequencyFunction(u,v,Uvar,Vvar,Uh);

        //if (abs(i)<3&&abs(j)<3)
        //  printf("(%.3f %.3f):%.3f\n",u,v,Gr->data[i+(int)(height/2)][j+(int)(width/2)]);
      }

  //if flag==NNO,remove DC coefficient
  if (flag==NNO)
    Gr->data[height/2][width/2]=0;
}

double GaborMotherFrequencyFunction(double u,double v,double Uvar,double Vvar,double Shift)
{
  double result;

  result=exp(-0.5*(pow((u-Shift)/Uvar,2.0)+pow(v/Vvar,2.0)));

  return result;
}

void GaborSpatialFunction(Matrix *Gr, Matrix *Gi, int Scale,int Orientation, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag)
{
  double Xvar,Yvar;
  double Uvar,Vvar,a;
  double tmp;
  double x,y;
  int i,j,height,width;  
  double dilation,theta;
  double real,imag;

  a=pow(Uh/Ul,1.0/(double)(ScaleNUM-1));
  
  Uvar=(a-1)*Uh/((a+1)*sqrt(2.0*log(2.0)));

  tmp=2.0*log(2.0)*Uvar/Uh;
  Vvar = tan(PI/(2*(OrientationNUM-1)))*(Uh-2.0*log(Uvar*Uvar/Uh))/sqrt(2.0*log(2.0)-tmp*tmp);
  //Vvar = tan(PI/(2*OrientationNUM))*(Uh-2.0*log(2.0)*Uvar*Uvar/Uh)/sqrt(2.0*log(2.0)-tmp*tmp);
  
  Xvar=1/(2*E_PI*Uvar);
  Yvar=1/(2*E_PI*Vvar);

  height=Gr->height;
  width=Gr->width;

  dilation=pow(a,-(double)Scale);
  theta=E_PI/(double)OrientationNUM*(double)Orientation;

  printf("dilation=%f,theta=%f Uvar=%f Vvar=%f Uh=%f\n",dilation,theta,Uvar,Vvar,Uh);

  for (i=-height/2;i<height/2;i++)
    for (j=-width/2;j<width/2;j++)
      {
        x=dilation*((double)i*cos(theta)+(double)j*sin(theta));
        y=dilation*(-(double)i*sin(theta)+(double)j*cos(theta));

        GaborMotherSpatialFunction(x,y,Xvar,Yvar,Uh,&real,&imag);
        Gr->data[i+(int)(height/2)][j+(int)(width/2)]=real;
        Gi->data[i+(int)(height/2)][j+(int)(width/2)]=imag;
        
      }

  //if flag==NNO,remove DC coefficient
  if (flag==NNO)
    {
      tmp=0;
      for (i=0;i<height;i++)
        for (j=0;j<width;j++)
          tmp += Gr->data[i][j];
    
      tmp/= height*width;
      
      for (i=0;i<height;i++)
        for (j=0;j<width;j++)
          Gr->data[i][j]-=tmp;
    }    
}

void GaborMotherSpatialFunction(double x,double y,double Xvar,double Yvar,double Shift,double* real,double* imag)
{
  double G;
  
  G=1.0/(2*E_PI*Xvar*Yvar)*exp(-0.5*(pow(x/Xvar,2.0)+pow(y/Yvar,2.0)));
  
  *real=G*cos(2*E_PI*Shift*x);
  *imag=G*sin(2*E_PI*Shift*x);
}


/* ------------------------------------------------------------------------------------------------------
The Gabor function generates a Gabor filter with the selected index 's' and 'n' (scale and orientation, 
respectively) from a Gabor filter bank. This filter bank is designed by giving the range of spatial 
frequency (Uh and Ul) and the total number of scales and orientations used to partition the spectrum. 

The returned filter is stored in 'Gr' (real part) and 'Gi' (imaginary part).
--------------------------------------------------------------------------------------------------------*/
void GaborFunction(Matrix *Gr, Matrix *Gi, int s, int n, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag)
{
  double base, a, u0, z, Uvar, Vvar, Xvar, Yvar, X, Y, G, t1, t2, m;
  int x, y, side;
  
  base = Uh/Ul;
  a = pow(base, 1.0/(double)(ScaleNUM-1));
  
  u0 = Uh/pow(a, (double) ScaleNUM-s);
  
  Uvar = (a-1.0)*u0/((a+1.0)*sqrt(2.0*log(2.0)));
  
  z = -2.0*log(2.0)*(Uvar*Uvar)/u0;
  Vvar = tan(PI/(2*(OrientationNUM-1)))*(u0+z)/sqrt(2.0*log(2.0)-z*z/(Uvar*Uvar));
  
  Xvar = 1.0/(2.0*PI*Uvar);
  Yvar = 1.0/(2.0*PI*Vvar);
  
  t1 = cos(PI/OrientationNUM*(n-1.0));
  t2 = sin(PI/OrientationNUM*(n-1.0));
  
  side = (int) (Gr->height-1)/2;
  
  for (x=0;x<2*side+1;x++) {
    for (y=0;y<2*side+1;y++) {
      X = (double) (x-side)*t1+ (double) (y-side)*t2;
      Y = (double) -(x-side)*t2+ (double) (y-side)*t1;
      G = 1.0/(2.0*PI*Xvar*Yvar)*pow(a, (double) ScaleNUM-s)*exp(-0.5*((X*X)/(Xvar*Xvar)+(Y*Y)/(Yvar*Yvar)));
      Gr->data[x][y] = G*cos(2.0*PI*u0*X);
      Gi->data[x][y] = G*sin(2.0*PI*u0*X);
    }
  }
  
  /* if flag = 1, then remove the DC from the filter */
  if (flag == 1) {
    m = 0;
    for (x=0;x<2*side+1;x++)
      for (y=0;y<2*side+1;y++)
        m += Gr->data[x][y];
    
    m /= pow((double) 2.0*side+1, 2.0);
    
    for (x=0;x<2*side+1;x++)
      for (y=0;y<2*side+1;y++)
        Gr->data[x][y] -= m;
  }       
}

void CalculateGaborRotationInvariants_ToProbeCenter(double **Features_Real, double **Features_Imginary, int image_X, int image_Y,  int ScaleNUM, int OrientationNUM )
{
  int    i, j ;
  float  theta ;
  double *Temp_Real, *Temp_Imginary ;
  float  ox, oy, ShiftedPhase ;
  int    s, n, ShtA, ShtB ;
  float  a, b ;

  printf("Invariant features...\n") ;

  /* Gabor features for the polarized image */
  Temp_Real     = Dalloc1d(2*OrientationNUM);
  Temp_Imginary = Dalloc1d(2*OrientationNUM);

  /* get Gabor features in XY domain */
  for(i=0; i<image_X; i++)
    for(j=0; j<image_Y; j++)
      {
        ox = i-CenterOfProbe.x ;
        oy = j-CenterOfProbe.y ;
        
        theta = atan2(oy, ox) ;
        
        if( theta<0 )  theta = 2*PI+theta ;
        
        ShiftedPhase = theta/(PI/OrientationNUM) ;
        a = 1 - (ShiftedPhase - (int)ShiftedPhase) ;
        b = 1 - a ;
        
        /* remember */
        for( s=0; s<ScaleNUM; s++)
          {
            for( n=0; n<OrientationNUM; n++)
              {
                Temp_Real[n]     = Features_Real[i*image_Y+j][s*OrientationNUM+n] ;
                Temp_Imginary[n] = Features_Imginary[i*image_Y+j][s*OrientationNUM+n] ;
                
                Temp_Real[n+OrientationNUM]     =  Temp_Real[n] ;
                Temp_Imginary[n+OrientationNUM] = -Temp_Imginary[n] ;
              }
            
            for( n=0; n<OrientationNUM; n++)
              {
                ShtA = (n+(int)ShiftedPhase  )%(2*OrientationNUM) ;
                ShtB = (n+(int)ShiftedPhase+1)%(2*OrientationNUM) ;
                
                Features_Real[i*image_Y+j][s*OrientationNUM+n]     = a*Temp_Real[ ShtA ] + b*Temp_Real[ ShtB ] ;
                Features_Imginary[i*image_Y+j][s*OrientationNUM+n] = a*Temp_Imginary[ ShtA ] + b*Temp_Imginary[ ShtB ] ;
              }
          }
      }
  
  
  /* free */
  free(Temp_Real) ;
  free(Temp_Imginary);
}

/* save the fetures as images to view the discreminant */
void Put2DGaborResponseInto3DArray(double **Features_Real, double **Features_Imginary, unsigned char **mask, int image_X, int image_Y,int StartScale,int EndScale,int StartOrient,int EndOrient,int ScaleNUM, int OrientationNUM,int SliceIndex, int method)
{
  int h, w, s, n ;
  unsigned char **Fea_img ;
  float mag ;
  int   value ;
  char  FileName[250] ;
  float minImaginary, maxImaginary, minReal, maxReal;
  float rangeImaginary, rangeReal;
  float translation, factor;  // added on 10/26/2011
  //translation=12.8/pow(1.25, 1.0-MAX((float)image_X, (float)image_Y)/256.0);
  //translation=12.8*MIN((float)image_X, (float)image_Y)/256.0;
  //factor=128.0/translation;
  translation=12.8; //12.8;
  factor=10.0; //10.0;
  //translation=128.0*pow(1.25, (MAX((float)image_X, (float)image_Y)/64.0-1.0));
  //factor=128.0/translation;
  //printf("translation=%f, factor=%f\n", translation, factor);
  
  
  
  
  //printf("Save Features to image...\n");
  /* apply for space */
  Fea_img   = UCalloc2d(image_X,image_Y);

  for(s=StartScale; s<=EndScale; s++)
    for(n=StartOrient; n<=EndOrient; n++)
      {// for...s...for...n...
          // Yangming commented this section on August 24, 2010 (begin)
           // for reasons:
           //      1) save computation time
           //      2) more importantly, convert the (Imaginary+12.8)*10 into (Imaginary-min(Imaginary))/(max(Imaginary)-min(Imaginary)), and the same thing for Real
           //      3) mask should not have been applied here, but later on when saving image files.
        
        if (method==1)
        { // if method==1
        // features  
        
        for(h=0;h<image_X;h++)
          for(w=0;w<image_Y;w++) 
            {
            if (mask[h][w]>0)
              {
              mag = Features_Imginary[h*image_Y+w][s*OrientationNUM+n] + translation ;  // changed on 10/26/2011
              //mag = Features_Imginary[h*image_Y+w][s*OrientationNUM+n] + 12.8 ;
              //mag = Features_Imginary[h*image_Y+w][s*OrientationNUM+n] + 8.0 ;
              //mag = fabs(Features_Imginary[h*image_Y+w][s*OrientationNUM+n]);
              if( mag<0 ) mag = 0 ;
              
              value = static_cast<int>(mag*factor);  // changed on 10/26/2011
              //value = mag*10.0;
              //value = mag * 16.0;
              if( value>255 ) value = 255 ;
              }
            else
              value=0;

            if (HoriOrVert==AlongZAxis)
              featureImage[(s*OrientationNUM+n)*2][SliceIndex][h][w]=value;
            else if (HoriOrVert==AlongYAxis)
              featureImage[(s*OrientationNUM+n)*2][h][w][SliceIndex]=value;
            else
              featureImage[(s*OrientationNUM+n)*2][h][SliceIndex][w]=value;
            }



        for(h=0;h<image_X;h++)
          for(w=0;w<image_Y;w++) 
            {
            if (mask[h][w]>0)
              {
              //mag = Features_Real[h*image_Y+w][s*OrientationNUM+n];
              mag = Features_Real[h*image_Y+w][s*OrientationNUM+n];
              if( mag<0 ) mag = 0 ;
              
              value = static_cast<int>(mag*factor);  // changed on 11/26/2011
              //value = mag*10.0;
              //value = mag*16.0;
              if( value>255 ) value = 255 ;
              }
            else
              value=0;
            
              //Fea_img[h][w] = value ;
            if (HoriOrVert==AlongZAxis)
              featureImage[(s*OrientationNUM+n)*2+1][SliceIndex][h][w]=value;
        else if (HoriOrVert==AlongYAxis)
              featureImage[(s*OrientationNUM+n)*2+1][h][w][SliceIndex]=value;
        else
              featureImage[(s*OrientationNUM+n)*2+1][h][SliceIndex][w]=value;
            }
          }  // if method==1
          // Yangming commented this section on August 24, 2010 (end)
        
        
        else if (method==2)
        { // if method==2
        // New part by Yangming on August 24, 2010 (start)
        minImaginary = 10000.0;
        minReal = 10000.0;
        maxImaginary = -10000.0;
        maxReal = -10000.0;
        for(h=0;h<image_X;h++)
          for(w=0;w<image_Y;w++) 
            { //for...h...for...w...
            //if (mask[h][w]>0)
              //{ // if in the foreground
              if (Features_Imginary[h*image_Y+w][s*OrientationNUM+n]>maxImaginary)
                maxImaginary = Features_Imginary[h*image_Y+w][s*OrientationNUM+n];
              if (Features_Imginary[h*image_Y+w][s*OrientationNUM+n]<minImaginary)
                minImaginary = Features_Imginary[h*image_Y+w][s*OrientationNUM+n];
              if (Features_Real[h*image_Y+w][s*OrientationNUM+n]>maxReal)
                maxReal = Features_Real[h*image_Y+w][s*OrientationNUM+n];
              if (Features_Real[h*image_Y+w][s*OrientationNUM+n]<minReal)
                minReal = Features_Real[h*image_Y+w][s*OrientationNUM+n];
                
              if (HoriOrVert==AlongZAxis)
                {
                featureImageFloat[(s*OrientationNUM+n)*2][SliceIndex][h][w]=Features_Imginary[h*image_Y+w][s*OrientationNUM+n];
                featureImageFloat[(s*OrientationNUM+n)*2+1][SliceIndex][h][w]=Features_Real[h*image_Y+w][s*OrientationNUM+n];
                }
              else if (HoriOrVert==AlongYAxis)
                {
                featureImageFloat[(s*OrientationNUM+n)*2][h][w][SliceIndex]=Features_Imginary[h*image_Y+w][s*OrientationNUM+n];
                featureImageFloat[(s*OrientationNUM+n)*2+1][h][w][SliceIndex]=Features_Real[h*image_Y+w][s*OrientationNUM+n];
                }
              else 
                {
                featureImageFloat[(s*OrientationNUM+n)*2][h][SliceIndex][w]=Features_Imginary[h*image_Y+w][s*OrientationNUM+n];
                featureImageFloat[(s*OrientationNUM+n)*2+1][h][SliceIndex][w]=Features_Real[h*image_Y+w][s*OrientationNUM+n];
                }
              //} // if in the foreground
            } //for...h...for...w...
        
        if (HoriOrVert==AlongZAxis)
          {
          minImaginaryHori[SliceIndex][s][n] = minImaginary;
          maxImaginaryHori[SliceIndex][s][n] = maxImaginary;
          minRealHori[SliceIndex][s][n] = minReal;
          maxRealHori[SliceIndex][s][n] = maxReal;
          }
        else  
          {
          minImaginaryVert[SliceIndex][s][n] = minImaginary;
          maxImaginaryVert[SliceIndex][s][n] = maxImaginary;
          minRealVert[SliceIndex][s][n] = minReal;
          maxRealVert[SliceIndex][s][n] = maxReal;
          }
        } // if method==2
        /*
        rangeImaginary = maxImaginary - minImaginary;
        rangeReal = maxReal - minReal;
        printf("s=%d, o=%d, imag=(%f-%f:%f), real=(%f-%f:%f)\n", s, n, minImaginary, maxImaginary, rangeImaginary, minReal, maxReal, rangeReal);
        
        for(h=0;h<image_X;h++)
          for(w=0;w<image_Y;w++) 
            { //for...h...for...w...
            if (mask[h][w]>0)
              value = (int)(255.0*(Features_Imginary[h*image_Y+w][s*OrientationNUM+n]-minImaginary)/rangeImaginary);
            else
              value = 0;
              
            if (HoriOrVert)
              featureImage[(s*OrientationNUM+n)*2][SliceIndex][h][w]=value;
            else
              featureImage[(s*OrientationNUM+n)*2][h][w][SliceIndex]=value;
             
            if (mask[h][w]>0)
              value = (int)(255.0*(Features_Real[h*image_Y+w][s*OrientationNUM+n]-minReal)/rangeReal);  
            else
              value = 0;
            if (HoriOrVert)
              featureImage[(s*OrientationNUM+n)*2+1][SliceIndex][h][w]=value;
            else
              featureImage[(s*OrientationNUM+n)*2+1][h][w][SliceIndex]=value;
            } //for...h...for...w...
        // New part by Yangming on August 24, 2010 (end)
        */
      } // for...s...for...n...
  
  // free 
  UCfree2d(Fea_img, image_X);
}

void Construct3dImg(unsigned char*** F_img3d,int Scale,int Orientation,int image_X,int image_Y,int SliceNum,int RealOrImag)
{
  FILE  *fp;
  int   i;
  char filename[500],command[500];
  unsigned char** UC_img2d;
  int h,w;  

  printf("\n\nConsturct3dFeatureImage...\n");
  UC_img2d=UCalloc2d(image_X,image_Y);

  for (i=0;i<SliceNum;i++)
    {
      if (RealOrImag)
        if (HoriOrVert==AlongZAxis)
          sprintf(filename, "%s_Hori_Freal.%.3d.%d_%d.raw", ResultDirectory,i,Scale,Orientation);
        else
          sprintf(filename, "%s_Vert_Freal.%.3d.%d_%d.raw", ResultDirectory,i,Scale,Orientation);
      else
        if (HoriOrVert==AlongZAxis)
          sprintf(filename, "%s_Hori_Fimag.%.3d.%d_%d.raw", ResultDirectory,i,Scale,Orientation);
        else
          sprintf(filename, "%s_Vert_Fimag.%.3d.%d_%d.raw", ResultDirectory,i,Scale,Orientation);
      
      printf("Reading %s...\n",filename);
      ReadImg2D(filename,UC_img2d,image_X,image_Y);

      sprintf(command,"rm %s",filename);
      system(command);
      
      if (HoriOrVert==AlongZAxis)
        {
          for (h=0;h<image_X;h++)
            for (w=0;w<image_Y;w++)          
              F_img3d[i][h][w]=UC_img2d[h][w];           
        }
      else
        {
          for (h=0;h<image_X;h++)
            for (w=0;w<image_Y;w++)          
              F_img3d[h][w][i]=UC_img2d[h][w];
        }
    }  

  
  UCfree2d(UC_img2d,image_X);

}

void ReadImg2D(const char* filename,unsigned char **data,int image_X,int image_Y)
{
  FILE  *fp;
  int   i;

  /* write the smoothed image */
  fp=fopen(filename,"r");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s!", filename);
    exit(1);
  }
  for(i=0;i<image_X;i++)
    fread(data[i],1,image_Y,fp);
  fclose(fp);
}

void WriteImg2D(const char* filename, unsigned char **data, int image_X, int image_Y)
{
  FILE  *fp;
  int   i;

  /* write the smoothed image */
  fp=fopen(filename,"w");
  if (!fp) {
    fprintf(stderr, "Failed to open file %s for writing!", filename);
    exit(1);
  }
  for(i=0;i<image_X;i++)
    fwrite(data[i],1,image_Y,fp);
  fclose(fp);
}



/* edge detection */
void DetectEdgeFromGaborFeatureVectors(unsigned char **Img_edge, int image_X, int image_Y, int FeaNUM, double **Features_Real, double **Features_Imginary, float CannyGuassianSigma)
{
  /* Canny */
  unsigned char *XY_magnitude, *XY_orient, *XY_derivative_mag; 
                                                    /* mag of del G before non-maximum suppression */
  int  i,j ;

  printf("\nDetect edge by Canny operator...\n");
  /* Canny part: memory applciation*/
  XY_magnitude     = (unsigned char *)calloc(image_X*image_Y,1) ;
  XY_orient        = (unsigned char *)calloc(image_X*image_Y,1) ;
  XY_derivative_mag= (unsigned char *)calloc(image_X*image_Y,1) ;


  /* edge detection*/
  /* Canny */
  canny_core(CannyGuassianSigma, image_Y, image_X, FeaNUM, Features_Real, Features_Imginary, XY_derivative_mag, XY_magnitude, XY_orient);

  for(i=0; i<image_X; i++)
    for(j=0; j<image_Y; j++)
      if (XY_magnitude[i*image_Y+j]>0.4)
        Img_edge[i][j] = 255;
      else
        Img_edge[i][j] = 0;
  
  
  /* 2D: free */
  free(XY_orient) ;
  free(XY_derivative_mag) ;
  free(XY_magnitude) ;
}

void CalculateAndSaveGaborFilters_ForVisualEffect(int side, int ScaleNUM, int OrientationNUM, int flag, double Ul, double Uh)
{
   Matrix *Gr, *Gi ;
   int s, n, i, j;
   double *GaborFilterBank, real, imaginary ;
   

   /* applying for space */
   CreateMatrix(&Gr, 2*side+1, 2*side+1);
   CreateMatrix(&Gi, 2*side+1, 2*side+1);   
   GaborFilterBank = (double *) calloc((2*side+1)*(2*side+1)*ScaleNUM*OrientationNUM*2, sizeof(double));
   
   /* action ... */
   for(s=0; s<ScaleNUM; s++)
     for(n=0; n<OrientationNUM; n++)
       {
         printf("GetGabor (%d, %d)\n", s, n) ;
     
         GaborFunction(Gr, Gi, s+1, n+1, Ul, Uh, ScaleNUM, OrientationNUM, flag) ;
         
         real = 0 ;
         imaginary = 0 ;
         for(i=0; i<(2*side+1); i++)
           for(j=0; j<(2*side+1); j++)
             {
               GaborFilterBank[2*(2*side+1)*(2*side+1)*(s*OrientationNUM+n) + i*(2*side+1)+j] = Gr->data[i][j]; 
               GaborFilterBank[2*(2*side+1)*(2*side+1)*(s*OrientationNUM+n) + (2*side+1)*(2*side+1)+i*(2*side+1)+j] = Gi->data[i][j]; 
               
               real += fabs(Gr->data[i][j]) ;
               imaginary += fabs(Gi->data[i][j]) ;
             }
       }   
   WriteGaborFiltersToFile(GaborFilterBank, side, ScaleNUM, OrientationNUM) ;
   
   /* free */
   FreeMatrix(Gr);
   FreeMatrix(Gi);
   free(GaborFilterBank);
}

void WriteGaborFiltersToFile(double *GaborFilterBank, int side, int ScaleNUM, int OrientationNUM)
{
  FILE *stream;
  int s, n, i, j ;
  
  stream = fopen( "Gabor.view.text", "w" );
  
  for(s=0; s<ScaleNUM; s++)
    for(n=0; n<OrientationNUM; n++)
      {
        fprintf(stream, "Real (S,O)=(%d, %d)\n", s, n) ;
        for(i=0; i<(2*side+1); i++)
          {
            for(j=0; j<(2*side+1); j++)
              fprintf(stream, "%7.4f  ", GaborFilterBank[2*(2*side+1)*(2*side+1)*(s*OrientationNUM+n)+i*(2*side+1)+j] ) ;
            fprintf(stream, "\n") ;
          }
        fprintf(stream, "\n") ;
        
        
        fprintf(stream, "Imginary (S,O)=(%d, %d)\n", s, n) ;
        for(i=0; i<(2*side+1); i++)
          {
            for(j=0; j<(2*side+1); j++)
              fprintf(stream, "%7.4f  ", GaborFilterBank[2*(2*side+1)*(2*side+1)*(s*OrientationNUM+n)+(2*side+1)*(2*side+1)+i*(2*side+1)+j]) ;
            fprintf(stream, "\n") ;
          }
        fprintf(stream, "\n") ;
      }
  
  fclose( stream );
}



void canny_core(float s, int cols, int rows, int FeaNUM, double **Features_Real,double **Features_Imginary,unsigned char *derivative_mag,unsigned char *magnitude,unsigned char *orientation)
{
  int      filter_width;            /* length of 1-D gaussian mask */
  float    *gsmooth_x,*gsmooth_y,mag,value;
  float    *derivative_x,*derivative_y;
  int      i,k,n;                   /* counters */
  int      t;                       /* temp. grad magnitude variable */
  double   a,b,c,d,g0;              /* mask generation intermediate vars*/
  double   ux,uy;
  double   t1,t2;
  double   *t1_Real, *t1_Imginary, *t2_Real, *t2_Imginary;
  double   G[20],dG[20],D2G[20];    /*gaussian & derivative filter masks*/
  double   gc,gn,gs,gw,ge,gnw,gne,gsw,gse;
  int      picsize,jstart,jlimit;
  int      ilimit;
  register int jfactor;
  int      kfactor1,kfactor2;
  int      kfactor;
  register int cindex,nindex,sindex,windex,eindex,nwindex,neindex,swindex,seindex;
  int      low=1,high=255;          /* tracker hysteresis parameters */
  int      mag_overflow_count=0;    /* used to measure how oft mag array overflows */
  int      mag_underflow_count=0;   /* used to measure how oft mag array underflows */

  double weightG, weightDG ;  /*normalizing ...*/
  int      forder ;


  picsize=cols*rows;    /* picture area */

  /* calc coeffs for 1-dimensional G, dG/dn and Delta-squared G filters */
  for(n=0; n<20; ++n)
    {
      a=gaussian((double)n,s);
      if (a>0.005 || n<2)
        {
          b=gaussian((double)n-0.5,s);
          c=gaussian((double)n+0.5,s);
          d=gaussian((double)n,s*0.5);
          /*fprintf(stderr,"a,b,c: %lf,%lf,%lf\n",a,b,c);*/
          G[n]=(a+b+c)/3/(6.283185*s*s);
          dG[n]=c-b;
          D2G[n]=1.6*d-a; /* DOG */
          /* fprintf(stderr,"G[%d]: %lf\n",n,G[n]);
             fprintf(stderr,"dG[%d]: %lf\n",n,dG[n]); 
             fprintf(stderr,"D2G[%d]: %lf\n",n,D2G[n]);*/     
        }
      else break;
    }
  filter_width=n; /* printf("n=%d\n", n); */
  

  /*normalizing ...*/
  weightG  = 0.5*G[0] ; 
  weightDG = 0 ; 
  for(k=1;k<filter_width;++k)
    {
      weightG  +=  G[k] ; 
      weightDG += fabs(dG[k]) ;
    }
  for(k=0;k<filter_width;++k)
    {
      G[k]  /= weightG ; 
      dG[k] /= weightDG ; 
      /*fprintf(stderr,"G[%d]: %lf\n",k,G[k]);
        fprintf(stderr,"dG[%d]: %lf\n",k,dG[k]);*/
    }
  
  
  /*fprintf(stderr,"canny_core: smooth pic\n");*/
  /* allocate space for gaussian smoothing arrays */
  if ((gsmooth_x=(float *)calloc(picsize,sizeof(float)))==(float *)NULL)
    {
      fprintf(stderr,"can't alloc gsmooth_x\n");
      exit(0);
    }
  if ((gsmooth_y=(float *)calloc(picsize,sizeof(float)))==(float *)NULL)
    {
      fprintf(stderr,"can't alloc gsmooth_y\n");
      exit(0);
    }
  
  
  /* apply memory for combined neighboring Gabor feature vectors*/
  t1_Real     = Dalloc1d(FeaNUM) ;
  t1_Imginary = Dalloc1d(FeaNUM) ;
  t2_Real     = Dalloc1d(FeaNUM) ;
  t2_Imginary = Dalloc1d(FeaNUM) ;
  
  
  
  /* produce x- and y- convolutions with gaussian */
  ilimit=cols-(filter_width-1);
  jstart=cols*(filter_width-1);
  jlimit=cols*(rows-(filter_width-1));
  for (i=filter_width-1;i<ilimit;++i)
    {
      for(jfactor=jstart;
          jfactor<jlimit;
          jfactor+=cols)
        {
          cindex=i+jfactor;
          
          for( forder=0; forder<FeaNUM; forder++)
            {
              t1_Real[forder]     = Features_Real[cindex][forder]*G[0];
              t1_Imginary[forder] = Features_Imginary[cindex][forder]*G[0];
              
              t2_Real[forder]     = t1_Real[forder];
              t2_Imginary[forder] = t1_Imginary[forder];
            }
          
          for(k=1,kfactor1=cindex-cols,kfactor2=cindex+cols;
              k<filter_width;
              k++,kfactor1-=cols,kfactor2+=cols)
            {
              for( forder=0; forder<FeaNUM; forder++)
                {
                  t1_Real[forder]     += G[k]*(Features_Real[kfactor1][forder]     + Features_Real[kfactor2][forder]);
                  t1_Imginary[forder] += G[k]*(Features_Imginary[kfactor1][forder] + Features_Imginary[kfactor2][forder]);
                  
                  t2_Real[forder]     += G[k]*(Features_Real[cindex-k][forder]     + Features_Real[cindex+k][forder]);
                  t2_Imginary[forder] += G[k]*(Features_Imginary[cindex-k][forder] + Features_Imginary[cindex+k][forder]);
        }
            }

          gsmooth_x[cindex] = 0 ;
          gsmooth_y[cindex] = 0 ;

          
          //Zhan modified 6.11.
          mag = t1_Imginary[3];
          if( mag<0 ) mag = 0 ;
              
          value = mag*10.0;
          if( value>255 ) value = 255 ;              
          gsmooth_x[cindex]=value;

          mag = t2_Imginary[3];
          if( mag<0 ) mag = 0 ;
              
          value = mag*10.0;
          if( value>255 ) value = 255 ;              
          gsmooth_y[cindex]=value;

          /*
          for( forder=0; forder<FeaNUM; forder++)
            {
              gsmooth_x[cindex] += t1_Real[forder]*t1_Real[forder] + t1_Imginary[forder]*t1_Imginary[forder] ; 
              gsmooth_y[cindex] += t2_Real[forder]*t2_Real[forder] + t2_Imginary[forder]*t2_Imginary[forder] ; 
            }
          */
        }
    }
  
  /* allocate space for gradient arrays */
  /*fprintf(stderr,"canny_core: find grad\n");*/

  
  if ((derivative_x=(float *)calloc(picsize,sizeof(float)))==(float *)NULL)
    {
      fprintf(stderr,"can't alloc x\n");
      exit(0);
    }
  /* produce x and y convolutions with derivative of gaussian */
  
  for (i=filter_width-1;i<ilimit;++i)
    {
      for(jfactor=jstart;jfactor<jlimit;jfactor+=cols)
        {
          t1=0;
          cindex=i+jfactor;
          for (k=1;k<filter_width;++k)
            t1+=dG[k]*(gsmooth_x[cindex-k]-gsmooth_x[cindex+k]);
          derivative_x[cindex] = (float)t1;
        }
    }
  free(gsmooth_x);
  if ((derivative_y=(float *)calloc(picsize,sizeof(float)))==(float *)NULL)
    {
      fprintf(stderr,"can't alloc y\n");
      exit(0);
    }
  
  for (i=n;i<cols-n;++i)
    {
      for(jfactor=jstart;jfactor<jlimit;jfactor+=cols)
        {
          t2=0;
          cindex=i+jfactor;
          for (k=1,kfactor=cols;k<filter_width;k++,kfactor+=cols)
            t2+=dG[k]*(gsmooth_y[cindex-kfactor]-gsmooth_y[cindex+kfactor]);
          derivative_y[cindex] = (float)t2;
        }
    }
  free(gsmooth_y);
    
  
  /* non-maximum suppression (4 cases for orientation of line of max slope) */
  /* fprintf(stderr,"canny_core: non-maximum suppression\n"); */
  ilimit=cols-filter_width;
  jstart=cols*filter_width;
  jlimit=cols*(rows-filter_width);
  
  for (i=filter_width;i<ilimit;++i)
    {
      for (jfactor=jstart;jfactor<jlimit;jfactor+=cols)
        {
          /* calculate current indeces */
          cindex=i+jfactor;
          nindex=cindex-cols;
          sindex=cindex+cols;
          windex=cindex-1;
          eindex=cindex+1;
          nwindex=nindex-1;
          neindex=nindex+1;
          swindex=sindex-1;
          seindex=sindex+1;
          ux=derivative_x[cindex];
          uy=derivative_y[cindex];
          gc=hypotenuse(ux,uy);
          /* scale gc to fit into an unsigned char array */
          t=(int)(gc*0.01); /* 1.5 Oct 2001; 3.0 July 26, 1999 */
          /*fprintf(stderr,"canny_core: i,j=(%d,%d), t=%lf\n",i,jfactor/cols,t);*/
          derivative_mag[cindex]=(t<256 ? t : 255);
          gn=hypotenuse(derivative_x[nindex],derivative_y[nindex]);
          gs=hypotenuse(derivative_x[sindex],derivative_y[sindex]);
          gw=hypotenuse(derivative_x[windex],derivative_y[windex]);
          ge=hypotenuse(derivative_x[eindex],derivative_y[eindex]);
          gne=hypotenuse(derivative_x[neindex],derivative_y[neindex]);
          gse=hypotenuse(derivative_x[seindex],derivative_y[seindex]);
          gsw=hypotenuse(derivative_x[swindex],derivative_y[swindex]);
          gnw=hypotenuse(derivative_x[nwindex],derivative_y[nwindex]);
          if (ux*uy>0)
            {
              if(fabs(ux)<fabs(uy))
                {
                  if((g0=fabs(uy*gc))< fabs(ux*gse+(uy-ux)*gs) ||g0<=fabs(ux*gnw+(uy-ux)*gn))
                    continue;
                }
              else
                {
                  if((g0=fabs(ux*gc))< fabs(uy*gse+(ux-uy)*ge)||g0<=fabs(uy*gnw+(ux-uy)*gw))
                    continue;
                }
            }
          else
            {
              if(fabs(ux)<fabs(uy))
                {
                  if((g0=fabs(uy*gc))< fabs(ux*gne-(uy+ux)*gn) ||g0<=fabs(ux*gsw-(uy+ux)*gs))
                    continue;
                }
              else
                {
                  if((g0=fabs(ux*gc))< fabs(uy*gne-(ux+uy)*ge) ||g0<=fabs(uy*gsw-(ux+uy)*gw))
                    continue;
                }
            }
          /* seems to be a good scale factor */
          magnitude[cindex]=derivative_mag[cindex];
          /* pi*40 ~= 128 - direction is (thought of as) a signed byte */
          orientation[cindex]=(unsigned char)((PI+atan2(ux, uy))*ORIENT_SCALE);
          /* if (orientation[cindex]>(3.1416*ORIENT_SCALE)) 
                orientation[cindex]=orientation[cindex];
                if (orientation[cindex]<0)
            orientation[cindex]=orientation[cindex];
                                                                                          float t1,t2,t3,t4;
                                                                                          unsigned char ut1,ut2,ut3,ut4;
                                                                                          
                                                                                          t1=atan2(1,2);
                                                                                          t2=atan2(-1,2);
                                                                                          t3=atan2(1,-2);
                                                                                          t4=atan2(-1,-2);
                                                                                          ut1=(unsigned char)(t1*ORIENT_SCALE);  
                                                                                          ut2=(unsigned char)(t2*ORIENT_SCALE); 
                                                                                          ut3=(unsigned char)(t3*ORIENT_SCALE); 
                                                                                          ut4=(unsigned char)(t4*ORIENT_SCALE);  
          */ 
        }
    } 

  free(derivative_x);
  free(derivative_y);
}

/*#define NUM_Ori  360*/
int     NUM_Rad, NUM_Ori ;
void PolarizeImg( unsigned char **UC_img, int image_X, int image_Y )
{
  int i, j, d0, d1, c, r, angle ;
  unsigned char **Polar_Img ;
  Fvector2d  *AxisOfEllipseWindow, CurrentP ;
  float mag, theta ;


  printf("Polariing image ... \n") ;

  /* apply for space */
  NUM_Rad = image_X ;
  NUM_Ori = image_Y ;
  Polar_Img = UCalloc2d(NUM_Rad, NUM_Ori);


  /* process */
  for( r=0; r<NUM_Rad; r++)
    for( angle=0; angle<NUM_Ori; angle++)
      {
        theta = (angle*180.0/NUM_Ori+90.0)/180.0*PI ;
        
        i = static_cast<int>(CenterOfProbe.x + r*cos( theta ));
        j = static_cast<int>(CenterOfProbe.y + r*sin( theta ));
        
        if( i<image_X && j<image_Y &&  i>=0 && j>=0 )
          Polar_Img[NUM_Rad-1-r][NUM_Ori-1-angle] = UC_img[i][j] ;
        else
          Polar_Img[NUM_Rad-1-r][NUM_Ori-1-angle] = 0 ;
      }
  WriteImg2D("PolarImg.raw", Polar_Img, NUM_Rad, NUM_Ori) ;

  printf("finished\n") ;

  /* free */
  UCfree2d(Polar_Img, NUM_Rad) ;
}

void GetGaborFeatureOfOneImage_ThroughPolarCoordinate(double **Features_Real, double **Features_Imginary, unsigned char **UC_img, int image_X, int image_Y,  int side, double Ul, double Uh, int ScaleNUM, int OrientationNUM, int flag)
{
  int i, j, d0, d1, c, r, angle ;
  unsigned char **Polar_Img ;
  float mag, theta ;
  double **PolarFeatures_Real, **PolarFeatures_Imginary ;
  float  ox, oy ;
  int    forder, FeaNUM ;
  
  
  printf("Polariing image ... \n") ;
  /* apply for space */
  NUM_Rad = image_X ;
  NUM_Ori = image_Y ;
  Polar_Img = UCalloc2d(NUM_Rad, NUM_Ori);
  
  
  /* process to get polar representation */
  for( r=0; r<NUM_Rad; r++)
    for( angle=0; angle<NUM_Ori; angle++)
      {
        theta = (angle*180.0/NUM_Ori+90.0)/180.0*PI ;
        
        i = static_cast<int>(CenterOfProbe.x + r*cos( theta ));
        j = static_cast<int>(CenterOfProbe.y + r*sin( theta ));
        
        if( i<image_X && j<image_Y &&  i>=0 && j>=0 )
          Polar_Img[NUM_Rad-1-r][NUM_Ori-1-angle] = UC_img[i][j] ;
        else
          Polar_Img[NUM_Rad-1-r][NUM_Ori-1-angle] = 0 ;
      }

  printf("finished\n") ;
  
  
  /* Gabor features for the polarized image */
  FeaNUM = ScaleNUM*OrientationNUM ;
  PolarFeatures_Real     = DDalloc2d(NUM_Rad*NUM_Ori, FeaNUM) ;
  PolarFeatures_Imginary = DDalloc2d(NUM_Rad*NUM_Ori, FeaNUM) ;
  
  /* Calculate the Gabor features of the polarized image */
  GetGaborFeatureOfOneImage(PolarFeatures_Real, PolarFeatures_Imginary, Polar_Img, NUM_Rad, NUM_Ori,  side, Ul, Uh, ScaleNUM, OrientationNUM, flag) ; 
  
  
  /* get Gabor features in XY domain */
  for(i=0; i<image_X; i++)
    for(j=0; j<image_Y; j++)
      {
        ox = i-CenterOfProbe.x ;
        oy = j-CenterOfProbe.y ;
        r = static_cast<int>(sqrt( ox*ox + oy*oy ));
        
        theta = atan2(oy, ox) ;
        
        if( theta<0 )  theta = 2*PI+theta ;
        
        angle = (int)( (theta/PI*180.0-90)*NUM_Ori/180.0 ) ;
        
        /*printf("ox=%f  angle=%f\n", ox, theta) ;*/
        /*printf("(%d,%d) x=%f  y=%f\n", i,j, CenterOfProbe.x+r*cos(theta), CenterOfProbe.y+r*sin(theta)) ;*/
        
        
        if( r<NUM_Rad && angle<NUM_Ori &&  r>=0 && angle>=0 )
          for( forder=0; forder<FeaNUM; forder++)
            {
              Features_Real[i*image_Y+j][forder] = PolarFeatures_Real[(NUM_Rad-1-r)*NUM_Ori+(NUM_Ori-1-angle)][forder] ;
              Features_Imginary[i*image_Y+j][forder] = PolarFeatures_Imginary[(NUM_Rad-1-r)*NUM_Ori+(NUM_Ori-1-angle)][forder] ;
            }
        else
          for( forder=0; forder<FeaNUM; forder++)
            {
              Features_Real[i*image_Y+j][forder] = 0 ;
              Features_Imginary[i*image_Y+j][forder] = 0 ;
            }
      }
  
  
  /* free */
  UCfree2d(Polar_Img, NUM_Rad) ;
  Dfree2d(PolarFeatures_Real, NUM_Rad*NUM_Ori) ;
  Dfree2d(PolarFeatures_Imginary, NUM_Rad*NUM_Ori) ;
}

void Write3DImage(const char* filename,unsigned char*** Image,int x_size,int y_size,int z_size)
{
  int i,j,k;
  FILE * fp;

  fp=fopen(filename,"w");
  
  for (k=0;k<z_size;k++)
    for (i=0;i<x_size;i++)
      fwrite(Image[k][i],1,y_size,fp);  
  fclose(fp);
}


double **DDalloc2d(int i_size,int j_size)
{
  double **array;
  int i,j;

  array=(double **) calloc(i_size,sizeof(double *));

  for(i=0;i<i_size;i++)
    array[i]=(double *) calloc(j_size,sizeof(double ));

  return(array);
}



void generateMask(unsigned char***image, Ivector3d imageSize, unsigned char***mask, int foregroundThre, int fillHoleOrNot, int regionGrowingOrNot)
{
  int i,j,k;
  int ii,jj,kk;
  int searchRangeX=(int)((float)imageSize.x/20.0);
  int searchRangeY=(int)((float)imageSize.y/20.0);
  int searchRangeZ=(int)((float)imageSize.z/20.0);
  int searchRangeXY=MAX(searchRangeX,searchRangeY);
  int maxCoorInRows[imageSize.z][imageSize.x];  // to record the maximum coordinate of the point in the foreground
  int minCoorInRows[imageSize.z][imageSize.x];
  int maxCoorInColumns[imageSize.z][imageSize.y];
  int minCoorInColumns[imageSize.z][imageSize.y];
  int maxCoorVertical[imageSize.x][imageSize.y];
  int minCoorVertical[imageSize.x][imageSize.y];
  int numTimes;
  int RegionPointNumThre;
  int left, right, up, down;
  
  printf("first, keep all voxels with intensity greater than %d\n", foregroundThre);
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
        { //for...for...for...
        mask[k][i][j]=0;
        if (image[k][i][j]>foregroundThre)
          mask[k][i][j]=255;
        }  // for...for...for...
    
    
  
  if (regionGrowingOrNot==YYES)
  {  // if (fillHoleOrNot)
  // detected connected regions by region growing, and remove regions that contain less than RegionPointNumThre voxels.
  if ( MIN(imageSize.x,imageSize.y)<65 )        RegionPointNumThre = 80;//*imageSize.z;
  else if ( MIN(imageSize.x,imageSize.y)<120 )  RegionPointNumThre = 160;//*imageSize.z;
  else if ( MIN(imageSize.x,imageSize.y)<260 )  RegionPointNumThre = 250;//*imageSize.z;
  else                                          RegionPointNumThre = 400;//*imageSize.z;
  
  RegionGrowingInBinaryImage(mask, imageSize, RegionPointNumThre);
  }
  
  
  if (fillHoleOrNot==YYES)
  {
  // then fill in holes
  printf("then fill in holes in the foreground\n");
  // initialize
  for (k=0;k<imageSize.z;k++)
    {
    for (i=0;i<imageSize.x;i++)
      {
      maxCoorInRows[k][i]=0;
      minCoorInRows[k][i]=imageSize.y-1;
      }
    for (j=0;j<imageSize.y;j++)
      {
      maxCoorInColumns[k][j]=0;
      minCoorInColumns[k][j]=imageSize.x-1;
      }
    }
  for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
      {
      maxCoorVertical[i][j]=0;
      minCoorVertical[i][j]=imageSize.z-1;
      } 
  
  
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
      for (j=0;j<imageSize.y;j++)
        { //for...for...for...
        if (image[k][i][j]>foregroundThre)
           {
           if (j>maxCoorInRows[k][i])      maxCoorInRows[k][i]=j;
           if (j<minCoorInRows[k][i])      minCoorInRows[k][i]=j;
           if (i>maxCoorInColumns[k][j])   maxCoorInColumns[k][j]=i;
           if (i<minCoorInColumns[k][j])   minCoorInColumns[k][j]=i;
           if (k>maxCoorVertical[i][j])    maxCoorVertical[i][j]=k;
           if (k<minCoorVertical[i][j])    minCoorVertical[i][j]=k;
           }
        }
  // expand search space a little bit
  for (k=0;k<imageSize.z;k++)
    {
    for (i=0;i<imageSize.x;i++)
       {
       if (minCoorInRows[k][i]-searchRangeY<0)  minCoorInRows[k][i]=0;
       if (maxCoorInRows[k][i]+searchRangeY>imageSize.y-1)   maxCoorInRows[k][i]=imageSize.y-1;    
       }
    for (j=0;j<imageSize.y;j++)
       {
       if (minCoorInColumns[k][j]-searchRangeX<0)  minCoorInColumns[k][j]=0;
       if (maxCoorInRows[k][j]+searchRangeX>imageSize.x-1)   maxCoorInColumns[k][j]=imageSize.x-1;
       }
    }
  for (i=0;i<imageSize.x;i++)
    for (j=0;j<imageSize.y;j++)
    {
    if (minCoorVertical[i][j]-searchRangeZ<0)  minCoorVertical[i][j]=0;
    if (maxCoorVertical[i][j]+searchRangeZ>imageSize.z-1)  maxCoorVertical[i][j]=imageSize.z-1;
    }
  
  
  // fill in holes inside the foreground
  for (numTimes=1;numTimes<=3;numTimes++)
   {  // for num
  for (k=0;k<imageSize.z;k++)
    { // for k
    for (i=0;i<imageSize.x;i++)
      for (j=minCoorInRows[k][i];j<=maxCoorInRows[k][i];j++)
        { // for i,j
        if (mask[k][i][j]==0)
          {
          left=0;right=0;up=0;down=0;
          for (jj=MAX(-j,-searchRangeY);jj<0;jj++)
            if (mask[k][i][j+jj]>0) {left=1; break;}
          for (jj=1;jj<MIN(imageSize.y-j-1, searchRangeY);jj++)
            if (mask[k][i][j+jj]>0) {right=1; break;}
          for (ii=MAX(-i,-searchRangeX);ii<0;ii++)
            if (mask[k][i+ii][j]>0) {up=1; break;}
          for (ii=1;ii<MIN(imageSize.x-i-1, searchRangeX);ii++)
            if (mask[k][i+ii][j]>0) {down=1; break;}
        
          if (i-searchRangeX<0) up=1;
          if (i+searchRangeX>imageSize.x-1) down=1;
          if (j-searchRangeY<0) left=1;
          if (j+searchRangeY>imageSize.y-1) right=1;
          
          if (left+right+up+down==4)
            mask[k][i][j]=255;
          }
        } // for i j
    } // for k
   }// for num
  } // if (fillHoleOrNot)  
}
   
   

float GetGreatestCommonDivisor(float aa, float bb)
{
  int aa_int = (int)floor(aa*2+0.5);
  int bb_int = (int)floor(bb*2+0.5);
  int ll = MIN(aa_int,bb_int);
  
  int temp;
  float gcd;
  for (temp=1;temp<=ll;temp++)
    {
    if ( (aa_int%temp==0) && (bb_int%temp==0) )
       gcd = (float)temp;
    }
  
  gcd /= 2.0;
  
  return gcd;
}
   

   
void zoomInImage(unsigned char **UC_Img, unsigned char **zoomedUC_Img, int origSizeX, int origSizeY, int newSizeX, int newSizeY, int zoomFactor_X, int zoomFactor_Y)
{
  int ii,jj;
  float coorX, coorY;

  Image image;
  image.region.nx    = origSizeX;
  image.region.ny    = origSizeY;
  image.region.nz    = 1;
  image.hdr.datatype = DT_UNSIGNED_CHAR;
  image.imgfmt       = Image::FORMAT_DRAMMS;
  image.img.uc       = &UC_Img;
  image.ownsimg      = false;
 
  for (ii=0;ii<newSizeX;ii++)
    for (jj=0;jj<newSizeY;jj++)  
      {
      coorX = (float)(ii+1)/(float)zoomFactor_X - 1.0;
      coorY = (float)(jj+1)/(float)zoomFactor_Y - 1.0;
      
      if ( (coorX<0.0) || (coorY<0.0) )
         zoomedUC_Img[ii][jj] = 0;
      else
         zoomedUC_Img[ii][jj] = static_cast<unsigned char>(image.value(coorX, coorY)); // bilinear interpolation
     }
}



void linearScalingGaborResponseIntoAttributeValues(int scaleIndex, int orientationIndex, int ScaleNUM, int OrientationNUM, int SliceNUM, int sizeX, int sizeY, int sizeZ, unsigned char ***mask)
{
  float minResponseImag=10000.0f;
  float minResponseReal=10000.0f;
  float maxResponseImag=-10000.0f;
  float maxResponseReal=-10000.0f;
  float rangeResponseImag, rangeResponseReal;
  int sliceIndex;
  
  
  printf("Converting Gabor response into attribute values by linear scaling, for attributes in scale %d and orientation %d\n", scaleIndex,orientationIndex);
  
  printf("scale=%d/%d, orientation=%d/%d\n", scaleIndex, ScaleNUM, orientationIndex, OrientationNUM);
  printf("sliceNUM=%d\n", SliceNUM);
  printf("HoriOrVert=%d\n", HoriOrVert);
  printf("sizeX=%d, sizeY=%d, sizeZ=%d\n", sizeX, sizeY, sizeZ);
  
  if (HoriOrVert==AlongZAxis)
    { // horizental attributes
    for (sliceIndex=0;sliceIndex<SliceNUM;sliceIndex++)
      { // for each slice
      // min and max in imaginary response
      if (minImaginaryHori[sliceIndex][scaleIndex][orientationIndex]<minResponseImag)
         minResponseImag = minImaginaryHori[sliceIndex][scaleIndex][orientationIndex];
      if (maxImaginaryHori[sliceIndex][scaleIndex][orientationIndex]>maxResponseImag)
         maxResponseImag = maxImaginaryHori[sliceIndex][scaleIndex][orientationIndex];
      
      // min and max in real response     
      if (minRealHori[sliceIndex][scaleIndex][orientationIndex]<minResponseReal)
         minResponseReal = minRealHori[sliceIndex][scaleIndex][orientationIndex];
      if (maxRealHori[sliceIndex][scaleIndex][orientationIndex]>maxResponseReal)
         maxResponseReal = maxRealHori[sliceIndex][scaleIndex][orientationIndex];    
      } // for each slice
    } // horizental attributes
  else  
    { // vertical attributes
    for (sliceIndex=0;sliceIndex<SliceNUM;sliceIndex++)
      { // for each slice
      // min and max in imaginary response
      if (minImaginaryVert[sliceIndex][scaleIndex][orientationIndex]<minResponseImag)
         minResponseImag = minImaginaryVert[sliceIndex][scaleIndex][orientationIndex];
      if (maxImaginaryVert[sliceIndex][scaleIndex][orientationIndex]>maxResponseImag)
         maxResponseImag = maxImaginaryVert[sliceIndex][scaleIndex][orientationIndex];
      
      // min and max in real response     
      if (minRealVert[sliceIndex][scaleIndex][orientationIndex]<minResponseReal)
         minResponseReal = minRealVert[sliceIndex][scaleIndex][orientationIndex];
      if (maxRealVert[sliceIndex][scaleIndex][orientationIndex]>maxResponseReal)
         maxResponseReal = maxRealVert[sliceIndex][scaleIndex][orientationIndex];    
      } // for each slice
    } // vertical attributes
    
  maxResponseImag/=2.0;
  minResponseImag/=2.0;
  maxResponseReal/=2.0;
  minResponseReal/=2.0;  
  rangeResponseImag = maxResponseImag - minResponseImag;
  rangeResponseReal = maxResponseReal - minResponseReal;
  printf("rangeResponseImag=%f, rangeResponseReal=%f\n", rangeResponseImag, rangeResponseReal);
  
  int xx,yy,zz;
  for (zz=0;zz<sizeZ;zz++)
    for (xx=0;xx<sizeX;xx++)
      for (yy=0;yy<sizeY;yy++)
        {
        featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy] = 0;
        featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy] = 0;
        if (mask[zz][xx][yy]>0)
          {
          // linear scaling for the whole 3D attribute image in scaleIndex and orientationIndex
          featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy] = (int)(255.0*(featureImageFloat[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy]-minResponseImag)/rangeResponseImag);
          featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy] = (int)(255.0*(featureImageFloat[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy]-minResponseReal)/rangeResponseReal);
          }
        
        /*
        if (featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy]<0)
            featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy]=0;
        if (featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy]>255)
            featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2][zz][xx][yy]=255;
        if (featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy]<0)
            featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy]=0;
        if (featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy]>255)
            featureImage[(scaleIndex*OrientationNUM+orientationIndex)*2+1][zz][xx][yy]=255;
        */
        }
}



void RegionGrowingInBinaryImage(unsigned char ***mask, Ivector3d imageSize, int RegionPointNumThre)
{
  unsigned char ***status;
  int  *PointNumInCertainOrder ;
  int  ***Detected_Region, order ;
  unsigned char ***VoxelInStack ;
  long int PointNum;
  int i,j,k;
  Ivector3d *Stack;
  
  status = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);
  Detected_Region   = Ialloc3d(imageSize.x, imageSize.y, imageSize.z);
  PointNumInCertainOrder = Ialloc1d(imageSize.x*imageSize.y*imageSize.z) ;
  VoxelInStack = UCalloc3d(imageSize.x, imageSize.y, imageSize.z);
  Stack = Ivector3dalloc1d(imageSize.x*imageSize.y*imageSize.z) ;
  
  for( i=0; i<imageSize.x; i++ )
    for( j=0; j<imageSize.y; j++ )
      for( k=0; k<imageSize.z; k++ )
        status[k][i][j] = NNO ;
        
        
  printf("Now, region growing -- only regions with more than %d voxels are kept...\n", RegionPointNumThre);
  order = 1 ;
  for( k=0; k<imageSize.z; k++ )
   for( i=0; i<imageSize.x; i++ )
    for( j=0; j<imageSize.y; j++ )      
    {
      if(mask[k][i][j]>0 && status[k][i][j]==NNO) 
        {
          PointNum = 0 ;

          //printf("(%d,%d,%d) = %d\n", i,j,k, mask[k][i][j]) ;

          NeighboringSimilarPointsSearch( &PointNum, Detected_Region, order, mask, status, VoxelInStack, Stack, i, j, k, 0, imageSize.x, imageSize.y, imageSize.z) ;
          
          
          if ( PointNum>RegionPointNumThre )
            printf("  Keep region #%d, source point=(%d,%d,%d), PointNum=%ld\n", order, i,j,k, PointNum) ;

          PointNumInCertainOrder[order] = PointNum ;
          order ++ ;
        }
    }

  // remove small connected regions 
  for( i=0; i<imageSize.x; i++ )
    for( j=0; j<imageSize.y; j++ )
      for( k=0; k<imageSize.z; k++ )
    {
      order = Detected_Region[k][i][j] ;
      if( order>=1 && PointNumInCertainOrder[order]<RegionPointNumThre )
        mask[k][i][j] = 0 ;
    }
  printf("done!\n");
  
  UCfree3d(status, imageSize.z, imageSize.x);
  UCfree3d(VoxelInStack, imageSize.z, imageSize.x);
  Ifree3d(Detected_Region, imageSize.z, imageSize.x);

}



void NeighboringSimilarPointsSearch( long int *PointNum, int ***Detected_Region, int order, unsigned char ***src, unsigned char ***status, unsigned char ***VoxelInStack, Ivector3d *Stack, int i, int j, int k, int Greythreshold, int x_size, int y_size, int z_size ) 
{
  int kk,ii,jj; 
  int k_plus_nbrSize, k_minus_nbrSize, i_plus_nbrSize, i_minus_nbrSize, j_plus_nbrSize, j_minus_nbrSize, x, y;
  long int   Pointer ;
  int nbrSize=1;
  //int considerNeighborSlice=YYES;
  int considerNeighborSlice=NNO;
  
  
  /* add the first searched point */
  Pointer = 0 ; 
  Stack[Pointer].x = i; 
  Stack[Pointer].y = j; 
  Stack[Pointer].z = k; 


  do{
    i = Stack[Pointer].x ;
    j = Stack[Pointer].y ;
    k = Stack[Pointer].z ;
    Pointer-- ;
    status[k][i][j] = YYES ;

    (*PointNum) ++ ;   
    Detected_Region[k][i][j] = order ;

    k_plus_nbrSize =k+nbrSize;
    k_minus_nbrSize=k-nbrSize;
    i_plus_nbrSize =i+nbrSize;
    i_minus_nbrSize=i-nbrSize;
    j_plus_nbrSize =j+nbrSize;
    j_minus_nbrSize=j-nbrSize;

    if (k_plus_nbrSize>=z_size)     k_plus_nbrSize=z_size-1;
    if (k_minus_nbrSize<0)          k_minus_nbrSize=0;
    if (j_plus_nbrSize>=y_size)     j_plus_nbrSize=y_size-1;
    if (j_minus_nbrSize<0)          j_minus_nbrSize=0;
    if (i_plus_nbrSize>=x_size)     i_plus_nbrSize=x_size-1;
    if (i_minus_nbrSize<0)          i_minus_nbrSize=0;

    // now we can check immediately adjacent points to see if they too could be added to the track 
    if (considerNeighborSlice)
      {   
      for(kk=k_minus_nbrSize; kk<=k_plus_nbrSize; kk++) 
       for(ii=i_minus_nbrSize; ii<=i_plus_nbrSize; ii++)
    for(jj=j_minus_nbrSize; jj<=j_plus_nbrSize; jj++)
        {
        if( !(kk==k && jj==j && ii==i) && src[kk][ii][jj]>Greythreshold && status[kk][ii][jj]==NNO && VoxelInStack[kk][ii][jj]==NNO )
          {
        Pointer ++ ;
        Stack[Pointer].x = ii ;
        Stack[Pointer].y = jj; 
        Stack[Pointer].z = kk; 

        VoxelInStack[kk][ii][jj]= YYES ;
          }
        }
      }
    else
      {
        kk=k;
        for(ii=i_minus_nbrSize; ii<=i_plus_nbrSize; ii++)
          for(jj=j_minus_nbrSize; jj<=j_plus_nbrSize; jj++)
          {
          if( !(kk==k && jj==j && ii==i) && src[kk][ii][jj]>Greythreshold && status[kk][ii][jj]==NNO && VoxelInStack[kk][ii][jj]==NNO )
          {
          Pointer ++ ;
          Stack[Pointer].x = ii ;
          Stack[Pointer].y = jj; 
          Stack[Pointer].z = kk; 
          
          VoxelInStack[kk][ii][jj]= YYES ;
          }
          }
      }
  }while(Pointer>=0) ;  
}




