/**
 * @file  MatchHistograms.cxx
 * @brief Checks histograms of two input images and matches them if necessary.
 *
 * Copyright (c) 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */
 
#include <iostream>

#include <common/imageio.h> 
#include <common/general.h>
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
    cout << "This program checks histogram of two input byte images and match them if necessary (always match from histogram of low entrophy to high entropy)" << endl;
    cout << "The output are two images with the exact same names of input and overwrites the input." << endl << endl;
    cout << "Please note: we assume byte input in this program. You need to convert other datatype into byte before using this program." << endl;
    cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <inputByteImg1> <inputByteImg2>" << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
    cout << "  <inputByteImg1>        : first input byte image." << endl;
    cout << "  <inputByteImg2>        : second input byte image." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
    cout << "  -t <int>               : threshold, to ignore the intesities less than this; default: 0" << endl;
    cout << "  -S                     : smooth histogram before matching. (default: off)" << endl;
    cout << "  -C                     : linear scale intensity to [0,255] after histogram matching (default: off)" << endl;
    cout << "  -M                     : force intensity matching (default: off)" << endl;
    cout << endl;
    cout << "Example:" << endl;
    cout << "  " << exec_name << " A.img  B.img" << endl;
    cout << endl;
    print_contact();
}

// ===========================================================================
// auxiliary functions
// ===========================================================================

// ---------------------------------------------------------------------------
void Vector_Normalization(float *vec, int len)
{
  int        i;
  float      total ;
  
  /* normalization ... */
  total = 0 ;
  for (i=0; i<len; i++)  
    total += vec[i] ; 

  for (i=0; i<len; i++)  
    vec[i] = vec[i]/total ;
}

// ---------------------------------------------------------------------------
void smoothHistogram(float *hist)
{
  double* gaussianFilter1D;
  int filterWidth=7;
  double smoothSigma=1.0;
  int halfWidth;
  int i, ii, intensity;
  float weight, sum;
  
  gaussianFilter1D=Dalloc1d(filterWidth);
  halfWidth=filterWidth/2;

  for (i=-halfWidth;i<=halfWidth;i++)
    gaussianFilter1D[i+halfWidth]=exp((-i*i)/(2*smoothSigma*smoothSigma))/(smoothSigma*sqrt(2*E_PI));

  
  float *histTemp;
  histTemp = (float *)malloc((MaxOfUC+1)*sizeof(float));
  for (i=0;i<(MaxOfUC+1);i++)
    histTemp[i]=hist[i];
    
  
  for (intensity=0;intensity<=(MaxOfUC+1);intensity++)
  {  
  sum=0.0;
  weight=0.0;
  hist[intensity]=0.0;
  
  for (ii=-halfWidth;ii<=halfWidth;ii++)
    {
      if ( (intensity+ii>=0) && (intensity+ii<(MaxOfUC+1)) )
        {
        sum += histTemp[intensity+ii]*gaussianFilter1D[ii+halfWidth];
        weight += gaussianFilter1D[ii+halfWidth];
        }
    }  
  
  if (weight!=0)   hist[intensity]=sum/weight;   
  }
    
  free(histTemp);
}

// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc,char *argv[])
{
    unsigned char ***imgA, ***imgB;
    int           i, j, k, t;
    float         *histoA, *histoB, *histoB_trans, *histoA_trans;
    float         *cumuHistA, *cumuHistB;
    float         *lookupTable;
    double        TotalIntensityA, TotalIntensityB, NonZeroPntsA, NonZeroPntsB ;
    int           sameImageSize=NNO;
    float         entropyA, entropyB;
    int           intA, intB;

    int threshold = 0 ;
    int smoothOrNot = NNO;
    int scalingIntensityOrNot = NNO;
    int forceMatching = NNO;

    // show usage if no arguments were given
    if (argc == 1) {
        print_help();
        exit(1);
    }

    // parse arguments
    int c = -1;
    while((c = getopt(argc, argv, "t:SCMh")) != -1) {
        switch(c) {
            case 't':
                sscanf(optarg, "%d", &threshold);
            break;

            case 'S':
                smoothOrNot = YYES;
                break;

            case 'C':
                scalingIntensityOrNot = YYES;
                break;

            case 'M':
                forceMatching = YYES;
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

  


  const char* imageAfilename   = argv[0];
  const char* imageBfilename   = argv[1];
  
  // read input images
  Image* imageA = ReadImage(imageAfilename);
  Image* imageB = ReadImage(imageBfilename);
  imgA = imageA->img.uc;
  imgB = imageB->img.uc;
    
    
  // see if their image dimensions are the same (image dimensions can be same or different)
  if ( (imageA->region.nx==imageB->region.nx)&&(imageA->region.ny==imageB->region.ny)&&(imageA->region.nz==imageB->region.nz) )
    sameImageSize=YYES;
    
  
  // create histogram 
  printf("\nStart to estimate histogram...\n");
  histoA = (float *)malloc((MaxOfUC+1)*sizeof(float));
  histoB = (float *)malloc((MaxOfUC+1)*sizeof(float));
  cumuHistA = (float *)malloc((MaxOfUC+1)*sizeof(float));
  cumuHistB = (float *)malloc((MaxOfUC+1)*sizeof(float));
  for (i=0; i<(MaxOfUC+1); i++) 
    {
    histoA[i] = 0.0;
    histoB[i] = 0.0;
    cumuHistA[i] = 0.0;
    cumuHistB[i] = 0.0;
    }

    
  TotalIntensityA = 0 ;
  TotalIntensityB = 0 ;
  NonZeroPntsA = 0 ;
  NonZeroPntsB = 0 ;
  if (sameImageSize)
  {
  for (k=0; k<imageA->region.nz; k++)
    for (i=0; i<imageA->region.nx; i++)
      for (j=0; j<imageA->region.ny; j++) 
        { 
        if (imgA[k][i][j]<threshold)
            imgA[k][i][j]=0;
        if (imgB[k][i][j]<threshold)
            imgB[k][i][j]=0;
        
        histoA[imgA[k][i][j]] +=1.0; 
        histoB[imgB[k][i][j]] +=1.0; 
        TotalIntensityA += imgA[k][i][j];       
        TotalIntensityB += imgB[k][i][j]; 
        
        if (imgA[k][i][j]>0)
            NonZeroPntsA += 1.0;
        if (imgB[k][i][j]>0)
            NonZeroPntsB += 1.0;
        }
  }
  else
  {
  for (k=0; k<imageA->region.nz; k++)
    for (i=0; i<imageA->region.nx; i++)
      for (j=0; j<imageA->region.ny; j++) 
        { 
        if (imgA[k][i][j]<threshold)
            imgA[k][i][j]=0;
        
        histoA[imgA[k][i][j]] +=1.0; 
        TotalIntensityA += imgA[k][i][j]; 
        
        if (imgA[k][i][j]>0)
            NonZeroPntsA += 1.0;
        }
  for (k=0; k<imageB->region.nz; k++)
    for (i=0; i<imageB->region.nx; i++)
      for (j=0; j<imageB->region.ny; j++) 
        {
        if (imgB[k][i][j]<threshold)
            imgB[k][i][j]=0;
        
        histoB[imgB[k][i][j]] +=1.0; 
        TotalIntensityB += imgB[k][i][j]; 
        if (imgB[k][i][j]>0)
            NonZeroPntsB += 1.0;
        }
  }     
        
  
  // derive histogram
  histoA[0]=0.0;    
  histoB[0]=0.0;    
  cumuHistA[0]=0.0;
  cumuHistB[0]=0.0;
  entropyA=0.0;
  entropyB=0.0;
  for (i=1; i<(MaxOfUC+1); i++) 
    { 
    histoA[i] /= NonZeroPntsA; 
    cumuHistA[i] = cumuHistA[i-1] + histoA[i];
    if (histoA[i]!=0)
        entropyA -= (histoA[i]*log(histoA[i]));
        
    histoB[i] /= NonZeroPntsB; 
    cumuHistB[i] = cumuHistB[i-1] + histoB[i];
    if (histoB[i]!=0)
        entropyB -= (histoB[i]*log(histoB[i]));
    }
  printf("TotalIntensityA=%f   NonZeroPntsA=%f\n", TotalIntensityA, NonZeroPntsA) ;
  printf("TotalIntensityB=%f   NonZeroPntsB=%f\n", TotalIntensityB, NonZeroPntsB) ;
  printf("cumulative histogram for A = %f, for B = %f\n", cumuHistA[MaxOfUC], cumuHistB[MaxOfUC]);
  printf("\nEntropy A = %f\nEntropy B = %f\n", entropyA, entropyB);

  
  // Ou added on 09/06/2010 (begin)
  // check if histogram matching is needed
  int int1 = (int)((float)MaxOfUC/5.0 + 0.5);
  int int2 = (int)((float)MaxOfUC*4.0/5.0 + 0.5);
  int MatchingAtoB = NNO;
  int MatchingBtoA = NNO;
  
  
  printf("cumulative histogram for image A: \n\t%f at intensity  12 (perhaps noise), \n\t%f at intensity  %d (lower  1/5 intensity range), \n\t%f at intensity %d (higher 1/5 intensity range).\n", cumuHistA[12], cumuHistA[int1], int1, cumuHistA[int2], int2);
  printf("cumulative histogram for image B: \n\t%f at intensity  12 (perhaps noise), \n\t%f at intensity  %d (lower  1/5 intensity range), \n\t%f at intensity %d (higher 1/5 intensity range).\n", cumuHistB[12], cumuHistB[int1], int1, cumuHistB[int2], int2);

  if (forceMatching==YYES)  // match histogram if forced to do so
    {
    if (entropyA>entropyB)  MatchingBtoA=YYES;
    if (entropyA<entropyB)  MatchingAtoB=YYES;
    }
  else if( ( (cumuHistA[int1]>0.8)||(cumuHistA[int1]<0.04)||(cumuHistA[int2]<0.2)||(cumuHistA[int2]>0.96)||(cumuHistB[int1]>0.8)||(cumuHistB[int1]<0.04)||(cumuHistB[int2]<0.2)||(cumuHistB[int2]>0.96) ) && (cumuHistA[12]<0.20 && cumuHistB[12]<0.20 ) )  // match histogram if 1) screwed histogram in either image 2) both images are not occurpied by noise.
    {
    if (entropyA>entropyB)  MatchingBtoA=YYES;
    if (entropyA<entropyB)  MatchingAtoB=YYES;
    }
  else
    {
    printf("No need to match histogram. \n\n");
    free(histoA);
    free(histoB);
    delete imageA;
    delete imageB;
    exit(0);
    }
  
  if (smoothOrNot)
  {
  printf("Histogram smoothing ...\n");
  smoothHistogram(histoA);
  smoothHistogram(histoB);

  printf("Histogram normalization ... \n");
  Vector_Normalization(histoA, (MaxOfUC+1)) ;
  Vector_Normalization(histoB, (MaxOfUC+1)) ;
  
  printf("Histogram cumulation...\n");
  for (i=1; i<(MaxOfUC+1); i++) 
    {
    cumuHistA[i] = histoA[i]+cumuHistA[i-1];
    cumuHistB[i] = histoB[i]+cumuHistB[i-1];
    }
  }
  
  
  // start histogram matching
  lookupTable = (float *)malloc((MaxOfUC+1)*sizeof(float));
  if (MatchingAtoB)
  {
  // matching imgA to imgB 
  printf("-----------------------------------\n");
  printf("Matching histogram A to histogram B\n");
  printf("-----------------------------------\n");
  
  // Yangming rewrote the histogram matching part here on 09/07/2010(begin)
  intB = 0;
  lookupTable[0]=0;
  for (intA=1;intA<(MaxOfUC+1);intA++)
    {
    if (cumuHistA[intA]<=cumuHistB[intB])
        lookupTable[intA] = intB;
    else
        {
        while ( (cumuHistA[intA]>cumuHistB[intB])&&(intB<MaxOfUC) )     intB++;
        if ( (cumuHistB[intB]-cumuHistA[intA])>(cumuHistA[intA]-cumuHistB[intB-1]) )
            lookupTable[intA] = intB--;
        else
            lookupTable[intA] = intB;
        }
    //printf("lookup table %d in A -> %d in B\n", intA, intB);
    }
    
  if (scalingIntensityOrNot==YYES)
    {
    int minIntA=static_cast<int>(lookupTable[0]);
    int maxIntA=static_cast<int>(lookupTable[MaxOfUC]);
    float rangeIntA = (float)(maxIntA-minIntA);
    for (intA=1;intA<(MaxOfUC+1);intA++)
      lookupTable[intA] = (int)(255.0*(lookupTable[intA]-minIntA)/rangeIntA+0.5);
    }
    
    
  TotalIntensityA = 0 ;
  for (k=0; k<imageA->region.nz; k++)
    for (i=0; i<imageA->region.nx; i++)
      for (j=0; j<imageA->region.ny; j++) 
        {
        imgA[k][i][j] = static_cast<int>(lookupTable[imgA[k][i][j]]);
        TotalIntensityA += imgA[k][i][j] ;
        }
  printf("new TotalIntensityA=%f\n", TotalIntensityA) ;
  // Yangming rewrote the histogram matching part here on 09/07/2010(end)
  
  
  // save result 
  WriteImage(imageAfilename, imageA);
  
  printf("\nOutput:\nimage %s has been re-written as a result of matching its histogram to %s\n\n", argv[0], argv[1]);
  }
  
  
  if (MatchingBtoA)
  {
  // matching imgB to imgA
  printf("-----------------------------------\n");
  printf("Matching histogram B to histogram A\n");
  printf("-----------------------------------\n");
  
  
  // Yangming rewrite the histogram matching part here on 09/07/2010(begin)
  intA = 0;
  lookupTable[0]=0;
  for (intB=1;intB<(MaxOfUC+1);intB++)
    {
    if (cumuHistB[intB]<=cumuHistA[intA])
        lookupTable[intB] = intA;
    else
        {
        while ( (cumuHistB[intB]>cumuHistA[intA])&(intA<MaxOfUC) )  intA++;
        if (cumuHistA[intA]-cumuHistB[intB]>cumuHistB[intB]-cumuHistA[intA-1])
            lookupTable[intB] = intA--;
        else
            lookupTable[intB] = intA;
        }
    
    //printf("%d in B with cdf %f => %d in A with cdf %f\n", intB, cumuHistB[intB], intA, cumuHistA[intA]);
    }
    
    
  if (scalingIntensityOrNot==YYES)
    {
    int minIntB=static_cast<int>(lookupTable[0]);
    int maxIntB=static_cast<int>(lookupTable[MaxOfUC]);
    float rangeIntB = (float)(maxIntB-minIntB);
    for (intB=1;intB<(MaxOfUC+1);intB++)
      lookupTable[intB] = (int)(255.0*(lookupTable[intB]-minIntB)/rangeIntB+0.5);
    }
    
    
  TotalIntensityB = 0 ;
  for (k=0; k<imageB->region.nz; k++)
    for (i=0; i<imageB->region.nx; i++)
      for (j=0; j<imageB->region.ny; j++) 
        {
        imgB[k][i][j] = static_cast<int>(lookupTable[imgB[k][i][j]]);
        TotalIntensityB += imgB[k][i][j] ;
        }
  printf("new TotalIntensityB=%f\n", TotalIntensityB) ;
  // Yangming rewrite the histogram matching part here on 09/07/2010(end)
  
  
  // save result 
  WriteImage(imageBfilename, imageB);
  
  printf("\nOutput:\nimage %s has been re-written as a result of matching its histogram to %s\n\n", argv[1], argv[0]);
  }

  free(histoA);
  free(histoB);
  delete imageA;
  delete imageB;

  exit(0);
}
