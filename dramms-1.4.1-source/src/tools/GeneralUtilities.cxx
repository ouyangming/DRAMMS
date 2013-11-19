#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include "common/image.h"

#include "GeneralUtilities.h"


#define MyMIN(x,y)  ((x) < (y) ? (x) : (y))
#define MyMAX(x,y)  ((x) > (y) ? (x) : (y))
#define SSD  0
#define CC   1

int myround(float x)
{
    int rounded = (int)(x+0.5);
    
	return rounded;
}

/*
float myTrilinearInterpolation(unsigned char ***Img, float x, float y, float z, Ivector3d imageSize)
{
	float intensity=0.0;
	int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
	
	xFloor = (int)MyMAX( MyMIN( floor(x), (float)(imageSize.x-1) ), 0.0);
	xCeil = (int)MyMAX( MyMIN( ceil(x), (float)(imageSize.x-1) ), 0.0);
	yFloor = (int)MyMAX( MyMIN( floor(y), (float)(imageSize.y-1) ), 0.0);
	yCeil = (int)MyMAX( MyMIN( ceil(y), (float)(imageSize.y-1) ), 0.0);
	zFloor = (int)MyMAX( MyMIN( floor(z), (float)(imageSize.z-1) ), 0.0);
	zCeil = (int)MyMAX( MyMIN( ceil(z), (float)(imageSize.z-1) ), 0.0);
		   
	if ( (xFloor==xCeil)&(xFloor!=(imageSize.x-1)) )  xCeil += 1;
	if ( (xFloor==xCeil)&(xFloor==(imageSize.x-1)) )  xFloor -= 1;
	if ( (yFloor==yCeil)&(yFloor!=(imageSize.y-1)) )  yCeil += 1;
	if ( (yFloor==yCeil)&(yFloor==(imageSize.y-1)) )  yFloor -= 1;
	if ( (zFloor==zCeil)&(zFloor!=(imageSize.z-1)) )  zCeil += 1;
	if ( (zFloor==zCeil)&(zFloor==(imageSize.z-1)) )  zFloor -= 1;
	
		
	intensity = (x-(float)xFloor)*(y-(float)yFloor)*(z-(float)zFloor)*(float)Img[zCeil][xCeil][yCeil] + 
				(x-(float)xFloor)*(y-(float)yFloor)*((float)zCeil-z)*(float)Img[MyMAX(zFloor,0)][xCeil][yCeil] + 
				(x-(float)xFloor)*((float)yCeil-y)*(z-(float)zFloor)*(float)Img[zCeil][xCeil][MyMAX(yFloor,0)] + 
				(x-(float)xFloor)*((float)yCeil-y)*((float)zCeil-z)*(float)Img[MyMAX(zFloor,0)][xCeil][MyMAX(yFloor,0)] + 
				((float)xCeil-x)*(y-(float)yFloor)*(z-(float)zFloor)*(float)Img[zCeil][MyMAX(xFloor,0)][yCeil] + 
				((float)xCeil-x)*(y-(float)yFloor)*((float)zCeil-z)*(float)Img[MyMAX(zFloor,0)][MyMAX(xFloor,0)][yCeil] + 
				((float)xCeil-x)*((float)yCeil-y)*(z-(float)zFloor)*(float)Img[zCeil][MyMAX(xFloor,0)][MyMAX(yFloor,0)] + 
				((float)xCeil-x)*((float)yCeil-y)*((float)zCeil-z)*(float)Img[MyMAX(zFloor,0)][MyMAX(xFloor,0)][MyMAX(yFloor,0)];
	  
	
    //if (intensity<0.0) intensity=0.0;
	
	return intensity;
}
*/


float featureCost(unsigned char ****TF, unsigned char ****SF, float ***confidenceMap, int x_trg, int y_trg,
				  int z_trg, float x_src, float y_src, float z_src, int numFeatures, int SimilarityMeasure, Ivector3d imageSize)
{
    //here you should put your own function to calculate the cost based on
    //the optimal features.
	int featureIndex;
	float energySquaredAtThisPoint = 0.0f;
	float cc=0.0;
	float totalEnergy = 0.0f;
	float feature;

	
	// trilinear interpolation
	int xFloor, xCeil, yFloor, yCeil, zFloor, zCeil;
	
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
	
	float distXtoXFloor = x_src-xFloor;
	float distXtoXCeil  = xCeil-x_src;
	float distYtoYFloor = y_src-yFloor;
	float distYtoYCeil  = yCeil-y_src;
	float distZtoZFloor = z_src-zFloor;
	float distZtoZCeil  = zCeil-z_src;
	
	float w1 = distXtoXFloor*distYtoYFloor*distZtoZFloor;
	float w2 = distXtoXFloor*distYtoYFloor*distZtoZCeil;
	float w3 = distXtoXFloor*distYtoYCeil*distZtoZFloor;
	float w4 = distXtoXFloor*distYtoYCeil*distZtoZCeil;
	float w5 = distXtoXCeil*distYtoYFloor*distZtoZFloor;
	float w6 = distXtoXCeil*distYtoYFloor*distZtoZCeil;
	float w7 = distXtoXCeil*distYtoYCeil*distZtoZFloor;
	float w8 = distXtoXCeil*distYtoYCeil*distZtoZCeil;
	
	
	
	int judgement=0;
	if (x_src<0 || x_src>(imageSize.x-1) || y_src<0 || y_src>(imageSize.y-1) || z_src<0 || z_src>(imageSize.z-1))
		judgement=1;
	
		
	if (SimilarityMeasure==SSD)
	{
	//#pragma omp parallel for num_threads(numFeatures) private(featureIndex,feature) reduction(+:energySquaredAtThisPoint) 
	for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
		  {
		  //printf("thread %d\t", omp_get_thread_num());
		  //feature = trilinearInterpolation(SF[featureIndex], x_src, y_src, z_src, imageSize);
		  feature = SF[featureIndex][zCeil][xCeil][yCeil]   * w1 +
					SF[featureIndex][zFloor][xCeil][yCeil]  * w2 +
					SF[featureIndex][zCeil][xCeil][yFloor]  * w3 +
					SF[featureIndex][zFloor][xCeil][yFloor] * w4 +
					SF[featureIndex][zCeil][xFloor][yCeil]  * w5 +
					SF[featureIndex][zFloor][xFloor][yCeil] * w6 +
					SF[featureIndex][zCeil][xFloor][yFloor] * w7 +
					SF[featureIndex][zFloor][xFloor][yFloor]* w8 ;
		  
		  //if ( (feature==0.0)&&(x_src<0 || x_src>=(imageSize.x-1) || y_src<0 || y_src>=(imageSize.y-1) || z_src<0 || z_src>=(imageSize.z-1)) )
		  //if ( (feature==0.0)&&(judgement==1) )
		    //feature = MyMAX(0.0f, SF[featureIndex][z_trg][x_trg][y_trg]);
		  //I AM NOT SO SURE IF THE PREVIOUS EXPRESSION SHOULD BE LIKE THAT
			
          //if (featureIndex%2==0)  // imaginary part			
			//energySquaredAtThisPoint += pow( (fabs(feature-128.0) - fabs((float)TF[featureIndex][z_trg][x_trg][y_trg] - 128.0)), 2.0f );
		  //else // real part
			energySquaredAtThisPoint += pow( (feature - (float)TF[featureIndex][z_trg][x_trg][y_trg]), 2.0f );
		  }
		  
		  
	//totalEnergy += energySquaredAtThisPoint;

	totalEnergy = energySquaredAtThisPoint * confidenceMap[z_trg][x_trg][y_trg] / (float)numFeatures;   // confidence map involved
	}
	
	
	if (SimilarityMeasure==CC)
	{
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
	float featureA, featureB;
   
	//#pragma omp parallel for num_threads(numFeatures) private(featureIndex,feature) reduction(+:energySquaredAtThisPoint) 
	for (featureIndex=0; featureIndex<numFeatures; featureIndex++)
		  {
		  //printf("thread %d\t", omp_get_thread_num());
		  //feature = trilinearInterpolation(SF[featureIndex], x_src, y_src, z_src, imageSize);
		  featureA= SF[featureIndex][zCeil][xCeil][yCeil]   * w1 +
					SF[featureIndex][zFloor][xCeil][yCeil]  * w2 +
					SF[featureIndex][zCeil][xCeil][yFloor]  * w3 +
					SF[featureIndex][zFloor][xCeil][yFloor] * w4 +
					SF[featureIndex][zCeil][xFloor][yCeil]  * w5 +
					SF[featureIndex][zFloor][xFloor][yCeil] * w6 +
					SF[featureIndex][zCeil][xFloor][yFloor] * w7 +
					SF[featureIndex][zFloor][xFloor][yFloor]* w8 ;
		  featureB=(float)TF[featureIndex][z_trg][x_trg][y_trg];
		  
		  sumA += featureA;
		  sumB += featureB;
		  sumAB += featureA*featureB;
		  sumAsq += featureA*featureA;
		  sumBsq += featureB*featureB;
		  }
	
	aveA = sumA/(float)numFeatures;
	aveB = sumB/(float)numFeatures;
   
	covAB = sumAB-aveB*sumA-aveA*sumB+(float)numFeatures*aveA*aveB;
	varA  = sumAsq-(float)numFeatures*aveA*aveA;
	varB  = sumBsq-(float)numFeatures*aveB*aveB;
   
    if (varA==0 && varB==0)
		cc=1.0;
	else if (varA==0 || varB==0)
		cc=0.0;
	else
		cc = fabsf(covAB/sqrt(varA*varB));
		//cc = fabs(covAB/sqrt(varA*varB));

	totalEnergy = -15000 * cc * confidenceMap[z_trg][x_trg][y_trg];   // confidence map involved
	
	}
	
	return totalEnergy;
}

void precomputeWeights(Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d resolutionRatio, int method,
					   float *weightsX, float *weightsY, float *weightsZ)
{
    //function that precomputes the spatial varying weights for the linear and cubic b-splines case. These weights are
	//used for the evaluation of the singleton potentials

    // X,Y coordinates of current pixel
    int x, y, z;
	int indexFirstEffectiveControlPointX, indexFirstEffectiveControlPointY, indexFirstEffectiveControlPointZ;
	int indexEffectiveControlPointX, indexEffectiveControlPointY, indexEffectiveControlPointZ;
    
    //help index
    int h, hI;
    
    // weighting variables
    float u, v, w;
	float *Bu, *Bv, *Bw;
	float sumWeight;

	int numControlPointsX, numControlPointsY, numControlPointsZ;

	// allocate memory
	if(method == 1) // linear weighting
	{
		Bu = (float *)calloc(2, sizeof(float));
		Bv = (float *)calloc(2, sizeof(float));
		Bw = (float *)calloc(2, sizeof(float));
	}
	else if(method == 2)  // B-spline coefficient
	{
		Bu = (float *)calloc(4, sizeof(float));
		Bv = (float *)calloc(4, sizeof(float));
		Bw = (float *)calloc(4, sizeof(float));
	}                
    
	numControlPointsX = (int)ceil((float)imageSize.x/(float)distBetweenControlPointsX);
	numControlPointsY = (int)ceil((float)imageSize.y/(float)distBetweenControlPointsY);
	numControlPointsZ = (int)ceil((float)(imageSize.z-1)/(float)distBetweenControlPointsZ);


    if(method == 2)//cubic b-splines case
    {                      
        for (x=0; x<imageSize.x; x++)
        {
            u = ((float)x/(float)distBetweenControlPointsX) - floor((float)x/(float)distBetweenControlPointsX);
		    indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPointsX) - 1;
		
			// Yangming added in each weight an if/else switch on August 16, 2010
		    indexEffectiveControlPointX = indexFirstEffectiveControlPointX;
			if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX) )
				Bu[0] = ((1-u)*(1-u)*(1-u))/6;   
			else
				Bu[0] = 0.0f;
				
			indexEffectiveControlPointX++;
			if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX) )
				Bu[1] = (((3*u-6)*u)*u+4)/6;
			else 
				Bu[1] = 0.0f;
				
			indexEffectiveControlPointX++;
			if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX) )
				Bu[2] = ((((-3*u+3)*u+3)*u)+1)/6;
			else
				Bu[2] = 0.0f;
			
			indexEffectiveControlPointX++;
			if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX) )
				Bu[3] = (u*u*u)/6;
			else 
				Bu[3] = 0.0f;
			
			
			// Yangming changed on August 13 2010 to accomodate anisotropic voxel size			
			Bu[0] = Bu[0]/resolutionRatio.x;   //added on August 13, 2010 by Yangming
			Bu[3] = Bu[3]/resolutionRatio.x;   //added on August 13, 2010 by Yangming
			sumWeight = Bu[0]+Bu[1]+Bu[2]+Bu[3];
			Bu[0] = Bu[0]/sumWeight;
			Bu[1] = Bu[1]/sumWeight;
			Bu[2] = Bu[2]/sumWeight;
			Bu[3] = Bu[3]/sumWeight;
			
			
            for (h=0; h<4; h++){
                hI = (int)(h*imageSize.x+x);
                weightsX[hI] = Bu[h];
            }
        }
        
        for (y=0; y<imageSize.y; y++)
        {
            v = ((float)y/(float)distBetweenControlPointsY) - floor((float)y/(float)distBetweenControlPointsY);
			indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPointsY) - 1;
			
			// Yangming added in each weight an if/else switch on August 16, 2010
		    indexEffectiveControlPointY = indexFirstEffectiveControlPointY;
			if ( (indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY) )
				Bv[0] = ((1-v)*(1-v)*(1-v))/6;
			else
				Bv[0] = 0.0f;
				
			indexEffectiveControlPointY++;
			if ( (indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY) )
				Bv[1] = (((3*v-6)*v)*v+4)/6;
			else
				Bv[1] = 0.0f;
				
			indexEffectiveControlPointY++;
			if ( (indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY) )				
				Bv[2] = ((((-3*v+3)*v+3)*v)+1)/6;
			else 
				Bv[2] = 0.0f;
				
			indexEffectiveControlPointY++;	
			if ( (indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY) )				
				Bv[3] = (v*v*v)/6;
			else
				Bv[3] = 0.0f;
			
			// Yangming changed on August 13 2010 to accomodate anisotropic voxel size			
			Bv[0] = Bv[0]/resolutionRatio.y;   //added on August 13, 2010 by Yangming
			Bv[3] = Bv[3]/resolutionRatio.y;   //added on August 13, 2010 by Yangming
			sumWeight = Bv[0]+Bv[1]+Bv[2]+Bv[3];
			Bv[0] = Bv[0]/sumWeight;
			Bv[1] = Bv[1]/sumWeight;
			Bv[2] = Bv[2]/sumWeight;
			Bv[3] = Bv[3]/sumWeight;
			
			
            for (h=0; h<4; h++)
            {
                hI = (int)(h*imageSize.y+y);
                weightsY[hI] = Bv[h];
            }
        }
        
        for (z=0; z<imageSize.z; z++)
        {
            w = ((float)(z-1)/(float)distBetweenControlPointsZ) - floor((float)(z-1)/(float)distBetweenControlPointsZ);
			indexFirstEffectiveControlPointZ = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ) - 1;
			
			// Yangming added in each weight an if/else switch on August 16, 2010
		    indexEffectiveControlPointZ = indexFirstEffectiveControlPointZ;
			if ( (indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				Bw[0] = ((1-w)*(1-w)*(1-w))/6;
			else
				Bw[0] = 0.0f;
				
			indexEffectiveControlPointZ++;
			if ( (indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				Bw[1] = (((3*w-6)*w)*w+4)/6;
			else
				Bw[1] = 0.0f;
				
			indexEffectiveControlPointZ++;
			if ( (indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				Bw[2] = ((((-3*w+3)*w+3)*w)+1)/6;
			else
				Bw[2] = 0.0f;
				
			indexEffectiveControlPointZ++;
			if ( (indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				Bw[3] = (w*w*w)/6;
			else
				Bw[3] = 0.0f;
		
		
			// Yangming changed on August 13 2010 to accomodate anisotropic voxel size			
			Bw[0] = Bw[0]/resolutionRatio.z;   //added on August 13, 2010 by Yangming
			Bw[3] = Bw[3]/resolutionRatio.z;   //added on August 13, 2010 by Yangming
			sumWeight = Bw[0]+Bw[1]+Bw[2]+Bw[3];
			Bw[0] = Bw[0]/sumWeight;
			Bw[1] = Bw[1]/sumWeight;
			Bw[2] = Bw[2]/sumWeight;
			Bw[3] = Bw[3]/sumWeight;
			
			
            for (h=0; h<4; h++)
            {
                hI = (int)(h*imageSize.z+z);
                weightsZ[hI] = Bw[h];
            }			
        }
    }
    else if(method == 1)//trilinear case
    {
        for (x=0; x<imageSize.x; x++)
        {
            u = ((float)x/(float)distBetweenControlPointsX) - floor((float)x/(float)distBetweenControlPointsX);
			indexFirstEffectiveControlPointX = (int)floor((float)x/(float)distBetweenControlPointsX);
			
			indexEffectiveControlPointX = indexFirstEffectiveControlPointX;
			if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX) )
				Bu[0] = 1-u;
			else
				Bu[0] = 0.0f;
				
			indexEffectiveControlPointX++;
			if ( (indexEffectiveControlPointX>=0)&&(indexEffectiveControlPointX<numControlPointsX) )
				Bu[1] = u;
			else
				Bu[1] = 0.0f;
            
			
			sumWeight = Bu[0]+Bu[1];
			if (sumWeight>0)
				{
				Bu[0] /= sumWeight;
				Bu[1] /= sumWeight;
				}
			else
				{
				if ( x<=(int)(imageSize.x/2) )  Bu[1]=1.0f;
				else  Bu[0]=1.0f;
				}
			
			
            for (h=0; h<2; h++)
            {
                hI = (int)(h*imageSize.x+x);
                weightsX[hI] = Bu[h];
            }
        }
        
        for (y=0; y<imageSize.y; y++)
        {
            v = ((float)y/(float)distBetweenControlPointsY) - floor((float)y/(float)distBetweenControlPointsY);
			indexFirstEffectiveControlPointY = (int)floor((float)y/(float)distBetweenControlPointsY);
			
			indexEffectiveControlPointY = indexFirstEffectiveControlPointY;
			if ( (indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY) )
				Bv[0] = 1-v;
			else 
				Bv[0] = 0.0f;
				
			indexEffectiveControlPointY++;	
			if ( (indexEffectiveControlPointY>=0)&&(indexEffectiveControlPointY<numControlPointsY) )
				Bv[1] = v;
			else
				Bv[1] = 0.0f;
				
            
			sumWeight = Bv[0]+Bv[1];
			if (sumWeight>0)
				{
				Bv[0] /= sumWeight;
				Bv[1] /= sumWeight;
				}
			else
				{
				if ( y<(int)(imageSize.y/2) )   Bv[1]=1.0f;
				else  Bv[0]=1.0f;
				}
			
			
            for (h=0; h<2; h++)
            {
                hI = (int)(h*imageSize.y+y);
                weightsY[hI] = Bv[h];
            }
        }
        
        
        for (z=0; z<imageSize.z; z++)
        {
            w = ((float)(z-1)/(float)distBetweenControlPointsZ) - floor((float)(z-1)/(float)distBetweenControlPointsZ);
			indexFirstEffectiveControlPointZ = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ);
			
			indexEffectiveControlPointZ = indexFirstEffectiveControlPointZ;
			if ( (indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				Bw[0] = 1-w;
			else
				Bw[0] = 0.0f;
				
			indexEffectiveControlPointZ++;	
			if ( (indexEffectiveControlPointZ>=0)&&(indexEffectiveControlPointZ<numControlPointsZ) )
				Bw[1] = w;
			else
				Bw[1] = 0.0f;
            
			
			sumWeight = Bw[0]+Bw[1];
			if (sumWeight>0)
				{
				Bw[0] /= sumWeight;
				Bw[1] /= sumWeight;
				}
			else
				{
				if ( (z-1)<(int)(imageSize.z/2) )  Bw[1]=1.0f;
				else Bw[0]=1.0f;
				}
				
			
            for (h=0; h<2; h++)
            {
                hI = (int)(h*imageSize.z+z);
                weightsZ[hI] = Bw[h];
            }
        }
    }   

}

