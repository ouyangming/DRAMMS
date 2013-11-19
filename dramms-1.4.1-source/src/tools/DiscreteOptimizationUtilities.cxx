/**
 * @file  DiscreteOptimizationUtilities.cxx
 * @brief Discrete optimization calling FastPD.
 *
 * The Copyright<br />
 * Department of Computer Science,<br />
 * University of Crete, Greece
 *
 * Mathématiques Appliquées aux Systèmes (MAS),<br />
 * Ecole Centrale de Paris, France
 *
 * Copyright (c) 2009. All rights reserved.<br />
 * See http://www.csd.uoc.gr/~komod/FastPD/ for more details.
 *
 * Copyright (c) 2011, 2012 University of Pennylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.<br />
 * Modification based on written authorization from FastPD's author.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */



#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <common/mvcd.h>
#include "GeneralUtilities.h"


void computeLabels(int numSamples, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d resolutionRatio, int levelIndex, int method, float *labels, float label_factor, int iter,  int numIterInEachResolution, int indexGridLevel, int numGridLevels, int fastApproximationOrNot)
{
	//function that computes the labels that correspond to the potential
	//deformations that have to be tested

	//numlabels number of labels are sampled along the principal directions and diagonals
	//(18-neighborhood in 3D). The first label correspond to do nothing. Maximum displacement
	//is constrained to 0.4 of the gridspacing in order for the deformation to be
	//diffeomorphic.
	int i; //index
	float max_x, max_y, max_z; //maximum displacement in each axis
	float step_x, step_y, step_z; //step made in each axis
	int numLabels;
	
	float maxRatio;   // maximum displacement (measured as maximum percent of the distance between two adjacent control points) at each control point to secure diffeomorephism in deformation
	float sumDiv=0.0;
	for (i=0; i<numIterInEachResolution; i++)
		sumDiv += pow(label_factor, (float)i);
		

	if (fastApproximationOrNot==0)
		numLabels = 18*numSamples + 1;
	else
		numLabels = 6*numSamples + 1;
		
	
	if (method==1)  // trilinear
	{
	//maxRatio = 0.78/sumDiv;
	//maxRatio = 0.8*pow(1.5, levelIndex-1)/sumDiv;  // tested, bad for brain and cardiac, 1/3/2011
	//maxRatio=0.42 * pow(0.9, 3-levelIndex);
	//maxRatio=0.40 * pow(1.025, 3-levelIndex);
	maxRatio = 0.40; // was 0.39 till 07/18/2011;
	printf("maxRatio = %f for method %d\n", maxRatio, method);
	max_x = maxRatio*((float)distBetweenControlPointsX) * pow(label_factor, iter + (numGridLevels-indexGridLevel)*0. + (indexGridLevel-1)*0.5);
	max_y = maxRatio*((float)distBetweenControlPointsY) * pow(label_factor, iter + (numGridLevels-indexGridLevel)*0. + (indexGridLevel-1)*0.5);
	max_z = maxRatio*((float)distBetweenControlPointsZ) * pow(label_factor, iter + (numGridLevels-indexGridLevel)*0. + (indexGridLevel-1)*0.5);
	printf("max_x,y,z=%f,%f,%f\n", max_x, max_y,max_z);
	}
	else  // cubic B-spline
	{
	//maxRatio = 0.735/sumDiv;
	//maxRatio = 0.76*pow(1.25, levelIndex-1)/sumDiv;
	maxRatio = 0.39;
	//printf("maxRatio = %f for method %d\n", maxRatio, method);
	max_x = maxRatio*((float)distBetweenControlPointsX) * pow(label_factor, iter+indexGridLevel-1);
	max_y = maxRatio*((float)distBetweenControlPointsY) * pow(label_factor, iter+indexGridLevel-1);
	max_z = maxRatio*((float)distBetweenControlPointsZ) * pow(label_factor, iter+indexGridLevel-1);
	}
	
	
	if (numSamples == 0)
	{
		//at least one sample should be taken in each direction
		numSamples = 1;
	}

	step_x = max_x / ((float)numSamples);
	step_y = max_y / ((float)numSamples);
	step_z = max_z / ((float)numSamples);

	//fill the label array
	//first label is do nothing
	labels[0] = 0; labels[numLabels] = 0; labels[2*numLabels] = 0;

	for(i=1; i<=numSamples; i++)
	{
		if (fastApproximationOrNot==0)
		{
		//18-neighborhood is considered	
		// (+x, -x, +y, -y, +z, -z, (+x)(+y), (+x)(-y), (-x)(+y), (-x)(-y), (+x)(+z), (+x)(-z), (-x)(+z), (-x)(-z), (+y)(+z), (+y)(-z), (-y)(+z), (-y)(-z) )
		labels[ i ] = step_x * i;
		labels[ i + numLabels ] = 0.0f;
		labels[ i + 2*numLabels ] = 0.0f;

		labels[ i + numSamples ] = 0.0f;
		labels[ i + numSamples + numLabels ] = step_y * i;
		labels[ i + numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 2*numSamples ] = 0.0f;
		labels[ i + 2*numSamples + numLabels ] = 0.0f;
		labels[ i + 2*numSamples + 2*numLabels ] = step_z * i;

		labels[ i + 3*numSamples ] = -step_x * i;
		labels[ i + 3*numSamples + numLabels ] = 0.0f;
		labels[ i + 3*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 4*numSamples ] = 0.0f;
		labels[ i + 4*numSamples + numLabels ] = -step_y * i;
		labels[ i + 4*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 5*numSamples ] = 0.0f;
		labels[ i + 5*numSamples + numLabels ] = 0.0f;
		labels[ i + 5*numSamples + 2*numLabels ] = -step_z * i;
		
		labels[ i + 6*numSamples ] = step_x * i;
		labels[ i + 6*numSamples + numLabels ] = step_y * i;
		labels[ i + 6*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 7*numSamples ] = step_x * i;
		labels[ i + 7*numSamples + numLabels ] = -step_y * i;
		labels[ i + 7*numSamples + 2*numLabels ] = 0.0f;
		
		labels[ i + 8*numSamples ] = -step_x * i;
		labels[ i + 8*numSamples + numLabels ] = step_y * i;
		labels[ i + 8*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 9*numSamples ] = -step_x * i;
		labels[ i + 9*numSamples + numLabels ] = -step_y * i;
		labels[ i + 9*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 10*numSamples ]  = step_x * i;
		labels[ i + 10*numSamples + numLabels ] = 0.0f;
		labels[ i + 10*numSamples + 2*numLabels ] = step_z * i;

		labels[ i + 11*numSamples ] = step_x * i;
		labels[ i + 11*numSamples + numLabels ] = 0.0f;
		labels[ i + 11*numSamples + 2*numLabels ] = -step_z * i;

		labels[ i + 12*numSamples ] = -step_x * i;
		labels[ i + 12*numSamples + numLabels ] = 0.0f;
		labels[ i + 12*numSamples + 2*numLabels ] = step_z * i;

		labels[ i + 13*numSamples ]  = -step_x * i;
		labels[ i + 13*numSamples + numLabels ] = 0.0f;
		labels[ i + 13*numSamples + 2*numLabels ] = -step_z * i;

		labels[ i + 14*numSamples ] = 0.0f;
		labels[ i + 14*numSamples + numLabels ] = step_y * i;
		labels[ i + 14*numSamples + 2*numLabels ] = step_z * i;

		labels[ i + 15*numSamples ] = 0.0f;
		labels[ i + 15*numSamples + numLabels ] = step_y * i;
		labels[ i + 15*numSamples + 2*numLabels ] = - step_z * i;

		labels[ i + 16*numSamples ] = 0.0f;
		labels[ i + 16*numSamples + numLabels ] = -step_y * i;
		labels[ i + 16*numSamples + 2*numLabels ] = step_z * i;

		labels[ i + 17*numSamples ] = 0.0f;
		labels[ i + 17*numSamples + numLabels ] = -step_y * i;
		labels[ i + 17*numSamples + 2*numLabels ] = -step_z * i;
		}
		else
		{
		// in fast approximation, only 6 directions are considered.
		// (+x, -x, +y, -y, +z, -z)
		labels[ i ] = step_x * i;
		labels[ i + numLabels ] = 0.0f;
		labels[ i + 2*numLabels ] = 0.0f;

		labels[ i + numSamples ] = 0.0f;
		labels[ i + numSamples + numLabels ] = step_y * i;
		labels[ i + numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 2*numSamples ] = 0.0f;
		labels[ i + 2*numSamples + numLabels ] = 0.0f;
		labels[ i + 2*numSamples + 2*numLabels ] = step_z * i;

		labels[ i + 3*numSamples ] = -step_x * i;
		labels[ i + 3*numSamples + numLabels ] = 0.0f;
		labels[ i + 3*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 4*numSamples ] = 0.0f;
		labels[ i + 4*numSamples + numLabels ] = -step_y * i;
		labels[ i + 4*numSamples + 2*numLabels ] = 0.0f;

		labels[ i + 5*numSamples ] = 0.0f;
		labels[ i + 5*numSamples + numLabels ] = 0.0f;
		labels[ i + 5*numSamples + 2*numLabels ] = -step_z * i;
		}
	}
}

void computePairs(int numControlPointsX, int numControlPointsY, int numControlPointsZ, int *pairs)
{
	//1-dimensional array of size 2*numpairs containing the nodes' indices for each MRF edge. The 
	//indices for the two nodes of the i-th node is assumed to be given by: pairs[2*i], pairs[2*i+1]
	
	//indexes for the 3D image
	int x, y, z, xyz, xy;
    
	int numControlPointsXY;
    int node1, node2, index;
    
	numControlPointsXY = numControlPointsX*numControlPointsY;

    index = 0;
    for(z=0; z<numControlPointsZ; z++)
    {
		xyz = z*numControlPointsXY;
        for(x=0; x<numControlPointsX; x++)
        {
			xy = x*numControlPointsY;
            for(y=0; y<numControlPointsY; y++)
            {
                node1 = y + xy + xyz;
                if(y < numControlPointsY - 1)
                {
                    node2 = node1 + 1;
                    pairs[2*index] = node1;
                    pairs[2*index + 1] = node2;
                    index += 1;
                }
                if(x < numControlPointsX - 1)
                {
                    node2 = node1 + numControlPointsY;
                    pairs[2*index] = node1;
                    pairs[2*index + 1] = node2;
                    index += 1;
                }
                if(z < numControlPointsZ - 1)
                {
                    node2 =  node1 + numControlPointsXY;
                    pairs[2*index] = node1;
                    pairs[2*index + 1] = node2;
                    index += 1;
                }
            }
        }
    }
	//printf("index = %d\n", index);
}

void computeWcosts(int numPairs, float regWeight, float *wcosts)
{
	//1-dimensional array of size number of pairs containing the weights w_{pq} used in the MRF pairwise
	//potentials. wcosts[i] is the weight corresponding to teh i-th MRF edge.
	int i;
	for (i=0; i<numPairs; i++)
	{
		wcosts[i] = regWeight; 
	}
}

void computeWcosts_Ou(int numPairs, int numPoints, int numLabels, float regWeight, float *lcosts, float *dist, float *wcosts)
{
	int i,num1;
	float base;
	float sumData, sumReg;
	float meanData, meanReg;
		
	num1 = numPoints * numLabels;
	sumData = 0.0f;
	for (i=0; i<numPoints; i++)
		sumData += lcosts[i*numLabels];
	meanData = sumData/(float)numPoints;
	
	num1 = numLabels * numLabels;
	sumReg = 0.0f;
	for (i=0; i<num1; i++)
	    {
		if ( ~isnan(dist[i]) )
			sumReg += dist[i];
		}
	meanReg = sumReg / (float)num1;
		
	if ( isnan(sumReg) || (sumReg==0.0f) )
		base = 1.0;
	else
		base = meanData/meanReg;
		
	
	if (base<8000.0)
		base=8000.0;
		
	printf("sumData = %f (mean=%f), sumReg = %f (mean=%f), base = %f\n", sumData, meanData, sumReg, meanReg, base);
	for (i=0;i<numPairs;i++)
		wcosts[i] = 3500.0*regWeight;  
		//wcosts[i] = regWeight*base;		
		//wcosts[i] = 8000.0*regWeight;  
	
}


void computeDist(float *labels, int numLabels, Fvector3d voxelSize, int distMethod, float threshold, float *dist, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ)
{
	//the distance function used for defining the MRF pairwise potentials. This is an 1-dimensional array 
	//of size numlabels*numlabels. The distance D(x_p, x_q) is assumed to be given by: dist[x_q*numlabels+x_p].
	//Three different methods are supported. The firs one is the Potts distance, the second one is the truncated
	//quadratic distance while the third one is the truncated absolute difference.

    //help index
    int hp, hq;
    
    //variable for cost
    float cost, dx, dy, dz;
    
    if(distMethod == 0)
    {		
        for(hp=0; hp<numLabels; hp++)
        {
            for(hq=0; hq<numLabels; hq++)
            {
                if(hp == hq)
                {
                    dist[hq*numLabels + hp] = 0;
                }
                else
                {
                    dist[hq*numLabels + hp] = threshold;
                }
            }
        }
    }
    else if(distMethod == 1)
    {
        for(hp=0; hp<numLabels; hp++)
        {
            for(hq=0; hq<numLabels; hq++)
            {
				// -------------------------------
				// Note: on December 8, 2010, Yangming normalized dx by the distBetweenControlPoints, according to Aris' suggestion
				// -------------------------------
				dx = voxelSize.x * (labels[hp] - labels[hq])/(float)distBetweenControlPointsX;
                dy = voxelSize.y * (labels[hp + numLabels] - labels[hq + numLabels])/(float)distBetweenControlPointsY;
                dz = voxelSize.z * (labels[hp + 2*numLabels] - labels[hq + 2*numLabels])/(float)distBetweenControlPointsZ;
				// -------------------------------
				// Note: on August 12, 2010, Yangming added the multiplications with voxel sizes to turn distances in image space into distances in physical space
				// -------------------------------
                //dx = (labels[hp] - labels[hq])*voxelSize.x;
                //dy = (labels[hp + numLabels] - labels[hq + numLabels])*voxelSize.y;
                //dz = (labels[hp + 2*numLabels] - labels[hq + 2*numLabels])*voxelSize.z;
				// -------------------------------
				// old version by Aris, in image space
				// -------------------------------
				//dx = labels[hp] - labels[hq];
                //dy = labels[hp + numLabels] - labels[hq + numLabels];
                //dz = labels[hp + 2*numLabels] - labels[hq + 2*numLabels];
                cost = sqrt(dx*dx + dy*dy + dz*dz);
                if(cost > threshold)
                    cost = threshold;
                dist[hq*numLabels+hp] = cost;
            }
        }
    }
    else if(distMethod == 2)
    {
        for(hp=0; hp<numLabels; hp++)
        {
            for(hq=0; hq<numLabels; hq++)
            {
				// -------------------------------
				// Note: on December 8, 2010, Yangming normalized dx by the distBetweenControlPoints, according to Aris' suggestion
				// -------------------------------
				dx = (labels[hp] - labels[hq])/(float)distBetweenControlPointsX;
                dy = (labels[hp + numLabels] - labels[hq + numLabels])/(float)distBetweenControlPointsY;
                dz = (labels[hp + 2*numLabels] - labels[hq + 2*numLabels])/(float)distBetweenControlPointsZ;
				// -------------------------------
				// Note: on August 12, 2010, Yangming added the multiplications with voxel sizes to turn distances in image space into distances in physical space
				// -------------------------------
                //dx = (labels[hp] - labels[hq])*voxelSize.x;
                //dy = (labels[hp + numLabels] - labels[hq + numLabels])*voxelSize.y;
                //dz = (labels[hp + 2*numLabels] - labels[hq + 2*numLabels])*voxelSize.z;
                cost = fabs(dx*dx + dy*dy + dz*dz);
				//printf("cost = %f, threshold = %f\n", cost, threshold);
                if(cost > threshold)
                    cost = threshold;
                dist[hq*numLabels+hp] = cost;
            }
        }
    }
}

void computeUnary(unsigned char ****SF, unsigned char ****TF, float *labels, int numControlPointsX, int numControlPointsY,
			int numControlPointsZ, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int method, 
			int distBasedWeight, Ivector3d imageSize, Fvector3d resolutionRatio, float ***confidenceMap, Fvector3d ***defField, int numLabels, 
			int numFeatures, int SimilarityMeasure, int chk, unsigned char ***mask, float *lcosts, int AdditionOrComposition, float *weightsX, float *weightsY, float *weightsZ)
{
	//function to compute the 1-dimensioanal array of size numlabels*numpoints containing the label costs (i.e.,
	//the MRF singleton-potentials). The label cost for the k-th label at the i-th node is assumed to be given
	//by lcosts[k*numpoints+i]

	//indexes used to index the graph, the output, the labels and a help one
    int node;
    int indexO;
    int h;
	int itlabels;
    
    //coordinates of the current pixel
    int x, y, z;
    
    //counter used for the normalization of the cost
    float *counter;
    
    //variables used to determine the positions of the nodes that influence
    //the current pixel.
    int u, v, w;
    int dw, dv, du; //(there are going to be used only in the case that a 
                    //cubic b-spline patch is used.
    
    //new position for the current voxel because of the applied deformation
    float nx, ny, nz;
    
    //number of nodes
    int numPoints;

    //in case that a spatial varying weighting is used factor gives the 
    //spatial weight for the current voxel
    float factor;
    //float *weightsX, *weightsY, *weightsZ;
	//int chk = 0;
	
	//cost calculated in every voxel
	float cost;


	/*
	//check if a weighting scheme based on the distance from the nodes will be 
	//used and in that case precompute the weights
	if(distBasedWeight == 1)
	{
		//nearest neighbor weighting scheme
		if(method == 0)
		{
			weightsX = NULL; weightsY = NULL; weightsZ = NULL; chk = 1;
		}
		//linear weighting scheme
		else if(method == 1)   
		{
			weightsX = (float *)malloc(2*imageSize.x*sizeof(float));
			weightsY = (float *)malloc(2*imageSize.y*sizeof(float));
			weightsZ = (float *)malloc(2*imageSize.z*sizeof(float));			
			precomputeWeights(imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, resolutionRatio, method, weightsX, weightsY, weightsZ);
			chk = 1;
		}
		// cubic B spline weighting scheme
		else if(method == 2)
		{
			weightsX = (float *)malloc(4*imageSize.x*sizeof(float));
			weightsY = (float *)malloc(4*imageSize.y*sizeof(float));
			weightsZ = (float *)malloc(4*imageSize.z*sizeof(float));
			precomputeWeights(imageSize, distBetweenControlPointsX, distBetweenControlPointsY, distBetweenControlPointsZ, resolutionRatio, method,	weightsX, weightsY, weightsZ);
			chk = 1;
		}
	}
	*/
	
	numPoints = numControlPointsX*numControlPointsY*numControlPointsZ;
	int numPointsLabels=numPoints*numLabels;
    counter = (float *)calloc(numPoints*numLabels,sizeof(float));    
	for (h=0;h<numPointsLabels;h++) {
	  counter[h]=0.0f;
	  lcosts[h]=0.0f;
	  }
	  
	  
    //calculate costs: the idea is that you traverse the target image, you
    //calculate the cost due to the deformation implied by the label that
    //is checked and then attribute it to the node that influence this
    //voxel depending on how the patch is defined.
    for(z = 0; z < imageSize.z; z++)
    {    
        //printf("  ... %d/%d ...  ", z+1, imageSize.z);	
        for(x = 0; x < imageSize.x; x++)
        {            
            for(y = 0; y < imageSize.y; y++)
            {    
//printf("[%d] node\n", omp_get_thread_num());            
                //check if the voxel belongs to foreground or not
                if(mask[z][x][y] > 0)
                {
                    //for all possible translations due to the application
                    //of the labels
//if ( (x==62)&&(y==28)&&(z==30)  )  printf("(x,y,z)=(%d,%d,%d)\n",x,y,z);
                    for(itlabels = 0; itlabels < numLabels; itlabels++)
                    {
                        //compute the new position of the current voxel due
                        //to the translation induced by the application of
                        //the labels.
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha0\n", itlabels+1, numLabels);			
						nx = (float)x + defField[z][x][y].x + labels[itlabels];
						ny = (float)y + defField[z][x][y].y + labels[itlabels + numLabels];
						nz = (float)z + defField[z][x][y].z + labels[itlabels + 2*numLabels];
						
                        cost = featureCost(TF, SF, confidenceMap, x, y, z, nx, ny, nz, numFeatures, SimilarityMeasure, imageSize);   // this cost already inccludes the role of confidenceMap
                        //attributing the cost to the appropriate graph node
                        u = (int)floor((float)x/(float)distBetweenControlPointsX);
                        v = (int)floor((float)y/(float)distBetweenControlPointsY);
                        w = (int)floor((float)(z-1)/(float)distBetweenControlPointsZ);
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha1\n", itlabels+1, numLabels);
                        //the attribution of the cost will vary depending 
                        //on how the patch is defined.
                        
                        //in the case that a nearest neighbor scheme has 
                        //been chosen
                        if (method==0)
                        {
                            u = (int)myround((float)x/(float)distBetweenControlPointsX);
                            v = (int)myround((float)y/(float)distBetweenControlPointsY);
                            w = (int)myround((float)(z-1)/(float)distBetweenControlPointsZ);
                            
                            node = u+v*numControlPointsX+w*numControlPointsX*numControlPointsY;
							indexO = node +itlabels*numPoints;
							if (indexO>=0){
								lcosts[indexO] += cost;
								counter[indexO] += 1;
								}
                        }      
                        //in the case that linear weighting scheme has been chosen
                        else if (method==1)
                        {
                            if(chk == 1)                          
                            {
                                factor = weightsX[x]*weightsY[y]*weightsZ[z];
                            }
							else
                            {
                                factor = 1;
                            }                                                                                    
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha2\n", itlabels+1, numLabels);
                            node = u+ v*numControlPointsX+ w*numControlPointsX*numControlPointsY;
                            indexO = node +itlabels*numPoints;
							if (indexO>=0){
								lcosts[indexO] += cost * factor;
								counter[indexO] += factor;
								}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d,haha3\n", itlabels+1, numLabels);
                            if(u+1 < numControlPointsX)
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha4\n", itlabels+1, numLabels);
							  if(chk == 1)
                                {
                                    factor = weightsX[x+imageSize.x]*weightsY[y]*weightsZ[z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = (u+1)+ v*numControlPointsX+ w*numControlPointsX*numControlPointsY;
                                indexO = node +itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor;
									counter[indexO] += factor;
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");								
                            }
                            
                            if(v+1 < numControlPointsY)
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha5...", itlabels+1, numLabels);
                                if(chk == 1)
                                {
                                    factor = weightsX[x]*weightsY[y+imageSize.y]*weightsZ[z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = u+ (v+1)*numControlPointsX+ w*numControlPointsX*numControlPointsY;
                                indexO = node +itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor;
									counter[indexO] += factor;
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");		
                            }
                            if(w+1 < numControlPointsZ)
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha6...", itlabels+1, numLabels);
                                if(chk == 1)
                                {
                                    factor = weightsX[x]*weightsY[y]*weightsZ[z+imageSize.z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = u+ v*numControlPointsX+ (w+1)*numControlPointsX*numControlPointsY;
                                indexO = node +itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor;
									counter[indexO] += factor;
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");										
                            }
                            if((u+1 < numControlPointsX)&&(v+1 < numControlPointsY))
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha7...", itlabels+1, numLabels);
                                if(chk == 1)
                                {
                                    factor = weightsX[x+imageSize.x]*weightsY[y+imageSize.y]*weightsZ[z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = (u+1)+ (v+1)*numControlPointsX+ w*numControlPointsX*numControlPointsY;
								indexO = node +	itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor * confidenceMap[z][x][y];
									counter[indexO] += factor * confidenceMap[z][x][y];
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");										
                            }
                            if((w+1 < numControlPointsZ)&&(v+1 < numControlPointsY))
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha8...", itlabels+1, numLabels);
                                if(chk == 1)
                                {
                                    factor = weightsX[x]*weightsY[y+imageSize.y]*weightsZ[z+imageSize.z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = u+ (v+1)*numControlPointsX+ (w+1)*numControlPointsX*numControlPointsY;
                                indexO = node +itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor * confidenceMap[z][x][y];
									counter[indexO] += factor * confidenceMap[z][x][y];
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");
                            }
                            if((u+1 < numControlPointsX)&&(w+1 < numControlPointsZ))
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha9...", itlabels+1, numLabels);
                                if(chk == 1)
                                {
                                    factor = weightsX[x+imageSize.x]*weightsY[y]*weightsZ[z+imageSize.z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = (u+1)+ v*numControlPointsX+ (w+1)*numControlPointsX*numControlPointsY;
                                indexO = node +itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor;
									counter[indexO] += factor;
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");										
                            }
                            if((u+1 < numControlPointsX)&&(v+1 < numControlPointsY)&&(w+1 < numControlPointsZ))
                            {
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("%d/%d, haha10...", itlabels+1, numLabels);
                                if(chk == 1)
                                {
                                    factor = weightsX[x+imageSize.x]*weightsY[y+imageSize.y]*weightsZ[z+imageSize.z];
                                }
                                else
                                {
                                    factor = 1;
                                }
                                node = (u+1)+ (v+1)*numControlPointsX+ (w+1)*numControlPointsX*numControlPointsY;
                                indexO = node +itlabels*numPoints;
								if (indexO>=0){
									lcosts[indexO] += cost * factor;
									counter[indexO] += factor;
									}
//if ( (x==62)&&(y==28)&&(z==30)  ) printf("done...\n");
                            }
                        }
                        //cubic B-splines patch
                        else if (method==2)
                        {                            
                            u -= 1; v -= 1; w -= 1;
                            
                            for (dw=0; dw<4; dw++)
                            {
                                for (dv=0; dv<4; dv++)
                                {
                                    for (du=0; du<4; du++)
                                    {
                                        if((u+du>=0)&&(u+du<numControlPointsX)
                                        && (v+dv>=0)&&(v+dv<numControlPointsY)
                                        && (w+dw>=0)&&(w+dw<numControlPointsZ))
                                        {
                                            if(chk == 1)
                                            {
                                                factor = weightsX[du*imageSize.x+x]*weightsY[dv*imageSize.y+y]*weightsZ[dw*imageSize.z+z];
                                            }
                                            else
                                            {
                                                factor = 1;
                                            }
                                            node = u+du+(v+dv)*numControlPointsX+(w+dw)*numControlPointsX*numControlPointsY;
                                            indexO = node +itlabels*numPoints;

											if (indexO>=0){
												lcosts[indexO] += cost * factor;
												counter[indexO] += factor;
												}
                                        }
                                    }// for
                                }//for 
                            }//for
                        }//else                        
                    }//for
                }//if
            } // for y
        }//for x
    }//for z
    
	
    for(h=0; h<numPointsLabels; h++)
    {
        if(counter[h]){
		    lcosts[h] /= counter[h];
		}
    }
	
//printf("test1\n");
    free(counter);
//printf("test2\n");	

	/*
	if(chk == 1)
	{
		free(weightsX);
		free(weightsY);
		free(weightsZ);
	}
	*/
//printf("test3\n");	
}

void labels2deformations(int *optimalLabels, float *labels, int numControlPointsX, int numControlPointsY, 
						 int numControlPointsZ, Fvector3d ***controlPointDisp, int numLabels)
{
	//function that writes the output of the optimization function in an appropriate form.
	//nodes indexes
	int u, v, w;
	//help index
	int h;
	int numControlPoints=numControlPointsX*numControlPointsY*numControlPointsZ;

	for(h=0; h<numControlPoints; h++)
	{
		w = (int)floor((float)h/((float)numControlPointsX*numControlPointsY));
		v = (int)floor(((float)(h - w*numControlPointsX*numControlPointsY))/((float)numControlPointsX));
		u = h - v*numControlPointsX - w*numControlPointsX*numControlPointsY;
		controlPointDisp[w][u][v].x = labels[optimalLabels[h]];
		controlPointDisp[w][u][v].y = labels[optimalLabels[h]+ numLabels];
 		controlPointDisp[w][u][v].z = labels[optimalLabels[h]+ numLabels*2];
	}
}


