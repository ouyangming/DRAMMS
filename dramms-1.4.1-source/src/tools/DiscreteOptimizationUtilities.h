/**
 * @file  DiscreteOptimizationUtilities.h
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



#ifndef __DISCRETEOPTUT_H__
#define __DISCRETEOPTUT_H__

#include "common/mvcd.h"

void computeLabels(int numSamples, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d resolutionRatio, int levelIndex, int method, float *labels, float label_factor, int iter,  int numIterInEachResolution, int indexGridLevel, int numGridLevels, int fastApproximationOrNot);
void computePairs(int numControlPointsX, int numControlPointsY, int numControlPointsZ, int *pairs);
void computeWcosts(int numPairs, float regWeight, float *wcosts);
void computeWcosts_Ou(int numPairs, int numPoints, int numLabels, float regWeight, float *lcosts, float *dist, float *wcosts);
void computeDist(float *labels, int numLabels, Fvector3d voxelSize, int distMethod, float threshold, float *dist, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ);
void computeUnary(unsigned char ****SF, unsigned char ****TF, float *labels, int numControlPointsX, int numControlPointsY,
					 int numControlPointsZ, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, int method, 
					 int distBasedWeight, Ivector3d imageSize, Fvector3d resolutionRatio, float ***confidenceMap, Fvector3d ***defField, int numLabels,
					 int numFeatures, int SimilarityMeasure, int chk, unsigned char ***mask, float *lcosts, int AdditionOrComposition, float *weightsX, float *weightsY, float *weightZ);
void labels2deformations(int *optimalLabels, float *labels, int numControlPointsX, int numControlPointsY,
						 int numControlPointsZ, Fvector3d ***controlPointDisp, int numLabels);
#endif


