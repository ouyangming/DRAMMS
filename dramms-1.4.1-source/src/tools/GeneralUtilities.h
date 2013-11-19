#ifndef __GENERALUTIL_H__
#define __GENERALUTIL_H__

#include "common/mvcd.h"

int myround(float x);
float featureCost(unsigned char ****trg, unsigned char ****src, float ***confidenceMap, int x_trg, int y_trg,
				  int z_trg, float x_src, float y_src, float z_src, int numFeatures, int SimilarityMeasure, Ivector3d imageSize);
//float myTrilinearInterpolation(unsigned char ***Img, float x, float y, float z, Ivector3d imageSize);
void precomputeWeights(Ivector3d imageSize, int distBetweenControlPointsX, int distBetweenControlPointsY, int distBetweenControlPointsZ, Fvector3d resolutionRatio, int method,
					   float *weightsX, float *weightsY, float *weightsZ);

#endif


