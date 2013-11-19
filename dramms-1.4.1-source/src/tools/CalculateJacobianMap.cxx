/**
 * @file  CalculateJacobianMap.cxx
 * @brief Calculate jacobian determinante map.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <string>
#include <iostream> // cout, cerr, endl
#include <math.h>   // M_PI, exp(), pow(), sqrt(),...
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <common/mvcd.h>
#include <common/imageio.h>
#include <common/cres.h>
#include <common/matrix.h>
#include <common/general.h>


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
    cout << "This program calculates jacobian determinants of the input deformation field. " << endl;
	cout << "-------------------------------------------------" << endl << endl;
    cout << "Usage:" << endl;
    cout << "  " << exec_name << " <input deformation> <output jacobian map>" << endl;
    cout << endl;
    cout << "Required arguments:" << endl;
	cout << "  <input deformation>    file of deformation field." << endl;
    cout << "  <output jacobian map>  file path of Jacobian image." << endl;
    cout << endl;
    cout << "Optional arguments:" << endl;
	cout << "  -W   : x-y switch off (default: on for default input deformation in DRAMMS format)" << endl;
	cout << "  -s   : smooth deformation field before calculating Jacobian map (default: off)" <<endl;
	cout << "  -S   : smooth calculated Jacobian map (default: off)" <<endl;
	cout << "  -C   : cut-off Jacobian values at 0 (set all negative Jacobians to be 0) (default: off)" << endl;
	cout << "  -L   : take log(JacobianDet) Note: all non-positive Jacobian determinant values will be mapped to -100.0 in log(JacDet) calculation. (default: off)" << endl;
	cout << endl;
	cout << "Example: " << endl;
	cout << "  " << exec_name <<" DF.def JacobianMap.img -s -S" << endl;
	cout << endl;
    print_contact();
    cout << endl;
	cout << " Note: the output Jacobian image is in float data type." << endl << endl;
}


// ---------------------------------------------------------------------------
double GaussianFunction(double input,double sigma)
{
  if (input<0)
    return 0;
  else
    return exp(-input/(2.*sigma*sigma));
}


// ---------------------------------------------------------------------------
void generateGaussianSmoothTemplate3D(int smoothSizeXY, int smoothSizeZ, float smoothKernelXY, float smoothKernelZ, Matrix** gaussianSmoothTemplate)
{
  int i,j,k;
  double x,y,z;
  double sumTemplate=0.0;
 
  for (k=0;k<smoothSizeZ;k++)
   for (i=0;i<smoothSizeXY;i++)
     for (j=0;j<smoothSizeXY;j++)
	   {
		z = (double)k-floor((float)smoothSizeZ/2.0);
		x = (double)i-floor((float)smoothSizeXY/2.0);
		y = (double)j-floor((float)smoothSizeXY/2.0);
		
		gaussianSmoothTemplate[k]->data[i][j] = GaussianFunction(pow(z,2.0), smoothKernelZ) * GaussianFunction(pow(x,2.0), smoothKernelXY) * GaussianFunction(pow(y,2.0), smoothKernelXY);
		
		sumTemplate += gaussianSmoothTemplate[k]->data[i][j];
	   }
  
  for (k=0;k<smoothSizeZ;k++)
   {
   for (i=0;i<smoothSizeXY;i++)
     {
     for (j=0;j<smoothSizeXY;j++)
	   {
		gaussianSmoothTemplate[k]->data[i][j] /= sumTemplate;
	   }
	 }
   }
}




// ---------------------------------------------------------
void smooth3DFloat(float ***FloatImg, Matrix **gaussianSmoothTemplate, Ivector3d imageSize, int smoothSizeXY, int smoothSizeZ)
{
  int i,j,k;
  int x,y,z;
  Ivector3d smoothCenter;
  Ivector3d gridPoint;
  float min=1000.0, max=0.0;
  bool breakOrNot;
  
  
  float ***FloatImg_backup;  
  FloatImg_backup = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++)
	    {
	    FloatImg_backup[k][i][j] = FloatImg[k][i][j];
		if (FloatImg_backup[k][i][j]>max)  max=FloatImg_backup[k][i][j];
		if (FloatImg_backup[k][i][j]<min)  min=FloatImg_backup[k][i][j];  
		}
  printf("\n(min,max) = (%f,%f) before smoothing;\n",min,max);
  
  
  // smoothing
  min=1000.0;
  max=0.0;
  printf("Start smoothing Jacobian Determinant Image...\n");
  for (k=0;k<imageSize.z;k++)
    {
	printf("z=%d...\t", k);
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++)
	    {
		  FloatImg[k][i][j] = 0.0;
		  breakOrNot=false;
		  
		  smoothCenter.x = i;
		  smoothCenter.y = j;
		  smoothCenter.z = k;
		  
		  for (z=0;z<smoothSizeZ;z++)
		    for (x=0;x<smoothSizeXY;x++)
			  for (y=0;y<smoothSizeXY;y++)
			    {
				  gridPoint.x = smoothCenter.x-(int)floor(smoothSizeXY/2)+x;
				  gridPoint.y = smoothCenter.y-(int)floor(smoothSizeXY/2)+y;
				  gridPoint.z = smoothCenter.z-(int)floor(smoothSizeZ/2)+z;
				  
				  if ( (gridPoint.x>=0)&&(gridPoint.x<=(imageSize.x-1))&&(gridPoint.y>=0)&&(gridPoint.y<=(imageSize.y-1))&&(gridPoint.z>=0)&&(gridPoint.z<=(imageSize.z-1)) )
	  			    FloatImg[k][i][j] += FloatImg_backup[gridPoint.z][gridPoint.x][gridPoint.y]*gaussianSmoothTemplate[z]->data[x][y];
				  else
				    breakOrNot=true;				
				}
		  if (breakOrNot)   FloatImg[k][i][j] = FloatImg_backup[k][i][j];
		  
		  
		  if (FloatImg[k][i][j]>max)   max=FloatImg[k][i][j];
		  if (FloatImg[k][i][j]<min)   min=FloatImg[k][i][j];
		}
	} // exhaustic visit to every point in 3D image
  printf("\n(min,max) = (%f,%f) after smoothing;\n", min,max);
		
  Ffree3d(FloatImg_backup, imageSize.z, imageSize.x);  
}


// ------------------------------------------------------------------
void smoothDeformationField(Fvector3d ***deformationField, Matrix **gaussianSmoothTemplate, Ivector3d imageSize, int smoothSizeXY, int smoothSizeZ)
{
  int i,j,k;
  int x,y,z;
  bool breakLabel;
  Ivector3d smoothCenter;
  Ivector3d gridPoint;
  
  // copy dfImg{X,Y,Z} to dfImg{X,Y,Z}_backup
  float ***dfImgX_backup, ***dfImgY_backup, ***dfImgZ_backup;
  dfImgX_backup = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
  dfImgY_backup = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
  dfImgZ_backup = Falloc3d(imageSize.x, imageSize.y, imageSize.z);
  for (k=0;k<imageSize.z;k++)
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++)
	    {
		  dfImgX_backup[k][i][j] = deformationField[k][i][j].x;
		  dfImgY_backup[k][i][j] = deformationField[k][i][j].y;
		  dfImgZ_backup[k][i][j] = deformationField[k][i][j].z;
		}
		
  // smoothing
  printf("Start smoothing deformation field...\n");
  for (k=0;k<imageSize.z;k++)
    {
    for (i=0;i<imageSize.x;i++)
	  for (j=0;j<imageSize.y;j++)
	    {
		  breakLabel = false;
		  
		  smoothCenter.x = i;
		  smoothCenter.y = j;
		  smoothCenter.z = k;
		  
		  deformationField[k][i][j].x=0;
		  deformationField[k][i][j].y=0;
		  deformationField[k][i][j].z=0;
		  
		  for (z=0;z<smoothSizeZ;z++)
		    for (x=0;x<smoothSizeXY;x++)
			  for (y=0;y<smoothSizeXY;y++)
			    {
				  gridPoint.x = smoothCenter.x-(int)floor(smoothSizeXY/2)+x;
				  gridPoint.y = smoothCenter.y-(int)floor(smoothSizeXY/2)+y;
				  gridPoint.z = smoothCenter.z-(int)floor(smoothSizeZ/2)+z;
				  
				  if ( (gridPoint.x>0)&(gridPoint.x<(imageSize.x-1))&(gridPoint.y>0)&(gridPoint.y<(imageSize.y-1))&(gridPoint.z>0)&(gridPoint.z<(imageSize.z-1)) )
				    {
					  deformationField[k][i][j].x += dfImgX_backup[gridPoint.z][gridPoint.x][gridPoint.y]*gaussianSmoothTemplate[z]->data[x][y];
					  deformationField[k][i][j].y += dfImgY_backup[gridPoint.z][gridPoint.x][gridPoint.y]*gaussianSmoothTemplate[z]->data[x][y];
					  deformationField[k][i][j].z += dfImgZ_backup[gridPoint.z][gridPoint.x][gridPoint.y]*gaussianSmoothTemplate[z]->data[x][y];
					}
				  else
				    breakLabel=true;				  
				}
		  
		  if (breakLabel)
		    {
		    deformationField[k][i][j].x = dfImgX_backup[k][i][j];
			deformationField[k][i][j].y = dfImgY_backup[k][i][j];
			deformationField[k][i][j].z = dfImgZ_backup[k][i][j];
			}
		}
	} // exhaustic visit to every point in 3D image
		
  Ffree3d(dfImgX_backup, imageSize.z, imageSize.x);
  Ffree3d(dfImgY_backup, imageSize.z, imageSize.x);
  Ffree3d(dfImgZ_backup, imageSize.z, imageSize.x);
  
}



// ---------------------------------------------------------------
void Determinant(Matrix *Jaco, float *result) 
{
  int i;
  Matrix *EigenVector ;
  float  *EigenValue;

  CreateMatrix(&EigenVector, Jaco->height, Jaco->height);
  EigenValue = vectorSHEN(0, Jaco->height-1) ;

  Mat_Calculate_EigenVectors_EigenValues(Jaco, EigenValue, EigenVector, FALSE) ;

  (*result) = 1.0 ;
  for(i=0; i<Jaco->height; i++)
    (*result) *= EigenValue[i] ;

  /* free */
  FreeMatrix(EigenVector) ;
  free_vectorSHEN(EigenValue, 0, Jaco->height-1) ;
}


// -------------------------------------------------
void Jacobians(  Fvector3d ***DeformFld, int DRAMMSDeformationOrNot, int x_size, int y_size, int z_size, Fvector3d voxelSize, float ***JacoDeterminant )
{
  int i, j, k;
  int ip, im, jp, jm, kp, km;
  float         max, min, current ;
  int singularityExistedOrNot=NNO;
  int counterSingularity=0;

  
  Matrix *Jaco ;
  CreateMatrix(&Jaco, 3, 3);

  max = 0 ;
  min = 10000.0 ;

  
  for(k=0; k<z_size; k++)
    {
      for(i=0; i<x_size; i++)
	    for(j=0; j<y_size; j++)
	    {
		JacoDeterminant[k][i][j]=0.0;


		//---------------
		// Ou rewrite this jacobian calculation on Feb. 20, 2010, to calculate this for the last column last row and last slice in the image space  (Begin)
		//---------------
		if (DRAMMSDeformationOrNot==NNO)
		{ // if deformation is not DRAMMS format
			if ( i==0 )
				{
				ip=MIN(i+1, x_size-1);
				Jaco->data[1][0] = 0 + (DeformFld[k][ip][j].x - DeformFld[k][i][j].x)*voxelSize.y/voxelSize.x ;   
				Jaco->data[0][0] = 1 + (DeformFld[k][ip][j].y - DeformFld[k][i][j].y)*voxelSize.x/voxelSize.x ;
				Jaco->data[2][0] = 0 + (DeformFld[k][ip][j].z - DeformFld[k][i][j].z)*voxelSize.z/voxelSize.x ;	
				}
			else if ( i==(x_size-1) )
				{
				im=MAX(i-1, 0);
				Jaco->data[1][0] = 0 + (DeformFld[k][i][j].x - DeformFld[k][im][j].x)*voxelSize.y/voxelSize.x ;   
				Jaco->data[0][0] = 1 + (DeformFld[k][i][j].y - DeformFld[k][im][j].y)*voxelSize.x/voxelSize.x ;
				Jaco->data[2][0] = 0 + (DeformFld[k][i][j].z - DeformFld[k][im][j].z)*voxelSize.z/voxelSize.x ;	
				}
			else 
				{
				ip=MIN(i+1, x_size-1);
				im=MAX(i-1, 0);
				Jaco->data[1][0] = 0 + (DeformFld[k][ip][j].x - DeformFld[k][im][j].x)*voxelSize.y/(2.0*voxelSize.x) ;   
				Jaco->data[0][0] = 1 + (DeformFld[k][ip][j].y - DeformFld[k][im][j].y)*voxelSize.x/(2.0*voxelSize.x) ;
				Jaco->data[2][0] = 0 + (DeformFld[k][ip][j].z - DeformFld[k][im][j].z)*voxelSize.z/(2.0*voxelSize.x) ;
				}	
		   
		   	   
			if ( j==0 )
				{
				jp=MIN(j+1, y_size-1);
				Jaco->data[1][1] = 1 + (DeformFld[k][i][jp].x - DeformFld[k][i][j].x)*voxelSize.y/voxelSize.y ;
				Jaco->data[0][1] = 0 + (DeformFld[k][i][jp].y - DeformFld[k][i][j].y)*voxelSize.x/voxelSize.y ;
				Jaco->data[2][1] = 0 + (DeformFld[k][i][jp].z - DeformFld[k][i][j].z)*voxelSize.z/voxelSize.y ;
				}
			else if ( j==(y_size-1) )
				{
				jm=MAX(j-1, 0);
				Jaco->data[1][1] = 1 + (DeformFld[k][i][j].x - DeformFld[k][i][jm].x)*voxelSize.y/voxelSize.y ;
				Jaco->data[0][1] = 0 + (DeformFld[k][i][j].y - DeformFld[k][i][jm].y)*voxelSize.x/voxelSize.y ;
				Jaco->data[2][1] = 0 + (DeformFld[k][i][j].z - DeformFld[k][i][jm].z)*voxelSize.z/voxelSize.y ;
				}
			else
				{
				jp=MIN(j+1, y_size-1);
				jm=MAX(j-1, 0);
				Jaco->data[1][1] = 1 + (DeformFld[k][i][jp].x - DeformFld[k][i][jm].x)*voxelSize.y/(2.0*voxelSize.y);
				Jaco->data[0][1] = 0 + (DeformFld[k][i][jp].y - DeformFld[k][i][jm].y)*voxelSize.x/(2.0*voxelSize.y);
				Jaco->data[2][1] = 0 + (DeformFld[k][i][jp].z - DeformFld[k][i][jm].z)*voxelSize.z/(2.0*voxelSize.y);
				}
		
		
		   
			if ( k==0 )
				{
				kp=MIN(k+1, z_size-1);
				Jaco->data[1][2] = 0 + (DeformFld[kp][i][j].x - DeformFld[k][i][j].x)*voxelSize.y/voxelSize.z;
				Jaco->data[0][2] = 0 + (DeformFld[kp][i][j].y - DeformFld[k][i][j].y)*voxelSize.x/voxelSize.z;
				Jaco->data[2][2] = 1 + (DeformFld[kp][i][j].z - DeformFld[k][i][j].z)*voxelSize.z/voxelSize.z;
				}
			else if ( k==(z_size-1) )
				{
				km=MAX(k-1, 0);
				Jaco->data[1][2] = 0 + (DeformFld[k][i][j].x - DeformFld[km][i][j].x)*voxelSize.y/voxelSize.z;
				Jaco->data[0][2] = 0 + (DeformFld[k][i][j].y - DeformFld[km][i][j].y)*voxelSize.x/voxelSize.z;
				Jaco->data[2][2] = 1 + (DeformFld[k][i][j].z - DeformFld[km][i][j].z)*voxelSize.z/voxelSize.z;
				}
			else
				{
				kp=MIN(k+1, z_size-1);
				km=MAX(k-1, 0);
				Jaco->data[1][2] = 0 + (DeformFld[kp][i][j].x - DeformFld[km][i][j].x)*voxelSize.y/(2.0*voxelSize.z);
				Jaco->data[0][2] = 0 + (DeformFld[kp][i][j].y - DeformFld[km][i][j].y)*voxelSize.x/(2.0*voxelSize.z);
				Jaco->data[2][2] = 1 + (DeformFld[kp][i][j].z - DeformFld[km][i][j].z)*voxelSize.z/(2.0*voxelSize.z);
				}
		} // if deformation is not DRAMMS format
		else
		{ // if deformation is DRAMMS format
			if ( i==0 )
				{
				ip=MIN(i+1, x_size-1);
				Jaco->data[1][1] = 1 + (DeformFld[k][ip][j].x - DeformFld[k][i][j].x)*voxelSize.x/voxelSize.x ;   
				Jaco->data[0][1] = 0 + (DeformFld[k][ip][j].y - DeformFld[k][i][j].y)*voxelSize.y/voxelSize.x ;
				Jaco->data[2][1] = 0 + (DeformFld[k][ip][j].z - DeformFld[k][i][j].z)*voxelSize.z/voxelSize.x ;	
				}
			else if ( i==(x_size-1) )
				{
				im=MAX(i-1, 0);
				Jaco->data[1][1] = 1 + (DeformFld[k][i][j].x - DeformFld[k][im][j].x)*voxelSize.x/voxelSize.x ;   
				Jaco->data[0][1] = 0 + (DeformFld[k][i][j].y - DeformFld[k][im][j].y)*voxelSize.y/voxelSize.x ;
				Jaco->data[2][1] = 0 + (DeformFld[k][i][j].z - DeformFld[k][im][j].z)*voxelSize.z/voxelSize.x ;	
				}
			else 
				{
				ip=MIN(i+1, x_size-1);
				im=MAX(i-1, 0);
				Jaco->data[1][1] = 1 + (DeformFld[k][ip][j].x - DeformFld[k][im][j].x)*voxelSize.x/(2.0*voxelSize.x) ;   
				Jaco->data[0][1] = 0 + (DeformFld[k][ip][j].y - DeformFld[k][im][j].y)*voxelSize.y/(2.0*voxelSize.x) ;
				Jaco->data[2][1] = 0 + (DeformFld[k][ip][j].z - DeformFld[k][im][j].z)*voxelSize.z/(2.0*voxelSize.x) ;
				}
		   
		   
		   
		   
			if ( j==0 )
				{
				jp=MIN(j+1, y_size-1);
				Jaco->data[1][0] = 0 + (DeformFld[k][i][jp].x - DeformFld[k][i][j].x)*voxelSize.x/voxelSize.y ;
				Jaco->data[0][0] = 1 + (DeformFld[k][i][jp].y - DeformFld[k][i][j].y)*voxelSize.y/voxelSize.y ;
				Jaco->data[2][0] = 0 + (DeformFld[k][i][jp].z - DeformFld[k][i][j].z)*voxelSize.z/voxelSize.y ;
				}
			else if ( j==(y_size-1) )
				{
				jm=MAX(j-1, 0);
				Jaco->data[1][0] = 0 + (DeformFld[k][i][j].x - DeformFld[k][i][jm].x)*voxelSize.x/voxelSize.y ;
				Jaco->data[0][0] = 1 + (DeformFld[k][i][j].y - DeformFld[k][i][jm].y)*voxelSize.y/voxelSize.y ;
				Jaco->data[2][0] = 0 + (DeformFld[k][i][j].z - DeformFld[k][i][jm].z)*voxelSize.z/voxelSize.y ;
				}
			else
				{
				jp=MIN(j+1, y_size-1);
				jm=MAX(j-1, 0);
				Jaco->data[1][0] = 0 + (DeformFld[k][i][jp].x - DeformFld[k][i][jm].x)*voxelSize.x/(2.0*voxelSize.y);
				Jaco->data[0][0] = 1 + (DeformFld[k][i][jp].y - DeformFld[k][i][jm].y)*voxelSize.y/(2.0*voxelSize.y);
				Jaco->data[2][0] = 0 + (DeformFld[k][i][jp].z - DeformFld[k][i][jm].z)*voxelSize.z/(2.0*voxelSize.y);
				}
		
		
		   
			if ( k==0 )
				{
				kp=MIN(k+1, z_size-1);
				Jaco->data[1][2] = 0 + (DeformFld[kp][i][j].x - DeformFld[k][i][j].x)*voxelSize.x/voxelSize.z;
				Jaco->data[0][2] = 0 + (DeformFld[kp][i][j].y - DeformFld[k][i][j].y)*voxelSize.y/voxelSize.z;
				Jaco->data[2][2] = 1 + (DeformFld[kp][i][j].z - DeformFld[k][i][j].z)*voxelSize.z/voxelSize.z;
				}
			else if ( k==(z_size-1) )
				{
				km=MAX(k-1, 0);
				Jaco->data[1][2] = 0 + (DeformFld[k][i][j].x - DeformFld[km][i][j].x)*voxelSize.x/voxelSize.z;
				Jaco->data[0][2] = 0 + (DeformFld[k][i][j].y - DeformFld[km][i][j].y)*voxelSize.y/voxelSize.z;
				Jaco->data[2][2] = 1 + (DeformFld[k][i][j].z - DeformFld[km][i][j].z)*voxelSize.z/voxelSize.z;
				}
			else
				{
				kp=MIN(k+1, z_size-1);
				km=MAX(k-1, 0);
				Jaco->data[1][2] = 0 + (DeformFld[kp][i][j].x - DeformFld[km][i][j].x)*voxelSize.x/(2.0*voxelSize.z);
				Jaco->data[0][2] = 0 + (DeformFld[kp][i][j].y - DeformFld[km][i][j].y)*voxelSize.y/(2.0*voxelSize.z);
				Jaco->data[2][2] = 1 + (DeformFld[kp][i][j].z - DeformFld[km][i][j].z)*voxelSize.z/(2.0*voxelSize.z);
				}
		} // if deformation is DRAMMS format
		
		
		
		//Determinant(Jaco, &current) ;
	    current = Jaco->data[0][0]*(Jaco->data[1][1]*Jaco->data[2][2]-Jaco->data[2][1]*Jaco->data[1][2]) - Jaco->data[1][0]*(Jaco->data[0][1]*Jaco->data[2][2]-Jaco->data[2][1]*Jaco->data[0][2]) + Jaco->data[2][0]*(Jaco->data[0][1]*Jaco->data[1][2]-Jaco->data[1][1]*Jaco->data[0][2]) ;
		JacoDeterminant[k][i][j] = current ;	  
		
			
		//---------------
		// Ou rewrite this jacobian calculation on Feb. 20, 2010, to calculate this for the last column las row and las slice in the image space  (Begin)
		//---------------
		if( JacoDeterminant[k][i][j]>max )
	        max = JacoDeterminant[k][i][j] ;
	    if( JacoDeterminant[k][i][j]<min )
	        min = JacoDeterminant[k][i][j] ;
			
		if ( JacoDeterminant[k][i][j]<0.0 )
			{
			singularityExistedOrNot = YYES;
			counterSingularity++;
			}
	  }
    }

	printf("\n");
	printf("SingularityExistedOrNot=%d\n", singularityExistedOrNot);
	printf("PercentOfSingularVoxels=%f\n", (float)counterSingularity/(float)(x_size*y_size*z_size));
	printf("MinJacobian=%f\n", min);
	printf("MaxJacobian=%f\n", max);


  FreeMatrix(Jaco) ;
}




int main(int argc,char *argv[])
{
    if (argc == 1) {
        print_help();
        exit(1);
    }

  Ivector3d     imageSize;
  Fvector3d     voxelSize;
  int           i,j,k;
  int 			JacobianSmoothOrNot=NNO;
  int			DeformationSmoothOrNot=NNO;
  int 			DRAMMSDeformationOrNot=YYES;
  int 			CutOffJacobianDetAtZeroOrNot=NNO;
  int 			LogJacobianDetOrNot=NNO;
  bool			ok = true;
  
  Image* deformation = NULL;
	
  
  int c=-1;
  while((c=getopt(argc,argv,"WsSCLhv")) != -1)
    {
      switch(c)
	{
	case 'W':
	    DRAMMSDeformationOrNot=NNO;
		break;
	
	case 's':
		DeformationSmoothOrNot=YYES;
		break;
		
	case 'S':
		JacobianSmoothOrNot=YYES;
		break;
		
	case 'C':
		CutOffJacobianDetAtZeroOrNot=YYES;
		break;
		
	case 'L':
		LogJacobianDetOrNot=YYES;
		break;

    case 'v':
        // ignore
        break;

    case 'h':
        print_help();
        exit(0);

	default:
        // error message printed by getopt() already
        exit(EXIT_FAILURE);
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

  const char* deformationfile = argv[0];
  const char* jacoimagefile   = argv[1];

  // read deformation field
  deformation = ReadDeformationField(deformationfile);
  if (deformation == NULL) {
    cerr << "Failed to read deformation field from file " << deformationfile << endl;
    exit(1);
  }
  
  
  imageSize.x = deformation->region.nx;
  imageSize.y = deformation->region.ny;
  imageSize.z = deformation->region.nz;
  voxelSize.x = deformation->hdr.pixdim[1];
  voxelSize.y = deformation->hdr.pixdim[2];
  voxelSize.z = deformation->hdr.pixdim[3];
  if (imageSize.z==0) imageSize.z=1;
  if (voxelSize.z==0) voxelSize.z=1;
  
  
  Fvector3d*** DeformFld;
  std::vector<float> v;
  DeformFld=Fvector3dalloc3d(imageSize.x, imageSize.y, imageSize.z);
  for (k=0; k<imageSize.z; k++)
    for (i=0; i<imageSize.x; i++)
	  for (j=0; j<imageSize.y; j++) 
		{
		deformation->get(i, j, k, v);
		if (v.size() == 2) {
			DeformFld[k][i][j].x = v[0];
			DeformFld[k][i][j].y = v[1];
			DeformFld[k][i][j].z = 0;
		} else {
			DeformFld[k][i][j].x = v[0];
			DeformFld[k][i][j].y = v[1];
			DeformFld[k][i][j].z = v[2];
		}
      }

  // smooth deformation field if needed
  if (DeformationSmoothOrNot)
  {
  int defSmoothSizeXY=5;   if (MIN(imageSize.x, imageSize.y)<150)   defSmoothSizeXY=3;
  int defSmoothSizeZ=3;
  int u;
  float defSmoothKernelXY=(float)defSmoothSizeXY*0.45;
  float defSmoothKernelZ=(float)defSmoothSizeZ*0.35;
  Matrix **gaussianSmoothTemplate;
  gaussianSmoothTemplate = (Matrix**)malloc(sizeof(Matrix*)*defSmoothSizeZ);
  for (u=0;u<defSmoothSizeZ;u++)
      {
		CreateMatrix(&gaussianSmoothTemplate[u],defSmoothSizeXY, defSmoothSizeXY);
	  }
  generateGaussianSmoothTemplate3D(defSmoothSizeXY, defSmoothSizeZ, defSmoothKernelXY, defSmoothKernelZ, gaussianSmoothTemplate);
  smoothDeformationField(DeformFld, gaussianSmoothTemplate, imageSize, defSmoothSizeXY, defSmoothSizeZ);
  }	
	
  // calculate jacobian and save the result 
  printf("Calculating Determinants of Jacobian Matrix at each voxel, imageSize=(%d, %d, %d)...\n", imageSize.y, imageSize.x, imageSize.z);
  Image* jacoimage = new Image(imageSize.x, imageSize.y, imageSize.z, DT_FLOAT, 1, Image::FORMAT_DRAMMS);
  jacoimage->CopyTransform(deformation);
  jacoimage->filefmt       = deformation->filefmt;
  jacoimage->compress      = deformation->compress;
  float ***JacoDeterminant = jacoimage->img.fl;
  Jacobians( DeformFld, DRAMMSDeformationOrNot, imageSize.x, imageSize.y, imageSize.z, voxelSize, JacoDeterminant ) ;
  
  if (JacobianSmoothOrNot==YYES)
     {
	 int JacoSmoothWindowSizeXY=3;
     int JacoSmoothWindowSizeZ=3;
	 if ( MIN(imageSize.x, imageSize.y) > 150 )   JacoSmoothWindowSizeXY=5;
	 float JacoSmoothKernelXY = (float)JacoSmoothWindowSizeXY*0.45;
	 float JacoSmoothKernelZ = (float)JacoSmoothWindowSizeZ*0.25;
	 Matrix **gaussianFilter3D;
	 gaussianFilter3D = (Matrix**)malloc(sizeof(Matrix*)*JacoSmoothWindowSizeZ);
	 int t;
	 for (t=0;t<JacoSmoothWindowSizeZ;t++)
       CreateMatrix(&gaussianFilter3D[t],JacoSmoothWindowSizeXY, JacoSmoothWindowSizeXY);
	 generateGaussianSmoothTemplate3D(JacoSmoothWindowSizeXY, JacoSmoothWindowSizeZ, JacoSmoothKernelXY, JacoSmoothKernelZ, gaussianFilter3D);
	 smooth3DFloat(JacoDeterminant, gaussianFilter3D, imageSize, JacoSmoothWindowSizeXY, JacoSmoothWindowSizeZ);
	 }
	
  if (CutOffJacobianDetAtZeroOrNot==YYES)
	{
	// by default: this part is not used. We honestly keep all original Jaociban determinant values
	for (k=0;k<imageSize.z;k++)
	  for (i=0;i<imageSize.x;i++)
	    for (j=0;j<imageSize.y;j++)
		  {
		  if (JacoDeterminant[k][i][j]<0)
			JacoDeterminant[k][i][j]=0.0f;
		  }
	}
	
	
  if (LogJacobianDetOrNot==YYES)
	{
	for (k=0;k<imageSize.z;k++)
	  for (i=0;i<imageSize.x;i++)
	    for (j=0;j<imageSize.y;j++)
		  {
		  if (JacoDeterminant[k][i][j]<=0)
			JacoDeterminant[k][i][j]=-100.0;
		  else
			JacoDeterminant[k][i][j] = log(JacoDeterminant[k][i][j]);
		  }
	}
	
	
  WriteImage(jacoimagefile, jacoimage);
  printf("\nJacobian Determinant image has been saved as %s in float data type.\n\n", jacoimagefile);
  
  // clean up
  if (deformation) delete deformation;
  if (jacoimage)   delete jacoimage;
}




