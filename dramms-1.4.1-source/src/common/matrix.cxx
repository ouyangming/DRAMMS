#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <math.h>
#include "matrix.h"

void nrerrorSHEN(const char *error_text)
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *dvectorSHEN(int nl, int nh)
{
	double *v;

	v = (double *) calloc((unsigned) (nh-nl+1), sizeof(double));
	if (!v) nrerrorSHEN("allocation failure in dvectorSHEN()");
	return v-nl;
}

void free_dvectorSHEN(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}


float *vectorSHEN(int nl, int nh)
{
	float *v;

	v = (float *) calloc((unsigned) (nh-nl+1), sizeof(float));
	if (!v) nrerrorSHEN("allocation failure in dvectorSHEN()");
	return v-nl;
}

void free_vectorSHEN(float *v, int nl, int nh)
{
	free((char*) (v+nl));
}

float **matrixSHEN(int nrl,int nrh,int ncl,int nch)
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerrorSHEN("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerrorSHEN("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_matrixSHEN(float **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


double **dmatrixSHEN(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m = (double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
	if (!m) nrerrorSHEN("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i] = (double *) calloc((unsigned) (nch-ncl+1), sizeof(double));
		if (!m[i]) nrerrorSHEN("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrixSHEN(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


double log2(double a)
{
	return log10(a)/log10(2.0);
}

void sort(double *Y, int *I, double *A, int length)
{
	int i, j;
	double max, *tmp;

	tmp = (double *) calloc(length, sizeof(double));

	for (i=0;i<length;i++) 
		tmp[i] = A[i];

	max = tmp[0];
	for (i=1;i<length;i++) {
		if (tmp[i] > max) 
			max = tmp[i];
	}

	max = fabs(10*max);

	for (i=0;i<length;i++) {
		Y[i] = tmp[0];
		I[i] = 0;
		for (j=1;j<length;j++) {
			if (tmp[j] < Y[i]) {
				Y[i] = tmp[j];
				I[i] = j;
			}
		}

		tmp[I[i]] = max;
	}

	free(tmp);
}

void minimun(double *Y, int *I, double *A, int length)
{
	int i, index;
	double min;

	min = A[0];
	index = 0;
	for (i=1;i<length;i++) 
		if (A[i] < min) {
			min = A[i];
			index = i;
		}

	*Y = min;
	*I = index;
}

void Mat_Abs(Matrix *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++) {
			if (A->data[h][w] < 0)
			    A->data[h][w] = -1.0*(A->data[h][w]);
		}
}

void Mat_Mean(double *mean, Matrix *A)
{
	int h, w;
	double tmp;

	tmp = 0.0;
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp += A->data[h][w];
		}
	}

	*mean = tmp/(double) (A->height*A->width);
}

void Mat_Variance(double *variance, Matrix *A)
{
	int h, w;
	double mean, tmp;

	Mat_Mean(&mean, A) ;

	tmp = 0.0;
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp += pow(A->data[h][w]-mean,2.0) ;
		}
	}

	*variance = sqrt(tmp/(double)(A->height*A->width));
}


void Mat_Vector(Matrix *A, float *a)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			a[h*A->width+w] = (float) A->data[h][w];
}

void Mat_Shift(Matrix *A, Matrix *B, int side)
{
	int h, w;

	for (h=side;h<B->height;h++)
		for (w=side;w<B->width;w++)
			A->data[h-side][w-side] = B->data[h][w];
		
	for (h=side;h<B->height;h++)
		for (w=0;w<side;w++)
			A->data[h-side][B->width-side+w] = B->data[h][w];

	for (h=0;h<side;h++)
		for (w=side;w<B->width;w++)
			A->data[B->height-side+h][w-side] = B->data[h][w];

	for (h=0;h<side;h++)
		for (w=0;w<side;w++)
			A->data[B->height-side+h][B->width-side+w] = B->data[h][w];
}

void Mat_Zeros(Matrix *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			A->data[h][w] = 0;
}

void Mat_Zeros_uc(uc_Matrix *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			A->data[h][w] = 0;
}

void Mat_Zeros_i(i_Matrix *A)
{
	int h, w;

	for (h=0;h<A->height;h++)
		for (w=0;w<A->width;w++)
			A->data[h][w] = 0;
}


void CreateMatrix(Matrix **M, int hei, int wid)
{
	int h;

	Matrix *tmp;

	tmp = (Matrix *) calloc(1, sizeof(Matrix)); 
	tmp->data = (double **) calloc(hei, sizeof(double *));

	if (!(tmp->data)) {
		nrerrorSHEN("allocation failure in CreateMatrix()");
		exit(1);
	}

	for (h=0; h<hei; h++) {
		tmp->data[h] = (double *) calloc(wid, sizeof(double));
		if (!(tmp->data[h])) {
			nrerrorSHEN("allocation failure in CreateMatrix()");
			exit(1);
		}
	}

	tmp->height = hei;
	tmp->width = wid;
	
	*M = tmp;
}


void FreeMatrix(Matrix *M)
{
	int h, hei = M->height;

	for (h=0; h<hei; h++) {
	     free(M->data[h]);
	}
	free(M->data);
	free(M);
}

void Create_i_Matrix(i_Matrix **M, int hei, int wid)
{
	int h;

	i_Matrix *tmp;

	tmp = (i_Matrix *) calloc(1, sizeof(i_Matrix));
	tmp->data = (int **) calloc(hei, sizeof(int *));
	if (!(tmp->data)) {
		nrerrorSHEN("allocation failure in Create_i_Matrix()");
		exit(1);
	}

	for (h=0; h<hei; h++) {
		tmp->data[h] = (int *) calloc(wid, sizeof(int));
		if (!(tmp->data[h])) {
			nrerrorSHEN("allocation failure in Create_i_Matrix()");
			exit(1);
		}
	}

	tmp->height = hei;
	tmp->width = wid;
	
	*M = tmp;
}

void Free_i_Matrix(i_Matrix *M)
{
	int h;

	for (h=0; h<M->height; h++) 
	     free(M->data[h]);
	free(M->data);
	free(M);
}

void Create_uc_Matrix(uc_Matrix **M, int hei, int wid)
{
	int h;

	uc_Matrix *tmp;

	tmp = (uc_Matrix *) calloc(1, sizeof(uc_Matrix));
	tmp->data = (unsigned char **) calloc(hei, sizeof(unsigned char *));
	if (!(tmp->data)) {
		nrerrorSHEN("allocation failure in Create_uc_Matrix()");
		exit(1);
	}

	for (h=0; h<hei; h++) {
		tmp->data[h] = (unsigned char *) calloc(wid, sizeof(unsigned char));
		if (!(tmp->data[h])) {
			nrerrorSHEN("allocation failure in Create_uc_Matrix()");
			exit(1);
		}
	}

	tmp->height = hei;
	tmp->width = wid;
	
	*M = tmp;
}

void Free_uc_Matrix(uc_Matrix *M)
{
	int h;

	for (h=0; h<M->height; h++) 
	     free(M->data[h]);
	free(M->data);
	free(M);
}


void Mat_FFT2(Matrix *Output_real, Matrix *Output_imag, Matrix *Input_real, Matrix *Input_imag)
{
	int xs, ys, i, j;
	double **R, **I, **Fr, **Fi;

	xs = Input_real->height;
	ys = Input_real->width;

	R  = dmatrixSHEN(1,xs,1,ys);
	I  = dmatrixSHEN(1,xs,1,ys);
	Fr = dmatrixSHEN(1,xs,1,ys);
	Fi = dmatrixSHEN(1,xs,1,ys);
			
	for (i=1;i<=Input_real->height;i++) 
	  for (j=1;j<=Input_real->width;j++) {
		R[i][j] = Input_real->data[i-1][j-1];
		I[i][j] = Input_imag->data[i-1][j-1];
	    }

	four2(Fr, Fi, R, I, xs, ys, 1);         /* 2-D FFT */

	for (i=1;i<=Input_real->height;i++) 
	    for (j=1;j<=Input_real->width;j++) {
		Output_real->data[i-1][j-1] = Fr[i][j];
		Output_imag->data[i-1][j-1] = Fi[i][j];
	    }

	free_dmatrixSHEN(R,1,xs,1,ys);
	free_dmatrixSHEN(I,1,xs,1,ys);   
	free_dmatrixSHEN(Fr,1,xs,1,ys);
	free_dmatrixSHEN(Fi,1,xs,1,ys);   
}

void Mat_IFFT2(Matrix *Output_real, Matrix *Output_imag, Matrix *Input_real, Matrix *Input_imag)
{
	int xs, ys, i, j;
	double **R, **I, **Fr, **Fi, NN;

	xs = Input_real->height;
	ys = Input_real->width;

	R  = dmatrixSHEN(1,xs,1,ys);
	I  = dmatrixSHEN(1,xs,1,ys);
	Fr = dmatrixSHEN(1,xs,1,ys);
	Fi = dmatrixSHEN(1,xs,1,ys);

	for (i=1;i<=Input_real->height;i++) 
	    for (j=1;j<=Input_real->width;j++) {
		R[i][j] = Input_real->data[i-1][j-1];
		I[i][j] = Input_imag->data[i-1][j-1];
	    }

	four2(Fr, Fi, R, I, xs, ys, -1);         /* 2-D IFFT */

	NN = (double) (xs*ys);

	for (i=1;i<=Input_real->height;i++) 
	    for (j=1;j<=Input_real->width;j++) {
		Output_real->data[i-1][j-1] = Fr[i][j]/NN;
		Output_imag->data[i-1][j-1] = Fi[i][j]/NN;
	    }

	free_dmatrixSHEN(R,1,xs,1,ys);
	free_dmatrixSHEN(I,1,xs,1,ys);   
	free_dmatrixSHEN(Fr,1,xs,1,ys);
	free_dmatrixSHEN(Fi,1,xs,1,ys);   
}

void four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign)
/************************************************************ 

   2-D fourier transform of data with real part stored in
   "rdata" and imaginary part in "idata" with size "rs" x
   "cs". The result is in "fftr" and "ffti". The isign is
   "isign" =  1 forward, and "isign" = -1 inverse 

*************************************************************/
{
	double **T, *tmp1, *tmp2;
	int i, j;

	tmp1 = dvectorSHEN(1,2*cs);
	tmp2 = dvectorSHEN(1,2*rs);
	T = dmatrixSHEN(1,2*rs,1,cs);

	for (i=1;i<=rs;i++) {
	    for (j=1;j<=cs;j++) {
		tmp1[j*2-1] = rdata[i][j];
		tmp1[j*2] = idata[i][j];
	    }
	    four1(tmp1, cs, isign);
	    for (j=1;j<=cs;j++) {
		T[i*2-1][j] = tmp1[j*2-1];
		T[i*2][j] = tmp1[j*2];
	    }
	}

	for (i=1;i<=cs;i++) {
	    for (j=1;j<=rs;j++) {
		tmp2[j*2-1] = T[j*2-1][i];
		tmp2[j*2] = T[j*2][i];
	    }
	    four1(tmp2,rs,isign);
	    for (j=1;j<=rs;j++) {
		fftr[j][i] = tmp2[j*2-1];
		ffti[j][i] = tmp2[j*2];
	    }
	}
	free_dvectorSHEN(tmp1, 1, 2*cs);
	free_dvectorSHEN(tmp2, 1, 2*rs);
	free_dmatrixSHEN(T, 1, 2*rs, 1, cs); 
}

void four1(double *data, int nn, int isign)
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) {
		istep = 2*mmax;
		theta = 6.28318530717959/(isign*mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j = i+mmax;
				tempr = wr*data[j]-wi*data[j+1];
				tempi = wr*data[j+1]+wi*data[j];
				data[j] = data[i]-tempr;
				data[j+1] = data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp=wr)*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}
		mmax = istep;
	}
}

void Mat_Copy(Matrix *A, Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end)
{
	int i, j, h, w, h_done, w_done;

	if ((h_target >= 0)&&(h_target < A->height)&&(w_target >= 0)&&(w_target < A->width)) {
		if ((h_begin >= 0)&&(h_begin < B->height)&&(w_begin >= 0)&&(w_begin < B->width)) {
			h = h_end-h_begin+1;
			w = w_end-w_begin+1;
			if ((h >= 1)&&(w >= 1)) {
				h_done = h_target+h-1;
				w_done = w_target+w-1;
				if ((h_done < A->height)&&(w_done < A->width)) {
					for (i=0;i<h;i++) {
						for (j=0;j<w;j++) {
							A->data[i+h_target][j+w_target] = B->data[i+h_begin][j+w_begin];
						}
					}
				}
			}
		}
	}
	else {
		printf("matrix dimension error!\n");
		exit(1);
	}
}

void Mat_uc_Copy(uc_Matrix *A, uc_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end)
{
	int i, j, h, w, h_done, w_done;

	if ((h_target >= 0)&&(h_target < A->height)&&(w_target >= 0)&&(w_target < A->width)) {
		if ((h_begin >= 0)&&(h_begin < B->height)&&(w_begin >= 0)&&(w_begin < B->width)) {
			h = h_end-h_begin+1;
			w = w_end-w_begin+1;
			if ((h >= 1)&&(w >= 1)) {
				h_done = h_target+h-1;
				w_done = w_target+w-1;
				if ((h_done < A->height)&&(w_done < A->width)) {
					for (i=0;i<h;i++) {
						for (j=0;j<w;j++) {
							A->data[i+h_target][j+w_target] = B->data[i+h_begin][j+w_begin];
						}
					}
				}
			}
		}
	}
	else {
		printf("matrix dimension error!\n");
		exit(1);
	}
}
	
void Mat_i_Copy(i_Matrix *A, i_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end)
{
	int i, j, h, w, h_done, w_done;

	if ((h_target >= 0)&&(h_target < A->height)&&(w_target >= 0)&&(w_target < A->width)) {
		if ((h_begin >= 0)&&(h_begin < B->height)&&(w_begin >= 0)&&(w_begin < B->width)) {
			h = h_end-h_begin+1;
			w = w_end-w_begin+1;
			if ((h >= 1)&&(w >= 1)) {
				h_done = h_target+h-1;
				w_done = w_target+w-1;
				if ((h_done < A->height)&&(w_done < A->width)) {
					for (i=0;i<h;i++) {
						for (j=0;j<w;j++) {
							A->data[i+h_target][j+w_target] = B->data[i+h_begin][j+w_begin];
						}
					}
				}
			}
		}
	}
	else {
		printf("matrix dimension error!\n");
		exit(1);
	}
}
				
void Mat_Product(Matrix *A, Matrix *B, Matrix *C)
{
	int h, w;


	if(A->height!=B->height || A->height!=C->height ||
	   A->width!=B->width || A->width!=C->width )
	  nrerrorSHEN("Mat_Substract fail!");


	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = B->data[h][w]*C->data[h][w];
		}
	}
}

void Mat_Sum(Matrix *A, Matrix *B, Matrix *C)
{
	int h, w;


	if(A->height!=B->height || A->height!=C->height ||
	   A->width!=B->width || A->width!=C->width )
	  nrerrorSHEN("Mat_Substract fail!");


	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = B->data[h][w]+C->data[h][w];
		}
	}
}

void Mat_Substract(Matrix *A, Matrix *B, Matrix *C)
{
	int h, w;

	if(A->height!=B->height || A->height!=C->height ||
	   A->width!=B->width || A->width!=C->width )
	  nrerrorSHEN("Mat_Substract fail!");


	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = B->data[h][w]-C->data[h][w];
		}
	}
}

void Mat_Fliplr(Matrix *A)
{
	Matrix *tmp;
	int h, w;

	CreateMatrix(&tmp, A->height, A->width);
	
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[h][(A->width)-w-1];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	FreeMatrix(tmp);
}

void Mat_Flipud(Matrix *A)
{
	Matrix *tmp;
	int h, w;

	CreateMatrix(&tmp, A->height, A->width);
	
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[(A->height)-h-1][w];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	FreeMatrix(tmp);
}


void Mat_uc_Fliplr(uc_Matrix *A)
{
	uc_Matrix *tmp;
	int h, w;

	Create_uc_Matrix(&tmp, A->height, A->width);
	
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[h][(A->width)-w-1];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	Free_uc_Matrix(tmp);
}

void Mat_uc_Flipud(uc_Matrix *A)
{
	uc_Matrix *tmp;
	int h, w;

	Create_uc_Matrix(&tmp, A->height, A->width);
	
	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			tmp->data[h][w] = A->data[(A->height)-h-1][w];
		}
	}

	for (h=0;h<A->height;h++) {
		for (w=0;w<A->width;w++) {
			A->data[h][w] = tmp->data[h][w];
		}
	}

	Free_uc_Matrix(tmp);
}


/*add by SHEN when in JHU*/
int *ivectorSHEN(long nl, long nh)
{ 
	int *v; 

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int))); 
	if (!v) nrerrorSHEN("allocation failure in ivector()"); 
	return v-nl+NR_END; 
} 
void free_ivectorSHEN(int *v, long nl, long nh) 
{ 
	free((FREE_ARG) (v+nl-NR_END)); 
} 


int gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv;

	indxc=ivectorSHEN(1,n);
	indxr=ivectorSHEN(1,n);
	ipiv=ivectorSHEN(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1)//nrerrorSHEN("gaussj: Singular Matrix-1");
					  {
					    printf("gaussj: Singular Matrix-1");
					    return -1;
					  }
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) 
		  {
		    printf("gaussj: Singular Matrix-2");
		    return -1;
		  }//nrerrorSHEN("gaussj: Singular Matrix-2");

		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivectorSHEN(ipiv,1,n);
	free_ivectorSHEN(indxr,1,n);
	free_ivectorSHEN(indxc,1,n);

	return 1;
}

float pythag(float a, float b)
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmpSHEN(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vectorSHEN(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerrorSHEN("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vectorSHEN(rv1,1,n);
}

void Mat_A_equal_BxC(Matrix *A, Matrix *B, Matrix *C)
{
	int h, w, k;
	double sum ;

	if(B->width!=C->height || A->height!=B->height || A->width!=C->width ) 
	  {
	    printf("Matrix operation error!\n"); 
	    exit(1);
	  }


	for (h=0; h<A->height; h++)
	  {
	    for (w=0; w<A->width; w++)
	      {
		sum = 0 ;
		for(k=0; k<B->width; k++)
		  sum += (B->data[h][k] * C->data[k][w]) ;

		A->data[h][w] = sum ;
	      }
	  }

}

void Mat_Print(Matrix *A)
{
  int k,l ;

  printf("\n");

  for(k=0; k<A->height; k++) 
    {
      for (l=0; l<A->width; l++) 
	printf("%12.6f ", A->data[k][l]) ;
      printf("\n");
    }
}


void Mat_Calculate_EigenVectors_EigenValues(Matrix *C, float *EigenValue, Matrix *EigenVector, int PRNorNOT)
{
  int   j, k, l,   np, mp,  position ;
  float tmp, max ;

  /* for getting eigenVectors (u) and eigenValues (w)*/
  float *w,**u,**v;


  /* Calculate eigenvectors and eigenvalues */
   mp = C->height ;  /* size */
   np = C->width ;
   /* apply memory */
   /* Notice here the beginning index for u (v or w) is 1, not 0.  */
   w = vectorSHEN(1,np);
   u = matrixSHEN(1,mp,1,np);
   v = matrixSHEN(1,np,1,np);
   
   /* set matrices: u <- C */
   for(k=1; k<=mp; k++)
     for(l=1; l<=np; l++) 
       u[k][l] = C->data[k-1][l-1] ;


  if(PRNorNOT==TRUE)
  {  
    printf("\nCheck product against original matrix:\n");
    printf("Original matrix:\n");
    for (k=1;k<=mp;k++) 
      {
	for (l=1;l<=np;l++)
	  printf("%12.3f ",u[k][l]);
	printf("\n");
      }
  }		

  /* perform decomposition */
  svdcmpSHEN(u,mp,np,w,v);


  if(PRNorNOT==TRUE)
  {
    printf("Product u*w*(v-transpose):\n");
    for (k=1;k<=mp;k++) 
      {
	for (l=1;l<=np;l++) 
	  {
	    tmp=0.0;
	    for (j=1;j<=np;j++)
	      tmp += u[k][j]*w[j]*v[l][j];
	    printf("%12.3f ",tmp) ;
	  }
	printf("\n");
      }
      
    /* write results */
    printf("Decomposition matrices:\n");
    printf("Matrix u\n");
    for (k=1; k<=mp; k++)
      {
	for (l=1; l<=np; l++)
	  printf("%12.3f ", u[k][l]);
	printf("\n");
      }

    printf("Diagonal of matrix w\n");
    for (k=1;k<=np;k++)
      printf("%12.3f ", w[k]);
    
    printf("\nMatrix v-transpose\n");
    for (k=1;k<=np;k++) 
      {
	for (l=1;l<=np;l++)
	  printf("%12.3f ", v[l][k]);
	printf("\n");
      }
  }


  /* record for returning */
   for(k=1; k<=np; k++)
     {
       max = -100000.0 ;
      for(j=1; j<=np; j++)
       {
	 if( w[j]>max ) 
	   {
	     max = w[j] ;
	     position = j ;
	   }
       }  
      /*eigenvalue*/
      EigenValue[k-1] = w[position] ;
      /*eigenvector*/
      for(l=1;l<=np;l++) 
	  EigenVector->data[l-1][k-1] = u[l][position] ;
      w[position] = -1000000.0 ;
     }


   free_vectorSHEN(w,1,np); 
   free_matrixSHEN(u,1,mp,1,np);
   free_matrixSHEN(v,1,np,1,np);
}


int Mat_Inverse(Matrix *A, Matrix *B) /*A in, B out*/
{
  int    i, j, m, n1, n2 ;
  float  **a,**b;
  Matrix *At, *AtA, *AtAInv ;

  n1 = A->height ;
  n2 = A->width ;
  m = 1 ;

  if(n1==n2)
    {
      if( B->height != n1 || B->width != n2 )
	nrerrorSHEN("marix size not matching") ;

      /* apply memory */
      /* Notice here the beginning index for u (v or w) is 1, not 0.  */
      a = matrixSHEN(1,n1,1,n2);
      b = matrixSHEN(1,n1,1,m);
      
      /* set matrices: a <- A */
      for(i=1; i<=n1; i++)
	for(j=1; j<=n2; j++) 
	  a[i][j] = A->data[i-1][j-1] ;

      /*cheat b*/
     for(i=1; i<=n1; i++)
       for(j=1; j<=m; j++) 
	  b[i][j] = A->data[i-1][j-1] ;

     if (-1==gaussj(a, n1, b, m))
       return -1; /* get inverse of a */

      /* set matrices: A <- a */
      for(i=1; i<=n1; i++)
	for(j=1; j<=n2; j++) 
	  B->data[i-1][j-1] = a[i][j] ;

     free_matrixSHEN(a,1,n1,1,n2);
     free_matrixSHEN(b,1,n1,1,m);
    }
  else /* n1 != n2*/
    {     
      if( B->height != n2 || B->width != n1 )
	nrerrorSHEN("marix size not matching") ;

      /*printf("n1=%d n2=%d\n", n1, n2) ;*/

      CreateMatrix(&At,  n2, n1);
      CreateMatrix(&AtA, n2, n2);
      CreateMatrix(&AtAInv, n2, n2); /* used as the inverse matrix of C*/

      /* B=A(t)*/
      for(i=0; i<n2; i++)
	for(j=0; j<n1; j++) 
	  At->data[i][j] = A->data[j][i] ;

printf("At\n") ; Mat_Print(At) ;
printf("A\n") ; Mat_Print(A) ;

      Mat_A_equal_BxC(AtA, At,  A) ; /*A=BxC*/
printf("AtA\n") ; Mat_Print(AtA) ;

      Mat_Inverse(AtA, AtAInv) ;     /*A in, B out*/
printf("AtAI\n") ; Mat_Print(AtAInv) ;

      Mat_A_equal_BxC(B, AtAInv,  At) ; /*A=BxC*/


      FreeMatrix(At) ;  
      FreeMatrix(AtA) ;
      FreeMatrix(AtAInv) ;
    }
  return 1;
}


void Mat_times_Vector(float *Vout, Matrix *A, float *Vin) 
{
  int    i, j, n1, n2 ;
  double sum ;

  n1 = A->height ;
  n2 = A->width ;

  /* set matrices: a <- A */
  for(i=0; i<n1; i++)
    {
      sum = 0 ;
      for(j=0; j<n2; j++) 
	sum += A->data[i][j]*Vin[j] ;
      Vout[i] = sum ;
    }
}


void vector_Print(float *v, int size)
{
  int k ;

  printf("\n");

  for(k=0; k<size; k++) 
    printf("%12.3f ", v[k]) ;
}


void Mat_EqualCopy(Matrix *A, Matrix *B)
{
	int i, j, h, w;

	if ( A->height==B->height && A->width==B->width )
	  {
	    h = A->height ;
	    w = A->width  ;
	    
	    for (i=0; i<h; i++) 
	      for (j=0; j<w; j++) 
		A->data[i][j] = B->data[i][j];
	  }
	else 
	  {
	    printf("Matrix copy: matrix dimension error!\n");
	    exit(1);
	  }
}



#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


/* June 2001  */
ImgAttribute ****ImgAttributealloc4d(int i_size,int j_size,int k_size, int t_size)
{
  ImgAttribute ****array;
  int i,k, t;

  array=(ImgAttribute ****) calloc(t_size,sizeof(ImgAttribute ***));

  for(t=0;t<t_size;t++)
    array[t]=(ImgAttribute ***) calloc(k_size,sizeof(ImgAttribute **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(ImgAttribute **) calloc(i_size,sizeof(ImgAttribute *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(ImgAttribute *) calloc(j_size,sizeof(ImgAttribute ));
	
  return(array);
}

ImgAttribute ***ImgAttributealloc3d(int i_size,int j_size,int k_size)
{
  ImgAttribute ***array;
  int i,k;

  array=(ImgAttribute ***) calloc(k_size,sizeof(ImgAttribute **));

  for(k=0;k<k_size;k++)
    array[k]=(ImgAttribute **) calloc(i_size,sizeof(ImgAttribute *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(ImgAttribute *) calloc(j_size,sizeof(ImgAttribute ));
	
  return(array);
}

ImgAttribute *ImgAttributealloc1d(int k_size)
{
  ImgAttribute *array;
  
  array=(ImgAttribute *) calloc(k_size,sizeof(ImgAttribute));

  return(array);
}

/*free*/
void ImgAttributefree4d(ImgAttribute ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}

void ImgAttributefree3d(ImgAttribute ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}


/* Nov 2001, for mapping volumetric head images */
HeadImgAttribute ***HeadImgAttributealloc3d(int i_size,int j_size,int k_size)
{
  HeadImgAttribute ***array;
  int i,k;

  array=(HeadImgAttribute ***) calloc(k_size,sizeof(HeadImgAttribute **));

  for(k=0;k<k_size;k++)
    array[k]=(HeadImgAttribute **) calloc(i_size,sizeof(HeadImgAttribute *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(HeadImgAttribute *) calloc(j_size,sizeof(HeadImgAttribute ));
	
  return(array);
}

HeadImgAttribute *HeadImgAttributealloc1d(int k_size)
{
  HeadImgAttribute *array;
  
  array=(HeadImgAttribute *) calloc(k_size,sizeof(HeadImgAttribute));

  return(array);
}

/*free*/
void HeadImgAttributefree3d(HeadImgAttribute ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}



/* Feb 2002, for warping DTI */
DTIattribute ***DTIattributealloc3d(int i_size,int j_size,int k_size)
{
  DTIattribute ***array;
  int i,k;

  array=(DTIattribute ***) calloc(k_size,sizeof(DTIattribute **));

  for(k=0;k<k_size;k++)
    array[k]=(DTIattribute **) calloc(i_size,sizeof(DTIattribute *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(DTIattribute *) calloc(j_size,sizeof(DTIattribute ));
	
  return(array);
}

DTIattribute *DTIattributealloc1d(int k_size)
{
  DTIattribute *array;
  
  array=(DTIattribute *) calloc(k_size,sizeof(DTIattribute));

  return(array);
}


/*free*/
void DTIattributefree3d(DTIattribute ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

//Yiqiang added Aug. 2005
void Mat_Transpose(Matrix* A,Matrix*B)
{
  int h,w;

  if ((A->height!=B->width)||(A->width!=B->height))
    {    
      printf("Matrix dimension doesn't match!\n");
      return;
    }
  
  for (h=0;h<A->height;h++)
    for (w=0;w<A->width;w++)
      B->data[w][h]=A->data[h][w];
  
}

void SaveMatrix(char textX[80],Matrix* A)
{
  FILE *f;
  int i;
  int hei,wei;

  if( (f=fopen(textX,"wb"))==NULL )  	{ printf("Fail\n") ; return ; }
  fseek(f,0L,SEEK_SET);

  hei=A->height;
  wei=A->width;

  fwrite(&hei,sizeof(int),1,f);
  fwrite(&wei,sizeof(int),1,f);
  
  for (i=0;i<hei;i++)      
    fwrite(A->data[i],sizeof(double),wei,f);
       
  fclose(f);
}

void OpenMatrix1(char textX[80],int *hei,int *wei)
{
  FILE *f;
  
  if( (f=fopen(textX,"rb"))==NULL )  	{ printf("Fail\n") ; return ; }
  fseek(f,0L,SEEK_SET);

  fread(hei,sizeof(int),1,f);
  fread(wei,sizeof(int),1,f);
       
  fclose(f);
}

void OpenMatrix2(char textX[80],Matrix* A)
{
  FILE *f;
  int i;
  int hei,wei;

  if( (f=fopen(textX,"rb"))==NULL )  	{ printf("Fail\n") ; return ; }
  fseek(f,2*sizeof(int),SEEK_SET);

  hei=A->height;
  wei=A->width;

  for (i=0;i<hei;i++)
    fread(A->data[i],sizeof(double),wei,f);
  
  fclose(f);
}

void Mat_Trace(double* trace,Matrix* A)
{

  int i;
  double sum;

  if (A->height!=A->width)
    printf("Matrix dimension error in function Mat_Trace!\n");

  sum=0;

  for (i=0;i<A->height;i++)
    sum+=A->data[i][i];

  *trace=sum;

}
