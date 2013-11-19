#ifndef MATRIX_H
#define MATRIX_H

#include <math.h> // required to avoid error regarding log2()

#ifdef __cplusplus
 extern "C" {
 #endif


/* Feb 2002, for warping DTI */
typedef struct 
{  
  float v1 ;
  float v2 ;
  float v3 ;
  float v4 ;
  float v5 ;
  float v6 ;
} DTIattribute ;  



/* Nov 2001, for Head brain image */
typedef struct 
{  
  unsigned char Edge;  
  unsigned char Tiss;  
  unsigned char Geom;   
  unsigned char BGvlm;  
  unsigned char CSFvlm;  
  unsigned char VNvlm;  
} HeadImgAttribute ;  



/* June 2001, for skull-stripped brain image */
typedef struct
{
  unsigned char Edge;
  unsigned char Tiss;
  unsigned char Geom; 
  unsigned char VNvlm;
  unsigned char CSFBG;
} ImgAttribute ;



typedef struct MatrixStruct 
{
	double **data;
	int height, width;
} Matrix;

typedef struct uc_MatrixStruct 
{
	unsigned char **data;
	int height, width;
} uc_Matrix;

typedef struct i_MatrixStruct 
{
	int **data;
	int height, width;
} i_Matrix;

#define TRUE     1
#define FREE_ARG char*
#define NR_END 1
#define PI 3.141592653589793115997963468544185161590576171875

static float tempr;
namespace { inline float dummy_use_of_tempr() { return tempr; } }
#define SWAP(a,b) {tempr=(a);(a)=(b);(b)=tempr;}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1,maxarg2;
namespace { inline float dummy_use_of_maxarg() { return maxarg1 + maxarg2; } }
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static float minarg1,minarg2;
namespace { inline float dummy_use_of_minarg() { return minarg1 + minarg2; } }
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static int iminarg1,iminarg2;
namespace { inline int dummy_use_of_iminarg() { return iminarg1 + iminarg2; } }
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?  (iminarg1) : (iminarg2))

static float sqrarg;
namespace { inline float dummy_use_of_sqrarg() { return sqrarg; } }
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)




void   nrerrorSHEN(const char *error_text);
double *dvectorSHEN(int nl, int nh);
void   free_dvectorSHEN(double *v, int nl, int nh);
float  *vectorSHEN(int nl, int nh);
void   free_vectorSHEN(float *v, int nl, int nh);
float  **matrixSHEN(int nrl,int nrh,int ncl,int nch);
void   free_matrixSHEN(float **m,int nrl,int nrh,int ncl,int nch);
double **dmatrixSHEN(int nrl, int nrh, int ncl, int nch);
void   free_dmatrixSHEN(double **m, int nrl, int nrh, int ncl, int nch);
double log2(double a);
void   sort(double *Y, int *I, double *A, int length);
void   minimun(double *Y, int *I, double *A, int length);
void   Mat_Abs(Matrix *A);
void   Mat_Mean(double *mean, Matrix *A);
void   Mat_Variance(double *variance, Matrix *A) ;
void   Mat_Vector(Matrix *A, float *a);
void   Mat_Shift(Matrix *A, Matrix *B, int side);
void   Mat_Zeros(Matrix *A);
void   Mat_Zeros_uc(uc_Matrix *A);
void   Mat_Zeros_i(i_Matrix *A);
void   CreateMatrix(Matrix **M, int hei, int wid);
void   FreeMatrix(Matrix *M);
void   Create_i_Matrix(i_Matrix **M, int hei, int wid);
void   Free_i_Matrix(i_Matrix *M);
void   Create_uc_Matrix(uc_Matrix **M, int hei, int wid);
void   Free_uc_Matrix(uc_Matrix *M);
void   Mat_FFT2(Matrix *Output_real, Matrix *Output_imag, Matrix *Input_real, Matrix *Input_imag);
void   Mat_IFFT2(Matrix *Output_real, Matrix *Output_imag, Matrix *Input_real, Matrix *Input_imag);
void   four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign);
void   four1(double *data, int nn, int isign);
void   Mat_Copy(Matrix *A, Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end,int w_end);
void   Mat_uc_Copy(uc_Matrix *A, uc_Matrix *B, int h_target, int w_target, int h_begin, int w_begin,int h_end, int w_end);
void   Mat_i_Copy(i_Matrix *A, i_Matrix *B, int h_target, int w_target, int h_begin, int w_begin, int h_end, int w_end);
void   Mat_Product(Matrix *A, Matrix *B, Matrix *C);
void   Mat_Sum(Matrix *A, Matrix *B, Matrix *C);
void   Mat_Substract(Matrix *A, Matrix *B, Matrix *C);
void   Mat_Fliplr(Matrix *A);
void   Mat_Flipud(Matrix *A);
void   Mat_uc_Fliplr(uc_Matrix *A);
void   Mat_uc_Flipud(uc_Matrix *A);

/*by SHEN in JHU*/
int    *ivectorSHEN(long nl, long nh) ;
void   free_ivectorSHEN(int *v, long nl, long nh) ;
int   Mat_Inverse(Matrix *A, Matrix *B) ; /*A in, B out*/
void   Mat_times_Vector(float *Vout, Matrix *A, float *Vin) ;
void   Mat_A_equal_BxC(Matrix *A, Matrix *B, Matrix *C) ;
void   Mat_EqualCopy(Matrix *A, Matrix *B) ;
void   Mat_Print(Matrix *A) ;
void   Mat_Calculate_EigenVectors_EigenValues(Matrix *C, float *EigenValue, Matrix *EigenVector, int PRNorNOT) ;
void   svdcmp(float **a, int m, int n, float w[], float **v) ;
void   vector_Print(float *v, int size) ;
float  gasdev(long *idum) ;

/* June 2001 */
ImgAttribute ****ImgAttributealloc4d(int i_size,int j_size,int k_size, int t_size);
ImgAttribute ***ImgAttributealloc3d(int i_size,int j_size,int k_size) ;
ImgAttribute *ImgAttributealloc1d(int k_size) ;
void ImgAttributefree4d(ImgAttribute ****array,int t_size,int k_size,int i_size);
void ImgAttributefree3d(ImgAttribute ***array,int k_size,int i_size) ;

/* June 2001 */
HeadImgAttribute ***HeadImgAttributealloc3d(int i_size,int j_size,int k_size) ;
HeadImgAttribute *HeadImgAttributealloc1d(int k_size) ;
void HeadImgAttributefree3d(HeadImgAttribute ***array,int k_size,int i_size) ;

/* Feb 2002 */
DTIattribute ***DTIattributealloc3d(int i_size,int j_size,int k_size) ;
DTIattribute *DTIattributealloc1d(int k_size) ;
void DTIattributefree3d(DTIattribute ***array,int k_size,int i_size) ;

//Yiqiang added on Aug. 2005
void Mat_Transpose(Matrix* A,Matrix*B);

void Mat_Trace(double* trace,Matrix* A);
void SaveMatrix(char textX[80],Matrix* A);
void OpenMatrix1(char textX[80],int *hei,int *wei);
void OpenMatrix2(char textX[80],Matrix* A);

#ifdef __cplusplus
 }
 #endif

#endif
