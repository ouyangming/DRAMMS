
#pragma once
#ifndef _DRAMMS_CRES_H
#define _DRAMMS_CRES_H

#define TINY 1.0e-20;


#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
  float H;
  float K;
  float k1;
  float k2;
  float error;
} Skappa;

void nrerror(const char []);
double *vector(int,int);
int *ivector(int,int);
double *devecotr(int,int);
double **matrix(int,int,int,int);
double **dmatrix(int,int,int,int);
int **imatrix(int,int,int,int);
double **submatrix(double **,int,int,int,int,int,int);
void free_vector(double *,int,int);
void free_ivector(int *,int,int);
void free_dvector(double *,int,int);
void free_matrix(double **,int,int,int,int);
void free_dmatrix(double **,int,int,int,int);
void free_imatrix(int **,int,int,int,int);
void ludcmp(double **,int,int *,double *);
void lubksb(double **,int,int *,double *);


#ifdef __cplusplus
}
#endif


#endif // _DRAMMS_CRES_H
