#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int initecal(double ex[],double ey[],double ez[],int strides[],int shapes[]);
int initbcal(double bx[],double by[],double bz[],int strides[],int shapes[]);
int inite2cal(double ex[],double ey[],double ez[],int strides[],int shapes[]);
int initb2cal(double bx[],double by[],double bz[],int strides[],int shapes[]);
int intarry(double adet);
int cal(int n);
double det=2;
int S0=0;
int S1=0;
int S2=0;
double * ex1;
double * ey1;
double * ez1;
double * bx1;
double * by1;
double * bz1;
int I,J,K;
int flog=0;
double * ex2;
double * ey2;
double * ez2;
double * bx2;
double * by2;
double * bz2;