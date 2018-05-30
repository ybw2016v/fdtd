#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// 函数声明
int initecal(double ex[],double ey[],double ez[],int strides[],int shapes[]);
int initbcal(double bx[],double by[],double bz[],int strides[],int shapes[]);
int inite2cal(double ex[],double ey[],double ez[],int strides[],int shapes[]);
int initb2cal(double bx[],double by[],double bz[],int strides[],int shapes[]);
int intarry(double adet);
int cal(int n);
int initc(double *cca,double * ccb,double * ccp,double * ccq);
int calmur1(int n);
int calmur2(int n);
int initmur1(double cc,double cdt,double cdx);



