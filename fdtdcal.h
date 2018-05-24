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
int initmur1(double cc,double cdt,double cdx);
// 全局变量
double det=1;//网格细度
int S0=0;//矩阵行列标识
int S1=0;//矩阵行列标识
int S2=0;//矩阵行列标识
int I,J,K;//矩阵长度极限
double * ex1;//步进电磁场矩阵
double * ey1;
double * ez1;
double * bx1;
double * by1;
double * bz1;

int flog=0;
// double * ex2;
// double * ey2;
// double * ez2;
// double * bx2;
// double * by2;
// double * bz2;
double * ca;//介质参数
double * cb;
double * cp;
double * cq;
// mur参数
double c;
double dt;
double dx;

