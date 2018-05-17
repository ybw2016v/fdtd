#include <stdio.h>
#include <stdlib.h>
#define N 60
// 网格尺寸
struct stdgrid 
{
    float data[N][N][N];
    int flog;
    struct stdgrid * p1;
    struct stdgrid * p2;
};
struct stdgridp
{
    struct stdgrid * ex;
    struct stdgrid * ey;
    struct stdgrid * ez;
    struct stdgrid * bx;
    struct stdgrid * by;
    struct stdgrid * bz;    
};

