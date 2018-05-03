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

