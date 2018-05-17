#include "fdtdcal.h"
int initecal(double ex[],double ey[],double ez[],int strides[],int shapes[])
{
// 对电场矩阵进行初始化。
// ex,ey,ez为电场矩阵，strides为矩阵行列标识，shape为矩阵形状向量
    I=shapes[0];
    J=shapes[1];
    K=shapes[2];
    S0 = strides[0] / sizeof(double);
    S1 = strides[1] / sizeof(double);
    S2 = strides[2] / sizeof(double);
    ex1=ex;
    ey1=ey;
    ez1=ez;
    flog++;
    return 0;
}
int initbcal(double bx[],double by[],double bz[],int strides[],int shapes[])
{
    // 对磁场矩阵进行初始化
    bx1=bx;
    by1=by;
    bz1=bz;
    flog++;
    return 0;
}
int inite2cal(double ex[],double ey[],double ez[],int strides[],int shapes[])
{
    // 对电场矩阵进行初始化。
    ex2=ex;
    ey2=ey;
    ez2=ez;
    flog++;
    return 0;
}
int initb2cal(double bx[],double by[],double bz[],int strides[],int shapes[])
{
// 对磁场矩阵进行初始化
    bx2=bx;
    by2=by;
    bz2=bz;
    flog++;
    return 0;
}
int cal(int n)
{
// 计算核心函数
// n为计算次数
    double * tem;//暂存指针。
    // printf("%d %d %d\n",I,J,K);
    for (int ii = 0; ii < n; ii++)
    {
        // printf("%d \n",ii);
        for(int i = 2; i < I-2; i++)
        {
            // printf("%s \n","777");
            for(int j = 2; j < J-2; j++)
            {
                //printf("%s","++");
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    // 迭代核心
                    ex2[S0*i+S1*j+S2*k]=ca[0*i+S1*j+S2*k]*ex1[S0*i+S1*j+S2*k]+cb[0*i+S1*j+S2*k]*(1/det)*((bz1[S0*i+S1*j+S2*k]-bz1[S0*i+S1*(j-1)+S2*k])-by1[S0*i+S1*j+S2*k]+by1[S0*i+S1*j+S2*(k-1)]);
                    ey2[S0*i+S1*j+S2*k]=ca[0*i+S1*j+S2*k]*ey1[S0*i+S1*j+S2*k]+cb[0*i+S1*j+S2*k]*(1/det)*((bx1[S0*i+S1*j+S2*k]-bx1[S0*i+S1*j+S2*(k-1)])-bz1[S0*i+S1*j+S2*k]+bz1[S0*(i-1)+S1*j+S2*(k)]);
                    ez2[S0*i+S1*j+S2*k]=ca[0*i+S1*j+S2*k]*ez1[S0*i+S1*j+S2*k]+cb[0*i+S1*j+S2*k]*(1/det)*((by1[S0*i+S1*j+S2*k]-by1[S0*(i-1)+S1*j+S2*(k)])-bx1[S0*i+S1*j+S2*k]+by1[S0*(i)+S1*(j-1)+S2*(k)]);
                    bx2[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bx1[S0*i+S1*j+S2*k]-cq[0*i+S1*j+S2*k]*(1/det)*((ez1[S0*i+S1*j+S2*k]-ez1[S0*i+S1*(j-1)+S2*k])-ey1[S0*i+S1*j+S2*k]+ey1[S0*i+S1*j+S2*(k-1)]);
                    by2[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*by1[S0*i+S1*j+S2*k]-cq[0*i+S1*j+S2*k]*(1/det)*((ex1[S0*i+S1*j+S2*k]-ex1[S0*i+S1*j+S2*(k-1)])-ez1[S0*i+S1*j+S2*k]+ez1[S0*(i-1)+S1*j+S2*(k)]);
                    bz2[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bz1[S0*i+S1*j+S2*k]-cq[0*i+S1*j+S2*k]*(1/det)*((ey1[S0*i+S1*j+S2*k]-ey1[S0*(i-1)+S1*j+S2*(k)])-ex1[S0*i+S1*j+S2*k]+ey1[S0*(i)+S1*(j-1)+S2*(k)]);
                    // if (i==32&&j==32&&k==15)
                    // {
                    //     printf("%f \n",ez1[S0*i+S1*j+S2*k]);
                    // }
                }
                
            }
            
        }
    // printf("%s \n","66");
    // 指针交换。
    printf("%f \n",ez1[S0*32+S1*32+S2*15]);
    tem=ex1;
    ex1=ex2;
    ex2=tem;

    tem=ey1;
    ey1=ey2;
    ey2=tem;

    tem=ez1;
    ez1=ez2;
    ez2=tem;

    tem=bx1;
    bx1=bx2;
    bx2=tem;

    tem=by1;
    by1=by2;
    by2=tem;

    tem=bz1;
    bz1=bz2;
    bz2=tem;
    printf("%f \n",ez1[S0*32+S1*32+S2*15]);
    }
    
    return 0;
}
int intarry(double adet)
{
    // 初始化网格细度。
    det=adet;
    return 0;
}
int initc(double *cca,double * ccb,double * ccp,double * ccq)
{
    // 对电磁介质进行初始化；
    ca=cca;
    cb=ccb;
    cp=ccp;
    cq=ccq;
    return 0;
}

