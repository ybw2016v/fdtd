#include "fdtdcal.h"

// OMP_SCHEDULE=8
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
double * ex2;
double * ey2;
double * ez2;
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
// int inite2cal(double ex[],double ey[],double ez[],int strides[],int shapes[])
// {
//     // 对电场矩阵进行初始化。
//     ex2=ex;
//     ey2=ey;
//     ez2=ez;
//     flog++;
//     return 0;
// }
// int initb2cal(double bx[],double by[],double bz[],int strides[],int shapes[])
// {
// // 对磁场矩阵进行初始化
//     bx2=bx;
//     by2=by;
//     bz2=bz;
//     flog++;
//     return 0;
// }
int cal(int n)
{
// 计算核心函数
// n为计算次数
    for (int ii = 0; ii < n; ii++)
    {
        for(int i = 2; i < I-2; i++)
        {
            for(int j = 2; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    ex1[S0*i+S1*j+S2*k]=ca[0*i+S1*j+S2*k]*ex1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((bz1[S0*i+S1*j+S2*(k-1)]-bz1[S0*i+S1*(j-1)+S2*(k-1)])-by1[S0*i+S1*(j-1)+S2*k]+by1[S0*i+S1*(j-1)+S2*(k-1)]);
                }
                
            }
            
        }
    for(int i = 2; i < I-2; i++)
        {
            for(int j = 2; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    ey1[S0*i+S1*j+S2*k]=ca[S0*i+S1*j+S2*k]*ey1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((bx1[S0*(i-1)+S1*j+S2*k]-bx1[S0*(i-1)+S1*j+S2*(k-1)])-bz1[S0*i+S1*j+S2*(k-1)]+bz1[S0*(i-1)+S1*j+S2*(k-1)]);
                }
                
            }
            
        }

    for(int i = 2; i < I-2; i++)
        {
            for(int j = 2; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    ez1[S0*i+S1*j+S2*k]=ca[S0*i+S1*j+S2*k]*ez1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((by1[S0*i+S1*(j-1)+S2*k]-by1[S0*(i-1)+S1*(j-1)+S2*(k)])-bx1[S0*(i-1)+S1*j+S2*k]+bx1[S0*(i-1)+S1*(j-1)+S2*(k)]);
                }
                
            }
            
        }
    for(int i = 2; i < I-2; i++)
        {
            for(int j = 2; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    bx1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bx1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ez1[S0*(i+1)+S1*(j+1)+S2*k]-ez1[S0*(i+1)+S1*j+S2*k])-ey1[S0*(i+1)+S1*j+S2*(k+1)]+ey1[S0*(i+1)+S1*j+S2*k]);
                }
                
            }
            
        }
    for(int i = 2; i < I-2; i++)
        {
            for(int j = 2; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    by1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*by1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ex1[S0*i+S1*(j+1)+S2*(k+1)]-ex1[S0*i+S1*(j+1)+S2*k])-ez1[S0*(i+1)+S1*(j+1)+S2*k]+ez1[S0*i+S1*(j+1)+S2*(k)]);
                }
                
            }
            
        }
    for(int i = 2; i < I-2; i++)
        {
            for(int j = 2; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 2; k < K-2; k++)
                {
                    // 迭代核心
                    bz1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bz1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ey1[S0*(i+1)+S1*j+S2*(k+1)]-ey1[S0*i+S1*j+S2*(k+1)])-ex1[S0*i+S1*(j+1)+S2*(k+1)]+ex1[S0*(i)+S1*(j)+S2*(k+1)]);
                }
                
            }
        }
    }

    
    return 0;
}
int intarry(double adet)
{
    // 初始化网格细度。
    det=1;
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
int initmur1(double cc,double cdt,double cdx)
{
    // mur边界初始化
    c=cc;
    dt=cdt;
    dx=cdx;
    return 0;
}
int calmur1(int n)
{
    double x1ey0[J][K];
    double xIey0[J][K];
    double x1ez0[J][K];
    double xIez0[J][K];


    double y1ex0[I][K];
    double yIex0[I][K];
    double y1ez0[I][K];
    double yIez0[I][K];

    double z1ex0[I][J];
    double zIex0[I][J];
    double z1ey0[I][J];
    double zIey0[I][J];
    //带有mur边界的计算函数
    for (int ii = 0; ii < n; ii++)
    {
        // 边界层的取值

        // x方向边界。
        for(int j = 1; j < J-2; j++)
    {
        int i=1;
        for (int k = 1; k < K-2; k++)
        {
            x1ey0[j][k]=ey1[S0*i+S1*j+S2*k];
            
        }
    }
        for(int j = 1; j < J-2; j++)
    {
        int i=I-1;
        for (int k = 1; k < K-2; k++)
        {
            xIey0[j][k]=ey1[S0*i+S1*j+S2*k];
            
        }
    }



        for(int j = 1; j < J-2; j++)
    {
        int i=1;
        for (int k = 1; k < K-2; k++)
        {
            x1ez0[j][k]=ez1[S0*i+S1*j+S2*k];
            
        }
    }
        for(int j = 1; j < J-2; j++)
    {
        int i=I-1;
        for (int k = 1; k < K-2; k++)
        {
            xIez0[j][k]=ez1[S0*i+S1*j+S2*k];
            
        }
    }



    // y方向
    for(int i = 1; i < I-2; i++)
    {
        int j=1;
        for (int k = 1; k < K-2; k++)
        {
            y1ex0[i][k]=ex1[S0*i+S1*j+S2*k];
            
        }
    }
        for(int i = 1; i < I-2; i++)
    {
        int j=J-2;
        for (int k = 1; k < K-2; k++)
        {
            yIex0[i][k]=ex1[S0*i+S1*j+S2*k];
            
        }
    }



        for(int i = 1; i < I-2; i++)
    {
        int j=1;
        for (int k = 1; k < K-2; k++)
        {
            y1ez0[i][k]=ez1[S0*i+S1*j+S2*k];
            
        }
    }
        for(int i = 1; i < I-2; i++)
    {
        int j=J-2;
        for (int k = 1; k < K-2; k++)
        {
            yIez0[i][k]=ez1[S0*i+S1*j+S2*k];
            
        }
    }
        // z方向


    for(int i = 1; i < I-2; i++)
    {
        int k=1;
        for (int j = 0; j < J-2; j++)
        {
            z1ex0[i][j]=ex1[S0*i+S1*j+S2*k];
        }
    }
    for(int i = 1; i < I-2; i++)
    {
        int k=K-3;
        for (int j = 0; j < J-2; j++)
        {
            z1ex0[i][j]=ex1[S0*i+S1*j+S2*k];
        }
    }
    for(int i = 1; i < I-2; i++)
    {
        int k=1;
        for (int j = 0; j < J-2; j++)
        {
            z1ey0[i][j]=ey1[S0*i+S1*j+S2*k];
        }
    }
    for(int i = 1; i < I-2; i++)
    {
        int k=K-2;
        for (int j = 1; j < J-2; j++)
        {
            z1ey0[i][j]=ey1[S0*i+S1*j+S2*k];
        }
    }



        // 计算核心
        for(int i = 1; i < I-2; i++)
        {
            for(int j = 1; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-2; k++)
                {
                    ex1[S0*i+S1*j+S2*k]=ca[0*i+S1*j+S2*k]*ex1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((bz1[S0*i+S1*j+S2*(k-1)]-bz1[S0*i+S1*(j-1)+S2*(k-1)])-by1[S0*i+S1*(j-1)+S2*k]+by1[S0*i+S1*(j-1)+S2*(k-1)]);
                }
                
            }
            
        }
    for(int i = 1; i < I-2; i++)
        {
            for(int j = 1; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-2; k++)
                {
                    ey1[S0*i+S1*j+S2*k]=ca[S0*i+S1*j+S2*k]*ey1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((bx1[S0*(i-1)+S1*j+S2*k]-bx1[S0*(i-1)+S1*j+S2*(k-1)])-bz1[S0*i+S1*j+S2*(k-1)]+bz1[S0*(i-1)+S1*j+S2*(k-1)]);
                }
                
            }
            
        }

    for(int i = 1; i < I-2; i++)
        {
            for(int j = 1; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-2; k++)
                {
                    ez1[S0*i+S1*j+S2*k]=ca[S0*i+S1*j+S2*k]*ez1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((by1[S0*i+S1*(j-1)+S2*k]-by1[S0*(i-1)+S1*(j-1)+S2*(k)])-bx1[S0*(i-1)+S1*j+S2*k]+bx1[S0*(i-1)+S1*(j-1)+S2*(k)]);
                }
                
            }
            
        }
    for(int i = 1; i < I-2; i++)
        {
            for(int j = 1; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-2; k++)
                {
                    bx1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bx1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ez1[S0*(i+1)+S1*(j+1)+S2*k]-ez1[S0*(i+1)+S1*j+S2*k])-ey1[S0*(i+1)+S1*j+S2*(k+1)]+ey1[S0*(i+1)+S1*j+S2*k]);
                }
                
            }
            
        }
    for(int i = 1; i < I-2; i++)
        {
            for(int j = 1; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-2; k++)
                {
                    by1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*by1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ex1[S0*i+S1*(j+1)+S2*(k+1)]-ex1[S0*i+S1*(j+1)+S2*k])-ez1[S0*(i+1)+S1*(j+1)+S2*k]+ez1[S0*i+S1*(j+1)+S2*(k)]);
                }
                
            }
            
        }
    for(int i = 1; i < I-2; i++)
        {
            for(int j = 1; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-2; k++)
                {
                    // 迭代核心
                    bz1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bz1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ey1[S0*(i+1)+S1*j+S2*(k+1)]-ey1[S0*i+S1*j+S2*(k+1)])-ex1[S0*i+S1*(j+1)+S2*(k+1)]+ex1[S0*(i)+S1*(j)+S2*(k+1)]);
                }
                
            }
        }
    




    for(int j = 0; j < J-2; j++)
    {
        int i=0;
        for (int k = 0; k < K-2; k++)
        {
            
            ey1[S0*i+S1*j+S2*k]=x1ey0[j][k]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*(i+1)+S1*j+S2*k]-ey1[S0*i+S1*j+S2*k]);
        }
    }
    for(int j = 0; j < J-2; j++)
    {
        int i=I-2;
        for (int k = 0; k < K-2; k++)
        {
            
            ey1[S0*i+S1*j+S2*k]=xIey0[j][k]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*(i+1)+S1*j+S2*k]-ey1[S0*i+S1*j+S2*k]);
        }
    }

        for(int j = 0; j < J-2; j++)
    {
        int i=0;
        for (int k = 0; k < K-2; k++)
        {
            
            ez1[S0*i+S1*j+S2*k]=x1ez0[j][k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*(i+1)+S1*j+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    }
    for(int j = 0; j < J-2; j++)
    {
        int i=I-2;
        for (int k = 0; k < K-2; k++)
        {
            
            ez1[S0*i+S1*j+S2*k]=xIez0[j][k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*(i-1)+S1*j+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    }


    // y方向

    for(int i = 0; i < I-2; i++)
    {
        int j=0;
        for (int k = 0; k < K-2; k++)
        {
            
            ex1[S0*i+S1*j+S2*k]=y1ex0[j][k]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*(i+0)+S1*(j+1)+S2*k]-ex1[S0*i+S1*j+S2*k]);
        }
    }
    for(int i = 2; i < I-2; i++)
    {
        int j=J-2;
        for (int k = 0; k < K-2; k++)
        {
            
            ex1[S0*i+S1*j+S2*k]=y1ex0[j][k]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*(i+0)+S1*(j-1)+S2*k]-ex1[S0*i+S1*j+S2*k]);
        }
    }

    for(int i = 0; i < I-2; i++)
    {
        int j=0;
        for (int k = 0; k < K-2; k++)
        {
            
            ez1[S0*i+S1*j+S2*k]=y1ez0[j][k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*(i+0)+S1*(j+1)+S2*(k+0)]-ez1[S0*i+S1*j+S2*k]);
        }
    }
    for(int i = 0; i < I-2; i++)
    {
        int j=J-2;
        for (int k = 0; k < K-2; k++)
        {
            
            ez1[S0*i+S1*j+S2*k]=y1ez0[j][k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*(i+0)+S1*(j-1)+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    
    }
    


// z方向

    for(int i = 0; i < I-2; i++)
    {
        int k=0;
        for (int j = 0; j < J-2; j++)
        {
            
            ex1[S0*i+S1*j+S2*k]=z1ex0[i][j]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*(i+0)+S1*j+S2*(k+1)]-ex1[S0*i+S1*j+S2*k]);
        }
    }

    for(int i = 0; i < I-2; i++)
    {
        int k=K-2;
        for (int j = 0; j < J-2; j++)
        {
            
            ex1[S0*i+S1*j+S2*k]=z1ex0[j][j]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*(i+0)+S1*j+S2*(k-1)]-ex1[S0*i+S1*j+S2*k]);
        }
    }
    for(int i = 0; i < I-2; i++)
    {
        int k=0;
        for (int j = 0; j < J-2; j++)
        {
            
            ey1[S0*i+S1*j+S2*k]=z1ey0[i][j]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*(i+0)+S1*j+S2*(k+1)]-ey1[S0*i+S1*j+S2*k]);
        }
    }

    for(int i = 0; i < I-2; i++)
    {
        int k=K-2;
        for (int j = 0; j < J-2; j++)
        {
            
            ey1[S0*i+S1*j+S2*k]=z1ey0[i][j]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*(i+0)+S1*j+S2*(k-1)]-ey1[S0*i+S1*j+S2*k]);
        }
    }



    }


return 0;
}

int initmur2(double ex[],double ey[],double ez[],int strides[],int shapes[])
{
    ex2=ex;
    ey2=ey;
    ez2=ez;
    return 0;
}

int calmur2(int n)
{
    // double 
    // double ex1p[S0*(I)+S1*J+S2*K];
    // double ey1p[S0*(I)+S1*J+S2*K];
    // double ez1p[S0*(I)+S1*J+S2*K];
    // ex1p=(double*) malloc(sizeof(double)*S0*(I)+S1*J+S2*K);
    // ey1p=(double*) malloc(sizeof(double)*S0*(I)+S1*J+S2*K);
    // ez1p=(double*) malloc(sizeof(double)*S0*(I)+S1*J+S2*K);
    double * ex1p=ex2;
    double * ey1p=ey2;
    double * ez1p=ez2;
    for (int i = 0;  i < I-1;i++)
    {
        for (int j = 0; j < J; j++)
        {
            #pragma omp parallel for
            for (int k = 0; k < K-1; k++)
            {
                ex1p[S0*i+S1*j+S2*k]=ex1[S0*i+S1*j+S2*k];
                ey1p[S0*i+S1*j+S2*k]=ey1[S0*i+S1*j+S2*k];
                ez1p[S0*i+S1*j+S2*k]=ez1[S0*i+S1*j+S2*k];
            }
        }
    }







    for (int ii = 0; ii < n; ii++)
    {
        for(int i = 0; i < I-1; i++)
        {
            for(int j = 1; j < J-1; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-1; k++)
                {
                    ex1[S0*i+S1*j+S2*k]=ca[0*i+S1*j+S2*k]*ex1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((bz1[S0*i+S1*j+S2*(k-1)]-bz1[S0*i+S1*(j-1)+S2*(k-1)])-by1[S0*i+S1*(j-1)+S2*k]+by1[S0*i+S1*(j-1)+S2*(k-1)]);
                }
                
            }
            
        }

    for (int i = 0,j=0; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ex1[S0*i+S1*j+S2*k]=ex1p[S0*i+S1*(j+1)+S2*k]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*i+S1*(j+1)+S2*k]-ex1p[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = 0,j=J-1; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ex1[S0*i+S1*j+S2*k]=ex1p[S0*i+S1*(j-1)+S2*k]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*i+S1*(j-1)+S2*k]-ex1p[S0*i+S1*j+S2*k]);
        }
    }

    for (int i = 0,k=0; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < J-1; j++)
        {
            ex1[S0*i+S1*j+S2*k]=ex1p[S0*i+S1*j+S2*(k+1)]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*i+S1*j+S2*(k+1)]-ex1p[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = 0,k=K-1; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < J-1; j++)
        {
            ex1[S0*i+S1*j+S2*k]=ex1p[S0*i+S1*j+S2*(k-1)]+((c*dt-dx)/(c*dt+dx))*(ex1[S0*i+S1*j+S2*(k-1)]-ex1p[S0*i+S1*j+S2*k]);
        }
    }




    for(int i = 1; i < I-1; i++)
        {
            for(int j = 0; j < J-1; j++)
            {
                #pragma omp parallel for
                for(int k = 1; k < K-1; k++)
                {
                    ey1[S0*i+S1*j+S2*k]=ca[S0*i+S1*j+S2*k]*ey1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((bx1[S0*(i-1)+S1*j+S2*k]-bx1[S0*(i-1)+S1*j+S2*(k-1)])-bz1[S0*i+S1*j+S2*(k-1)]+bz1[S0*(i-1)+S1*j+S2*(k-1)]);
                }
                
            }
            
        }

    for (int i = 0,j=0; j < J-1; j++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ey1[S0*i+S1*j+S2*k]=ey1p[S0*(i+1)+S1*j+S2*k]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*(i+1)+S1*j+S2*k]-ey1[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = I-1,j=0; j < J-1; j++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ey1[S0*i+S1*j+S2*k]=ey1p[S0*(i-1)+S1*j+S2*k]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*(i-1)+S1*j+S2*k]-ey1[S0*i+S1*j+S2*k]);
        }
    }

    for (int i = 0,k=0; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < J-1; j++)
        {
            ey1[S0*i+S1*j+S2*k]=ey1p[S0*i+S1*j+S2*(k+1)]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*i+S1*j+S2*(k+1)]-ey1[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = 0,k=K-1; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < J-1; j++)
        {
            ey1[S0*i+S1*j+S2*k]=ey1p[S0*i+S1*j+S2*(k-1)]+((c*dt-dx)/(c*dt+dx))*(ey1[S0*i+S1*j+S2*(k-1)]-ey1[S0*i+S1*j+S2*k]);
        }
    }






    for(int i = 1; i < I-1; i++)
        {
            for(int j = 1; j < J-1; j++)
            {
                #pragma omp parallel for
                for(int k = 0; k < K-1; k++)
                {
                    ez1[S0*i+S1*j+S2*k]=ca[S0*i+S1*j+S2*k]*ez1[S0*i+S1*j+S2*k]+cb[S0*i+S1*j+S2*k]*(1/det)*((by1[S0*i+S1*(j-1)+S2*k]-by1[S0*(i-1)+S1*(j-1)+S2*(k)])-bx1[S0*(i-1)+S1*j+S2*k]+bx1[S0*(i-1)+S1*(j-1)+S2*(k)]);
                }
                
            }
            
        }
    for (int i = 0,j=0; j < J-1; j++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ez1[S0*i+S1*j+S2*k]=ez1p[S0*(i+1)+S1*j+S2*k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*(i+1)+S1*j+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = I-1,j=0; j < J-1; j++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ez1[S0*i+S1*j+S2*k]=ez1p[S0*(i-1)+S1*j+S2*k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*(i-1)+S1*j+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = 0,j=0; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ez1[S0*i+S1*j+S2*k]=ez1p[S0*i+S1*(j+1)+S2*k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*i+S1*(j+1)+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    }
    for (int i = 0,j=J-1; i < I-1; i++)
    {
        #pragma omp parallel for
        for (int k = 0; k < K-1; k++)
        {
            ez1[S0*i+S1*j+S2*k]=ez1p[S0*i+S1*(j-1)+S2*k]+((c*dt-dx)/(c*dt+dx))*(ez1[S0*i+S1*(j-1)+S2*k]-ez1[S0*i+S1*j+S2*k]);
        }
    }







    for(int i = 0; i < I-2; i++)
        {
            for(int j = 0; j < J-1; j++)
            {
                #pragma omp parallel for
                for(int k = 0; k < K-2; k++)
                {
                    bx1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bx1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ez1[S0*(i+1)+S1*(j+1)+S2*k]-ez1[S0*(i+1)+S1*j+S2*k])-ey1[S0*(i+1)+S1*j+S2*(k+1)]+ey1[S0*(i+1)+S1*j+S2*k]);
                }
                
            }
            
        }
    for(int i = 0; i < I-2; i++)
        {
            for(int j = 0; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 0; k < K-1; k++)
                {
                    by1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*by1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ex1[S0*i+S1*(j+1)+S2*(k+1)]-ex1[S0*i+S1*(j+1)+S2*k])-ez1[S0*(i+1)+S1*(j+1)+S2*k]+ez1[S0*i+S1*(j+1)+S2*(k)]);
                }
                
            }
            
        }
    for(int i = 0; i < I-2; i++)
        {
            for(int j = 0; j < J-2; j++)
            {
                #pragma omp parallel for
                for(int k = 0; k < K-1; k++)
                {
                    // 迭代核心
                    bz1[S0*i+S1*j+S2*k]=cp[0*i+S1*j+S2*k]*bz1[S0*i+S1*j+S2*k]-cq[S0*i+S1*j+S2*k]*(1/det)*((ey1[S0*(i+1)+S1*j+S2*(k+1)]-ey1[S0*i+S1*j+S2*(k+1)])-ex1[S0*i+S1*(j+1)+S2*(k+1)]+ex1[S0*(i)+S1*(j)+S2*(k+1)]);
                }
                
            }
        }
    }
// free(ex1p);
// free(ey1p);
// free(ez1p);
// printf("%f\n",((c*dt-dx)/(c*dt+dx)));
return 0;
}