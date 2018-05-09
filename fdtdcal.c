#include "fdtdcal.h"
int initecal(double ex[],double ey[],double ez[],int strides[],int shapes[])
{
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
    bx1=bx;
    by1=by;
    bz1=bz;
    flog++;
    return 0;
}
int inite2cal(double ex[],double ey[],double ez[],int strides[],int shapes[])
{
    ex2=ex;
    ey2=ey;
    ez2=ez;
    flog++;
    return 0;
}
int initb2cal(double bx[],double by[],double bz[],int strides[],int shapes[])
{
    bx2=bx;
    by2=by;
    bz2=bz;
    flog++;
    return 0;
}
int cal(int n)
{
    double * tem;
    
    for (int ii = 0; ii < n; ii++)
    {
        
        for(int i = 1; i < I-1; i++)
        {
            
            for(int j = 1; i < J-1; j++)
            {
                
                for(int k = 1; k < K-1; k++)
                {
                    ex2[S0*i+S1*j+S2*k]=1*ex1[S0*i+S1*j+S2*k]+1*(1/det)*((bz1[S0*i+S1*j+S2*k]-bz1[S0*i+S1*(j-1)+S2*k])-by1[S0*i+S1*j+S2*k]+by1[S0*i+S1*j+S2*(k-1)]);
                    ey2[S0*i+S1*j+S2*k]=1*ey1[S0*i+S1*j+S2*k]+1*(1/det)*((bx1[S0*i+S1*j+S2*k]-bx1[S0*i+S1*j+S2*(k-1)])-bz1[S0*i+S1*j+S2*k]+bz1[S0*(i-1)+S1*j+S2*(k)]);
                    ez2[S0*i+S1*j+S2*k]=1*ez1[S0*i+S1*j+S2*k]+1*(1/det)*((by1[S0*i+S1*j+S2*k]-by1[S0*(i-1)+S1*j+S2*(k)])-bx1[S0*i+S1*j+S2*k]+by1[S0*(i)+S1*(j-1)+S2*(k)]);
                }
                
            }
            
        }
        

    }
    
    return 0;
}
int intarry(double adet)
{
    det=adet;
    return 0;
}

