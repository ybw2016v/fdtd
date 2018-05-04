#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
int main ()
{
    int i,j;
    #pragma omp parallel for
    for(i=1;i<=16;i++)
    {
        printf("%d\n",i);
    }
}
