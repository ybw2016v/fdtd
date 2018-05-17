#include "fdtdlib.h"
struct stdgridp calst();
// int calf(struct stdgrid * head)
// {

// }
struct stdgridp calst()//网格初始化。
{
    struct stdgridp * p1=NULL;//创建指针。
    p1=(struct stdgridp *)malloc(sizeof(struct stdgrid));//申请内存空间。
    p1->bx=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->by=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->bz=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->ex=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->ey=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->ez=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->ex->p1=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->ey->p1=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->ez->p1=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->bx->p1=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->by->p1=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    p1->bz->p1=(struct stdgrid *)malloc(sizeof(struct stdgrid));
    return *p1;
}
// struct name_t
// {
    
// };