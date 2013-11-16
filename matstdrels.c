#include "matrix.h"


int gen_gt(mtype a)
{
    if(a>0)
        return 1;
    else
        return 0;
}

int gen_lt(mtype a)
{
    if(a<0)
        return 1;
    else
        return 0;
}

int gen_eq(mtype a)
{
    if(a==0)
        return 1;
    else
        return 0;
}

mtype gen_abs_ceil(mtype a)
{
    if(a>0) return ceil(a);
    else return floor(a);
}

