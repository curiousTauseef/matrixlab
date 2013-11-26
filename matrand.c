#include <stdlib.h>
#include "matrix.h"


unsigned int MAT_SEED = 0;
int MAT_SET_SEED = 0;

MATRIX mat_rand(int n, int m, MATRIX result)
{
    int i, j;
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    if(!MAT_SET_SEED)mat_set_seed(0);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = ((mtype)rand())/((mtype)RAND_MAX+1.0);
        }
    }
    return result;
}

MATRIX mat_randn(int n, int m, MATRIX result)
{
    int i, j;
    mtype tmp0;
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    if(!MAT_SET_SEED)mat_set_seed(0);

    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = (mtype)((mtype)rand()/((mtype)RAND_MAX+1.0f)+EPS);
        }
    }
    if(MatNumel(result)>0) srand((unsigned int)(result[0][0]*7923+1));
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            tmp0 =((mtype)rand())/((mtype)RAND_MAX+1.0f);
            if(result[i][j]!=0) result[i][j] = (mtype)(sqrt(-2.0f*log(result[i][j]))*cos(2.0f*3.141592f*tmp0));
        }
    }
    return result;
}

MATRIX mat_randexp(int n, int m, mtype mu, MATRIX result)
{
    int i, j;
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    if(!MAT_SET_SEED)mat_set_seed(0);

    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = (mtype)(-mu*log(1.0f-((mtype)rand())/((mtype)RAND_MAX+1.0f)));
        }
    }
    return result;
}

MATRIX mat_randfun(int n, int m, mtype (*fun)(mtype), mtype xmin, mtype xmax, MATRIX result)
{
    int i, j;
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    if(!MAT_SET_SEED)mat_set_seed(0);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = __mat_randfun(fun, xmin, xmax);
        }
    }
    return result;
}

void mat_set_seed(int seed)
{
    if(seed==0)
    {
        MAT_SEED = (unsigned int)time(NULL)*37519+rand();
        srand((unsigned int)MAT_SEED);
    }
    else
    {
        MAT_SEED = seed;
        srand((unsigned int)MAT_SEED);
    }
    MAT_SET_SEED = 1;
}

mtype __mat_randfun(mtype(*fun)(mtype), mtype xmin, mtype xmax)
{
    static mtype (*Fun)(mtype) = NULL, YMin, YMax;
    mtype X, Y;
    int iX;
    if(!MAT_SET_SEED)mat_set_seed(0);
    if(fun!=Fun)
    {
        Fun = fun;
        YMin = 0;
        YMax = Fun(xmin);
        for(iX=1; iX<RAND_MAX; ++iX)
        {
            X = xmin+(xmax-xmin)*iX/((mtype)RAND_MAX+1.0);
            Y = Fun(X);
            YMax = Y>YMax?Y:YMax;
        }
    }
    X = xmin+(xmax-xmin)*((mtype)rand())/((mtype)RAND_MAX+1.0);
    Y = YMin+(YMax-YMin)*((mtype)rand())/((mtype)RAND_MAX+1.0);;
    return Y <= fun(X) ? X : __mat_randfun(Fun, xmin, xmax);
}

mtype __mat_rand(void)
{
    if(!MAT_SET_SEED)mat_set_seed(0);
    return ((mtype)rand())/((mtype)RAND_MAX+1.0);
}

mtype __mat_randn(void)
{
    mtype tmp0, tmp1;
    if(!MAT_SET_SEED)mat_set_seed(0);
    tmp0 = ((mtype)rand())/((mtype)RAND_MAX+1.0f);
    tmp1 = ((mtype)rand())/((mtype)RAND_MAX+1.0f);
    if(tmp0!=0 ) tmp0 = (mtype)(sqrt(-2.0f*log(tmp0))*cos(2.0f*3.141592f*tmp1));
    return tmp0;
}

mtype __mat_randexp(mtype mu)
{
    if(!MAT_SET_SEED)mat_set_seed(0);
    return (mtype)(-mu*log(1.0f-((mtype)rand())/((mtype)RAND_MAX+1.0f)));
}

MATRIX mat_randperm(int m, int n, MATRIX result)
{
    int i, j;
    MATRIX tmp = NULL;
    if(result==NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i)
    {
        tmp = mat_randperm_n(n, tmp);
        for(j=0; j<n; ++j) result[i][j] = tmp[0][j];
    }
    mat_free(tmp);
    return result;
}

MATRIX mat_randperm_n(int n, MATRIX result)
{
    int i, j;
    mtype t = 0.0;
    if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    if(!MAT_SET_SEED) mat_set_seed(0);
    for(i=0; i<n; ++i)
        result[0][i] = i;
    for(i=0; i<n; ++i)
    {
        j = rand()%(n-i)+i;
        t = result[0][j];
        result[0][j] =result[0][i];
        result[0][i] = t;
    }
    return result;
}

INT_VECTOR int_vec_randperm(int n, INT_VECTOR result)
{
    int i, j;
    int t = 0;
    if(result==NULL) if((result = int_vec_creat(n, UNDEFINED))==NULL)
            return int_vec_error(INT_VEC_MALLOC);
    if(!MAT_SET_SEED)mat_set_seed(0);
    for(i=0; i<n; ++i)
        result[i] = i;
    for(i=0; i<n; ++i)
    {
        j = rand()%(n-i)+i;
        t = result[j];
        result[j] =result[i];
        result[i] = t;
    }
    return result;
}

