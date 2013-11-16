#include "matrix.h"


MATRIX mat_evalpoly(MATRIX a, mtype x, int dir, MATRIX result)
{
    int m, n, i, j;
    mtype r;
    m = MatRow(a);
    n = MatCol(a);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for(j=m-1; j>=0; --j) r=r*x+a[j][i];
            result[0][i] = r;
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for(j=n-1; j>=0; --j) r=r*x+a[i][j];
            result[i][0] = r;
        }
    }
    return result;
}

MATRIX mat_dpoly(MATRIX a, int dir, MATRIX result)
{
    int m, n, i, j;
    m = MatRow(a);
    n = MatCol(a);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(m-1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            for(j=m-2; j>=0; --j) result[j][i] = (j+1)*a[j+1][i];
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, n-1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=n-2; j>=0; --j) result[i][j]= (j+1)*a[i][j+1];
        }
    }
    return result;
}

MATRIX mat_devalpoly(MATRIX a, mtype x, int dir, MATRIX result)
{
    int m, n, i, j;
    mtype r;
    m = MatRow(a);
    n = MatCol(a);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for(j=m-2; j>=0; --j) r=r*x+a[j+1][i]*(j+1);
            result[0][i] = r;
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for(j=n-2; j>=0; --j) r=r*x+a[i][j+1]*(j+1);
            result[i][0] = r;
        }
    }
    return result;
}

