#include <stdio.h>
#include <math.h>
#include "matrix.h"

MATRIX mat_evalpoly(MATRIX a, mtype x, int direction)
{
    int m, n, i, j;
    mtype r;
    MATRIX b;
    m = MatRow(a);
    n = MatCol(a);

    if (direction == 1)
    {
        b = mat_creat(1, n, UNDEFINED);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for (j=m-1; j>=0; --j) r=r*x+a[j][i];
            b[0][i] = r;
        }
    }
    else
    {
        b = mat_creat(m, 1, UNDEFINED);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for (j=n-1; j>=0; --j) r=r*x+a[i][j];
            b[i][0] = r;
        }
    }
    return (b);
}

MATRIX mat_dpoly(MATRIX a, int direction)
{
    int m, n, i, j;
    MATRIX b;
    m = MatRow(a);
    n = MatCol(a);

    if (direction == 1)
    {
        b = mat_creat(m-1, n, UNDEFINED);
        for(i=0; i<n; ++i)
        {
            for (j=m-2; j>=0; --j) b[j][i] = (j+1)*a[j+1][i];
        }
    }
    else
    {
        b = mat_creat(m, n-1, UNDEFINED);
        for(i=0; i<m; ++i)
        {
            for (j=n-2; j>=0; --j) b[i][j]= (j+1)*a[i][j+1];
        }
    }
    return (b);
}

MATRIX mat_devalpoly(MATRIX a, mtype x, int direction)
{
    int m, n, i, j;
    mtype r;
    MATRIX b;
    m = MatRow(a);
    n = MatCol(a);

    if (direction == 1)
    {
        b = mat_creat(1, n, UNDEFINED);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for (j=m-2; j>=0; --j) r=r*x+a[j+1][i]*(j+1);
            b[0][i] = r;
        }
    }
    else
    {
        b = mat_creat(m, 1, UNDEFINED);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for (j=n-2; j>=0; --j) r=r*x+a[i][j+1]*(j+1);
            b[i][0] = r;
        }
    }
    return (b);
}



