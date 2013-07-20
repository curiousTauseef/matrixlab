#include <stdio.h>
#include "matrix.h"
#include <math.h>
#include <limits.h>

__inline mtype huber_wt(mtype x, mtype k)
{
    if(fabs(x)<= k) return 1.0;
    else return (mtype)(k/fabs(x));
}

__inline mtype bisquare_wt(mtype x, mtype k)
{
    mtype a;
    if(fabs(x)<= k)
    {
        a = x/k;
        a = 1 - (a*a);
        return (a*a);
    }
    else return 0.0;
}

mtype arcsinh(mtype x)
{
    mtype y;
    if (fabs(x) > 1.0e10)
        return (mtype)((x > 0.0) ? 0.69314718055995+log(fabs(x)) :-0.69314718055995+log(fabs(x)));
    else
    {
        y=x*x;
        return (mtype)((x == 0.0f) ? 0.0f : ((x > 0.0f) ?
                                    logplusone((float)(fabs(x)+y/(1.0f+sqrt(1.0f+y)))) :
                                    -logplusone((float)(fabs(x)+y/(1.0f+sqrt(1.0f+y))))));
    }
}

mtype arccosh(mtype x)
{
    return (mtype)((x <= 1.0) ? 0.0 : ((x > 1.0e10) ?
                                0.69314718055995+log(x) :
                                log(x+sqrt((x-1.0)*(x+1.0)))));
}

mtype arctanh(mtype x)
{
    mtype ax;
    if (fabs(x) >= 1.0)return (mtype)((x > 0.0) ? DOUBLE_MAX : -DOUBLE_MAX);
    else
    {
        ax=(mtype)fabs(x);
        return (mtype)((x == 0.0) ? 0.0 : ((x > 0.0) ? 0.5*logplusone(2.0f*ax/(1.0f-ax)) :
                                    -0.5*logplusone(2.0f*ax/(1.0f-ax))));
    }
}

__inline mtype logplusone(mtype x)
{
    mtype y,z;
    if (x == 0.0) return 0.0;
    else if (x < -0.2928 || x > 0.4142) return (mtype)log(1.0f+x);
    else
    {
        z=x/(x+2.0f);
        y=z*z;
        return (mtype)(z*(2.0+y*(0.66666666663366+y*(0.400000001206045+y*
                                             (0.285714091590488+y*(0.22223823332791+y*(0.1811136267967+y*0.16948212488)))))));
    }
}

MATRIX mat_huber_wt( MATRIX A, mtype k, mtype sigma)
{
    int	i, j, m, n;
    MATRIX	B;
    m = MatCol(A);
    n = MatRow(A);

    if ((B = mat_creat( n, m, UNDEFINED )) == NULL)
        return mat_error(MAT_MALLOC);;

    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            B[i][j] = huber_wt(A[i][j]/sigma, k);
        }
    return (B);
}

MATRIX mat_bisquare_wt( MATRIX A, mtype k, mtype sigma)
{
    int	i, j, m, n;
    MATRIX	B;
    m = MatCol(A);
    n = MatRow(A);

    if ((B = mat_creat( n, m, UNDEFINED )) == NULL)
        return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            B[i][j] = bisquare_wt(A[i][j]/sigma, k);
        }
    return (B);
}

MATRIX mat_gfunc( MATRIX A, mtype (*pt2func)(mtype), MATRIX result)
{
    int	i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result== NULL|| result ==A)if ((result = mat_creat( n, m, UNDEFINED )) == NULL)
            return mat_error(MAT_MALLOC);


    #pragma omp parallel for private(j)
    for (i=0; i<n; ++i)
        for (j=0; j<m; ++j)
        {
            result[i][j] = (*pt2func)(A[i][j]);
        }
    return (result);
}





