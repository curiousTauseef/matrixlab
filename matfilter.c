#include <stdio.h>
#include <math.h>
#include "matrix.h"

MATRIX mat_conv2(MATRIX a, MATRIX mask, MATRIX scratch, MATRIX result)
{
    int i, j, k, l, m, n, o, p, ii, jj, mm, nn, flag = 0;
    mtype acc = 0.0;
    MATRIX b = scratch, c = result;
    m = MatCol(a);
    n = MatRow(a);
    o = MatCol(mask);
    p = MatRow(mask);
    if((o%2)!=1 ||(p%2)!=1) gen_error(GEN_SIZE_ERROR);
    ii = (o - 1)/2;
    jj = (p - 1)/2;
    mm = ii+ii+m;
    nn = jj+jj+n;
    k = ii+m;
    l = jj+n;
    if(scratch == NULL)
    {
        b = mat_creat(nn, mm, UNDEFINED);
        flag = 1;
        for(i=0; i<mm; ++i)
        {
            for(j=0; j<nn; ++j)
            {
                if(i<ii || j<jj || i>=k || j>=l ) b[j][i] = 0.0;
                else b[j][i] = a[j-jj][i-ii];
            }
        }
    }

    if(result == NULL)
    {
        c = mat_creat(n, m, UNDEFINED);
    }
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n; ++j)
        {
            acc = 0.0;
            for(k = -ii; k<=ii; ++k)
                for(l = -jj; l<=jj; ++l)
                    acc += b[j+jj+l][i+ii+k]*mask[jj-l][ii-k];
            c[j][i] = acc;
        }
    }
    if(flag == 1)mat_free(b);
    return c;
}

