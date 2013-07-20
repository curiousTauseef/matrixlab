#include <stdio.h>
#include "matrix.h"

mtype mat_get_middle(MATRIX A)
{
    int	i=0, j=0, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(m==1)
    {
        i= n/2;
        if (n%2==0)  return(A[i-1][0]+A[i][0])/2;
        else return(A[i][0]);
    }
    else  if(n==1)
    {
        j = m/2;
        if (m%2==0)  return(A[0][j-1]+A[0][j])/2;
        else return(A[0][j]);
    }
    else return -1;
}

mtype mat_median(MATRIX A)
{
    MATSTACK B;
    mtype med;
    B = mat_qsort( A, 1, NULL);
    if(MatRow(B[0])*MatCol(B[0])!=1) B[0] = mat_vectorize(B[0], B[0]);
    med =  mat_get_middle(B[0]);
    matstack_free(B);
    return med;
}

MATSTACK mat_qsort(MATRIX A, int dim, MATSTACK result)
{
    int m, n, i, j;
    MATRIX B = NULL, ind = NULL;
    m = MatCol(A);
    n = MatRow(A);
    if(result== NULL)
    {
        if ((result = matstack_creat(2)) == NULL)
            return matstack_error(MATSTACK_MALLOC);
    }

    if(dim==1)
    {
        result[0] = mat_copy(A, NULL);
        ind = mat_creat(n, m, UNDEFINED);
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j) ind[i][j] = j;
            __quicksort(result[0], 0, m-1, i, ind);
        }
    }
    else  if(dim==2)
    {
        B = mat_tran(A, NULL);
        ind = mat_creat(m, n, UNDEFINED);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<n; ++j) ind[i][j] = j;
            __quicksort(B, 0, n-1, i, ind);
        }
        result[0] = mat_tran(B, result[0]);
        mat_free(B);
    }
    else gen_error(GEN_BAD_TYPE);
    result[1] = ind;

    return result;
}

void __quicksort(MATRIX a, int l, int r, int offset, MATRIX ind)
{
    mtype t;
    int i=l, m;
    int j=r;
    mtype pivot=a[offset][l];
    mtype pivot_ind=0;
    if(ind!=NULL) pivot_ind = ind[offset][l];

    m = MatCol(a);
    if(l>=r) return;

    while ( j > i )
    {
        while ( a[offset][i] <= pivot && i<m) i++ ;
        while ( a[offset][j] > pivot && j>=0) j-- ;
        if ( j > i )
        {
            t = a[offset][i] ;
            a[offset][i] = a[offset][j] ;
            a[offset][j] = t ;
            if(ind!=NULL)
            {
                t = ind[offset][i] ;
                ind[offset][i] = ind[offset][j] ;
                ind[offset][j] = t ;
            }
        }
    }

    a[offset][l] = a[offset][j];
    a[offset][j] = pivot;
    if(ind!=NULL)
    {
        ind[offset][l] = ind[offset][j];
        ind[offset][j] = pivot_ind;
    }
    __quicksort(a,l,j-1, offset, ind);
    __quicksort(a,j+1,r, offset, ind);
}


