#include <stdio.h>
#include "matrix.h"

mtype mat_median(MATRIX A)
{
    MATRIX B;
    mtype med;
    int left = 0, right, pos, i, k;
    mtype t, pivot;
    B = mat_vectorize_tr(A, NULL);
    right = MatRow(B)*MatCol(B)-1;
    k = (right+1)/2;
 #define __swap(a, b) {\
		t = B[0][(a)];\
		B[0][(a)] = B[0][(b)];\
        B[0][(b)] = t;\
	}
	while(left<right)
    {
		pivot = B[0][k];
		__swap(k, right);
		for(i=pos=left; i<right; ++i)
		{
		    if(B[0][i]<pivot)
			{
			    __swap(i, pos);
				++pos;
			}
		}
		__swap(right, pos);
		if(pos==k) break;
		if(pos<k) left = pos+1;
		else right = pos-1;
	}
	med = B[0][k];
    mat_free(B);
    return med;
}

mtype mat_order_statistic(MATRIX A, int k)
{
    MATRIX B;
    mtype korder;
    int left = 0, right, pos, i;
    mtype t, pivot;
    B = mat_vectorize_tr(A, NULL);
    right = MatRow(B)*MatCol(B)-1;
 #define __swap(a, b) {\
		t = B[0][(a)];\
		B[0][(a)] = B[0][(b)];\
        B[0][(b)] = t;\
	}
	while(left<right)
    {
		pivot = B[0][k];
		__swap(k, right);
		for(i=pos=left; i<right; ++i)
		{
			if(B[0][i]<pivot)
			{
				__swap(i, pos);
				++pos;
			}
		}
		__swap(right, pos);
		if(pos==k) break;
		if(pos<k) left = pos+1;
		else right = pos-1;
	}
    korder = B[0][k];
    mat_free(B);
    return korder;
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


