#include <stdio.h>
#include "matrix.h"
#include <malloc.h>


INT_VECTOR int_vec_find (INT_VECTOR a, int rel_type, int n)
{
    int l, i, flag = 0, tmp0;
    INT_VECTOR indices = int_vec_creat(1, UNDEFINED);
    l = Int_VecLen(a);
    indices[0]= -1;
    for (i =0; i<l; ++i)
    {
        switch(rel_type)
        {
        case GEN_GREATER_THAN:
            tmp0 = a[i]>n;
            break;
        case GEN_LESS_THAN:
            tmp0 = a[i]<n;
            break;
        case GEN_NOT_EQUAL_TO:
            tmp0 = a[i]!=n;
            break;
        case GEN_GREATER_THAN_EQUAL_TO:
            tmp0 = a[i]>=n;
            break;
        case GEN_LESS_THAN_EQUAL_TO:
            tmp0 = a[i]<=n;
            break;
        case GEN_EQUAL_TO:
        default:
            tmp0 = a[i]==n;
            break;
        }
        if(tmp0==1)
        {
            if(flag ==0)
            {
                indices[0] = i;
                flag = 1;
            }
            else indices = int_vec_append(indices, i);
        }
    }
    return indices;
}

INT_VECSTACK mat_find(MATRIX a, int rel_type, mtype x)
{
    int m, n, i, j, flag = 0, tmp0;
    INT_VECSTACK indices = int_vecstack_creat(2);
    indices[0] = int_vec_creat(1, UNDEFINED);
    indices[1] = int_vec_creat(1, UNDEFINED);
    m = MatCol(a);
    n = MatRow(a);
    indices[0][0]= -1;
    indices[1][0]= -1;
    for (i = 0; i<n; ++i)
        for (j = 0; j<m; ++j)
        {
            switch(rel_type)
            {
            case GEN_GREATER_THAN:
                tmp0 = a[i][j]>x;
                break;
            case GEN_LESS_THAN:
                tmp0 = a[i][j]<x;
                break;
            case GEN_NOT_EQUAL_TO:
                tmp0 = a[i][j]!=x;
                break;
            case GEN_GREATER_THAN_EQUAL_TO:
                tmp0 = a[i][j]>=x;
                break;
            case GEN_LESS_THAN_EQUAL_TO:
                tmp0 = a[i][j]<=x;
                break;
            case GEN_EQUAL_TO:
            default:
                tmp0 = a[i][j]==x;
            }
            if(tmp0==1)
            {
                if(flag ==0)
                {
                    indices[0][0] = i;
                    indices[1][0] = j;
                    flag = 1;
                }
                else
                {
                    indices[0] = int_vec_append(indices[0], i);
                    indices[1] = int_vec_append(indices[1], j);
                }
            }
        }
    return indices;
}

