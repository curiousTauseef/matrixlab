#include "matrix.h"


MATRIX mat_submat(MATRIX A, int i, int j, MATRIX result)
{
    int m0, m1, p, p1, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result= mat_creat(n-1, m-1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(m0=m1=0; m<n; ++m)
    {
        if(m0==i) continue;
        for(p=p1=0; p<m; ++p)
        {
            if(p==j) continue;
            result[m1][p1] = A[m0][p];
            ++p1;
        }
        ++m1;
    }
    return result;
}

