#include <stdio.h>
#include "matrix.h"


MATRIX mat_submat(MATRIX A, int i, int j )
{
    int	m0, m1, p, p1, m, n;
    MATRIX	S;
    m = MatCol(A);
    n = MatRow(A);
    S = mat_creat(n-1, m-1, UNDEFINED);
    for (m0=m1=0; m<n; ++m)
    {
        if (m0==i) continue;
        for (p=p1=0; p<m; ++p)
        {
            if (p==j) continue;
            S[m1][p1] = A[m0][p];
            ++p1;
        }
        ++m1;
    }
    return (S);
}

