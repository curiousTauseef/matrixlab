#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

MATRIX mat_symtoeplz(MATRIX R)
{
    int i, j, n;
    MATRIX T;

    n = MatRow(R);
    T = mat_creat(n, n, UNDEFINED);

    for(i=0; i<n; ++i)
        for(j=0; j<n; ++j)
        {
            T[i][j] = R[abs(i-j)][0];
        }
    return (T);
}

