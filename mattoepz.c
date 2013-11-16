#include <stdlib.h>
#include "matrix.h"


MATRIX mat_symtoeplz(MATRIX R, MATRIX result)
{
    int i, j, n;
    n = MatRow(R);
    if(result==NULL) if((result = mat_creat(n, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
        for(j=0; j<n; ++j)
        {
            result[i][j] = R[abs(i-j)][0];
        }
    return result;
}

