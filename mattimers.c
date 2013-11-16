#include "matrix.h"


clock_t MAT_CLOCK_TIME;

void mat_tic(void)
{
    MAT_CLOCK_TIME = clock();
}

double mat_toc(void)
{
    double diff = ((double)clock() - (double)MAT_CLOCK_TIME)*1000/(double)CLOCKS_PER_SEC;
    return diff;
}

void mat_toc_print(void)
{
    double diff = ((double)clock() - (double)MAT_CLOCK_TIME)*1000/(double)CLOCKS_PER_SEC;
    printf("Running Time since last tic: %g ms.\n", diff);
}


