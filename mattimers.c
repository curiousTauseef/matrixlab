#include "matrix.h"

#if defined(__STDC_VERSION__) && __STDC_VERSION__>201112L
_Thread_local clock_t MAT_CLOCK_TIME;
#elif defined(__GNUC__) || defined (__clang__)
__thread clock_t MAT_CLOCK_TIME;
#elif defined(_MSC_VER) || defined (__INTEL_COMPILER)
__declspec(thread) clock_t MAT_CLOCK_TIME;
#else
clock_t MAT_CLOCK_TIME;
#endif

/** \brief Starts stopwatch timer
 *
 *
 */

void mat_tic(void)
{
    MAT_CLOCK_TIME = clock();
}

/** \brief Computes elapsed time from last start of timer
 *
 * \return Elapsed time
 *
 */

double mat_toc(void)
{
    double diff = ((double)clock() - (double)MAT_CLOCK_TIME)*1000/(double)CLOCKS_PER_SEC;
    return diff;
}

/** \brief Computes and prints elapsed time from last start of timer on the stdout
 *
 *
 */

void mat_toc_print(void)
{
    double diff = ((double)clock() - (double)MAT_CLOCK_TIME)*1000/(double)CLOCKS_PER_SEC;
    printf("Running Time since last tic: %g ms.\n", diff);
}


