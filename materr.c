#include <stdlib.h>
#include "matrix.h"


/** \brief Generates error message for general errors and exits
 *
 * \param[in] err Error type (GEN_NOT_CONVERGED/GEN_FNOTOPEN/
        GEN_FNOTGETMAT/GEN_SIZEMISMATCH/GEN_MATH_ERROR/GEN_MALLOC/GEN_NOT_FOUND/GEN_SIZE_ERROR/GEN_BAD_TYPE)
 *
 */

int gen_error(int err_)
{
    switch(err_)
    {
    case GEN_NOT_CONVERGED:
        fprintf(stderr, "gen: convergence  failure\n");
        break;
    case GEN_FNOTOPEN:
        fprintf(stderr, "mat: fileopen error\n");
        break;
    case GEN_FNOTGETMAT:
        fprintf(stderr, "fgetmat: matrix read error\n");
        break;
    case GEN_SIZEMISMATCH:
        fprintf(stderr, "gen: dimension size mismatch\n");
        break;
    case GEN_MATH_ERROR:
        fprintf(stderr, "gen: math error\n");
        break;
    case GEN_MALLOC:
        fprintf(stderr, "gen: malloc error\n");
        break;
    case GEN_NOT_FOUND:
        fprintf(stderr, "gen: not found\n");
        break;
    case GEN_SIZE_ERROR:
        fprintf(stderr, "gen: dimension size error\n");
        break;
    case GEN_BAD_TYPE:
        fprintf(stderr, "gen: bad type error\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for matrix errors and exits
 *
 * \param[in] err Error type (MAT_MALLOC/MAT_FNOTOPEN/MAT_FNOTGETMAT/MAT_SIZEMISMATCH/
            MAT_INVERSE_ILL_COND/MAT_INVERSE_NOT_SQUARE/MAT_CHOLESKY_FAILED)
 *
 */

MATRIX mat_error(int err_)
{
    switch(err_)
    {
    case MAT_MALLOC:
        fprintf(stderr, "mat: malloc error\n");
        break;
    case MAT_FNOTOPEN:
        fprintf(stderr, "mat: fileopen error\n");
        break;
    case MAT_FNOTGETMAT:
        fprintf(stderr, "fgetmat: matrix read error\n");
        break;
    case MAT_SIZEMISMATCH:
        fprintf(stderr, "mat: matrix size mismatch\n");
        break;
    case MAT_INVERSE_ILL_COND:
        fprintf(stderr, "mat: matrix inverse ill conditioned\n");
        break;
    case MAT_INVERSE_NOT_SQUARE:
        fprintf(stderr, "mat: matrix not square for inverse\n");
        break;
    case MAT_CHOLESKY_FAILED:
        fprintf(stderr, "mat: matrix cholesky failed due to failure in positive-definiteness\n");
    }
    exit(1);
}

/** \brief Generates error message for matrix stack errors and exits
 *
 * \param[in] err Error type (MATSTACK_MALLOC/MATSTACK_FNOTOPEN/MATSTACK_FNOTGETMAT/MATSTACK_SIZEMISMATCH/
            MATSTACK_INVERSE_ERROR)
 *
 */

MATSTACK matstack_error(int err_)
{
    switch(err_)
    {
    case MATSTACK_MALLOC:
        fprintf(stderr, "matstack: malloc error\n");
        break;
    case MATSTACK_FNOTOPEN:
        fprintf(stderr, "matstack: fileopen error\n");
        break;
    case MATSTACK_FNOTGETMAT:
        fprintf(stderr, "fgetmatstack: matrixstack read error\n");
        break;
    case MATSTACK_SIZEMISMATCH:
        fprintf(stderr, "matstack: matrixstack size mismatch\n");
        break;
    case MATSTACK_INVERSE_ERROR:
        fprintf(stderr, "matstack: matrixstack ill conditioned\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for integer vector errors and exits
 *
 * \param[in] err Error type (INT_VEC_MALLOC/INT_VEC_FNOTOPEN/INT_VEC_FNOTGETINT_VEC/INT_VEC_SIZEMISMATCH)
 *
 */

INT_VECTOR int_vec_error(int err_)
{
    switch(err_)
    {
    case INT_VEC_MALLOC:
        fprintf(stderr, "int_vec: malloc error\n");
        break;
    case  INT_VEC_FNOTOPEN:
        fprintf(stderr, "int_vec: fileopen error\n");
        break;
    case INT_VEC_FNOTGETINT_VEC:
        fprintf(stderr, "fgetint_vec: int_vector read error\n");
        break;
    case  INT_VEC_SIZEMISMATCH:
        fprintf(stderr, "int_vec: int_vector size mismatch\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for integer vector stack errors and exits
 *
 * \param[in] err Error type (INT_VECSTACK_MALLOC/INT_VECSTACK_FNOTOPEN/INT_VECSTACK_FNOTGETINT_VEC/INT_VECSTACK_SIZEMISMATCH)
 *
 */

INT_VECSTACK int_vecstack_error(int err_)
{
    switch(err_)
    {
    case INT_VECSTACK_MALLOC:
        fprintf(stderr, "int_vecstack: malloc error\n");
        break;
    case INT_VECSTACK_FNOTOPEN:
        fprintf(stderr, "int_vecstack: fileopen error\n");
        break;
    case INT_VECSTACK_FNOTGETINT_VEC:
        fprintf(stderr, "fgetint_vecstack: int_vectorstack read error\n");
        break;
    case INT_VECSTACK_SIZEMISMATCH:
        fprintf(stderr, "int_vecstack: int_vectorstack size mismatch\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for stack errors and exits
 *
 * \param[in] err Error type (STACK_MALLOC/STACK_EMPTY)
 *
 */

int stack_error(int err_)
{
    switch(err_)
    {
    case STACK_MALLOC:
        fprintf(stderr, "stack: malloc error\n");
        break;
    case  STACK_EMPTY:
        fprintf(stderr, "stack: stack empty\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for queue errors and exits
 *
 * \param[in] err Error type (QUEUE_MALLOC/QUEUE_EMPTY)
 *
 */

int queue_error(int err_)
{
    switch(err_)
    {
    case QUEUE_MALLOC:
        fprintf(stderr, "queue: malloc error\n");
        break;
    case  QUEUE_EMPTY:
        fprintf(stderr, "queue: queue empty\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for priority queue errors and exits
 *
 * \param[in] err Error type (PQ_MALLOC/PQ_EMPTY)
 *
 */

int pq_error(int err_)
{
    switch( err_)
    {
    case PQ_MALLOC:
        fprintf(stderr, "priority queue: malloc error\n");
        break;
    case  PQ_EMPTY:
        fprintf(stderr, "priority queue: queue empty\n");
        break;
    }
    exit(1);
}

/** \brief Generates error message for graph errors and exits
 *
 * \param[in] err Error type (GRAPH_MALLOC/GRAPH_READ/GRAPH_ELSE)
 *
 */

int graph_error(int err_)
{
    switch( err_)
    {
    case GRAPH_MALLOC:
        fprintf(stderr, "graph: malloc error\n");
        break;
    case GRAPH_READ:
        fprintf(stderr, "graph: read error\n");
        break;
    case GRAPH_ELSE:
        fprintf(stderr, "graph: else error\n");
        break;
    }
    exit(1);
}

