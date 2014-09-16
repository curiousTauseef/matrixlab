#include <string.h>
#include <stdlib.h>
#include "matrix.h"


/** \brief Checks if scalar is NaN
 *
 * \param[in] x Input scalar
 * \return Zero/non-zero
 *
 */

int mats_isnan(mtype x)
{
    volatile mtype temp = x;
    return temp != x;
}

/** \brief Checks if scalar is infinite
 *
 * \param[in] x Input scalar
 * \return Zero/non-zero
 *
 */

int mats_isinf(mtype x)
{
    volatile mtype temp = x;
    if ((temp == x) && ((temp - x) != 0.0))
        return (x < 0.0 ? -1 : 1);
    else return 0;
}

/** \brief Prints nextline to stdout
 *
 *
 */

void mat_nextline(void)
{
    mat_fnextline(stdout);
}

/** \brief Prints nextline to file
 *
 * \param[in] fp Pointer to opened file
 *
 */

void mat_fnextline(MAT_FILEPOINTER fp)
{
    fprintf(fp, "\n");
}

/** \brief Computes element-wise binary function for two matrices
 *
 * \param[in] A First matrix
 * \param[in] B Second matrix
 * \param[in] func Pointer to the function
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_bsxfun(MATRIX A, MATRIX B, mtype (*func)(mtype, mtype), MATRIX result)
{
    int m, n, o, p, i, j;
    m = MatRow(A);
    n = MatCol(A);

    o = MatRow(B);
    p = MatCol(B);
    if(m<o && n==p && m==1)
    {
        if(result== NULL) if((result = mat_creat(o, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<o; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*func)(A[0][j], B[i][j]);
            }
        }
    }
    else if(m>o && n==p && o==1)
    {
        if(result== NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*func)(A[i][j], B[0][j]);
            }
        }
    }
    else if(m==o && n<p && n==1)
    {
        if(result== NULL) if((result = mat_creat(m, p, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<p; ++j)
            {
                result[i][j] = (*func)(A[i][0], B[i][j]);
            }
        }
    }
    else if(m==o && n>p && p==1)
    {
        if(result== NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*func)(A[i][j], B[i][0]);
            }
        }
    }
    else if(m==1 && p==1)
    {
        if(result== NULL) if((result = mat_creat(o, n, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<o; ++i)
        {
            for(j=0; j<n; ++j)
            {
                result[i][j] = (*func)(A[0][j], B[i][0]);
            }
        }
    }
    else if(n==1 && o==1)
    {
        if(result== NULL) if((result = mat_creat(m, p, UNDEFINED))==NULL) return mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=0; j<p; ++j)
            {
                result[i][j] = (*func)(A[i][0], B[0][j]);
            }
        }
    }
    else mat_error(MAT_SIZEMISMATCH);
    return result;
}

/** \brief Computes a randomly permutation of first k positive integers
 *
 * \param[in] n Number of random permutations to make
 * \param[in] k Integer upto which it will consider
 * \param[in] result Vector to store the result
 * \return Permuted vector
 *
 */

INT_VECTOR int_vec_permute_vect(int n, int k, INT_VECTOR result)
{
    int i, j, d;
    if(result== NULL) if((result = int_vec_creat(k, SERIES_INT_VECTOR))==NULL)
            return int_vec_error(INT_VEC_MALLOC);
    for(i=0; i<n-1; i++)
    {
        j = result[i];
        d = i + (rand()%(n-i));
        result[i] = result[d];
        result[d] = j;
    }
    return result;
}

/** \brief Extracts sub-vector from an integer vector
 *
 * \param[in] data Input vector
 * \param[in] indices Indices to extracted
 * \return Extracted vector
 *
 */

INT_VECTOR mat_get_sub_vector(INT_VECTOR data, INT_VECTOR indices)
{
    int i, n;
    INT_VECTOR subvec;
    n = Int_VecLen(indices);
    subvec = int_vec_creat(n, UNDEFINED);
    for(i=0; i<n; ++i)subvec[i] = data[indices[i]];
    return subvec;
}

/** \brief Extracts sub-matrix from rows of a matrix
 *
 * \param[in] data Input matrix
 * \param[in] indices Rows to extract
 * \param[in] result Matrix to store the result
 * \return Extracted matrix
 *
 */

MATRIX mat_get_sub_matrix_from_rows(MATRIX data, INT_VECTOR indices, MATRIX result)
{
    int i, j, k, n;
    k = MatCol(data);
    n = Int_VecLen(indices);
    if(result==NULL) if((result = mat_creat(n, k, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<k; ++j) result[i][j] = data[indices[i]][j];
    }
    return result;
}

/** \brief Extracts sub-matrix from columns of a matrix
 *
 * \param[in] data Input matrix
 * \param[in] indices Columns to extract
 * \param[in] result Matrix to store the result
 * \return Extracted matrix
 *
 */

MATRIX mat_get_sub_matrix_from_cols(MATRIX data, INT_VECTOR indices, MATRIX result)
{
    int i, j, k, n;
    k = MatRow(data);
    n = Int_VecLen(indices);
    if(result==NULL) if((result = mat_creat(k, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i)
    {
        for(j=0; j<k; ++j) result[j][i] = data[j][indices[i]];
    }
    return result;
}

/** \brief Computes the Euclidean distances of  points from a given point
 *
 * \param[in] data Points matrix (d x N)
 * \param[in] curr_data Matrix point from which the distance to be computed (d x 1)
 * \param[in] result Matrix to store the result
 * \return Euclidean distance matrix
 *
 */

MATRIX mat_calc_dist_sq(MATRIX data, MATRIX curr_data, MATRIX result)
{
    int i, j, m, n;
    mtype dist;
    m = MatRow(data);
    n = MatCol(data);
    if(result==NULL) if((result = mat_creat(1, n, ZERO_MATRIX))==NULL) return mat_error(MAT_MALLOC);;
    for(i=0; i<n; ++i)
    {
        dist = 0.0;
        for(j=0; j<m; ++j)
        {
            dist += (data[j][i]-curr_data[j][0])*(data[j][i]-curr_data[j][0]);
        }
        result[0][i] = dist;
    }
    return result;
}

/** \brief Finds points within a neighborhood
 *
 * \param[in] data Points matrix (d x N)
 * \param[in] curr_data Matrix point from which the distance to be computed (d x 1)
 * \param[in] range Radius to search within
 * \return Indices Vector
 *
 */

INT_VECTOR mat_find_within_dist(MATRIX data, MATRIX curr_data, mtype range)
{
    int i, count = 0, n;
    INT_VECTOR indices = 0;
    MATRIX dist = NULL;
    mtype tmp0;
    n = MatCol(data);
    tmp0 = range*range;
    dist = mat_calc_dist_sq(data, curr_data, dist);
    for(i=0; i<n; ++i)
    {
        if(dist[0][i]<tmp0)
        {
            ++count;
            indices = (int*) realloc(indices, (count+1)*sizeof(int));
            indices[count] = i;
        }
    }
    indices[0] = count;
    indices = (int *)indices+1;
    mat_free(dist);
    return indices;
}

/** \brief Picks a row from a matrix
 *
 * \param[in] data Input matrix
 * \param[in] r Row index
 * \param[in] result Matrix to store the result
 * \return Row matrix
 *
 */

MATRIX mat_pick_row(MATRIX data, int r, MATRIX result)
{
    int m, i;
    m = MatCol(data);
    if(result==NULL) if((result = mat_creat(1, m, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i) result[0][i] = data[r][i];
    return result;
}

/** \brief Picks a column from a matrix
 *
 * \param[in] data Input matrix
 * \param[in] r Column index
 * \param[in] result Matrix to store the result
 * \return Column matrix
 *
 */

MATRIX mat_pick_col(MATRIX data, int c, MATRIX result)
{
    int n, i;
    n = MatRow(data);
    if(result==NULL) if((result = mat_creat(n, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(i=0; i<n; ++i) result[i][0] = data[i][c];
    return result;
}

void __mat_cart2pol(mtype x, mtype y, mtype *rho, mtype *th)
{
    *rho = (mtype)sqrt(x*x + y*y);
    *th = (mtype)atan2(y, x);
}

/**< \cond HIDDEN_SYMBOLS */
#define __cart2pol(x, y, rho, th)\
{\
    rho = (mtype)sqrt(x*x + y*y);\
    th = (mtype)atan2(y, x);\
}
/**< \endcond */

void __mat_pol2cart(mtype rho, mtype th, mtype *x, mtype *y)
{
    *x = (mtype)(rho*cos(th));
    *y = (mtype)(rho*sin(th));
}

/**< \cond HIDDEN_SYMBOLS */
#define __pol2cart(rho, th, x, y)\
{\
    x = (mtype)(rho*cos(th));\
    y = (mtype)(rho*sin(th));\
}
/**< \endcond */

/** \brief Converts Cartesian co-ordinates to polar co-ordinates
 *
 * \param[in] A Input matrix
 * \param[in] dim Data order ROWS/COLS
 * \return Polar co-ordinate matrix
 *
 */

MATRIX mat_cart2pol(MATRIX A, int dim, MATRIX result)
{
    int m, n, i;
    m = MatCol(A);
    n = MatRow(A);
    if(dim==0 && n>1)
    {
        if(result==NULL)if((result = mat_creat(2, m, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<m; ++i)
            __cart2pol(A[0][i], A[1][i], result[0][i], result[1][i]);
        return result;
    }
    if(dim==1 && m>1)
    {
        if(result==NULL)if((result = mat_creat(n, 2, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<n; ++i)
            __cart2pol(A[i][0], A[i][1], result[i][0], result[i][1]);
        return result;
    }
    return mat_error (MAT_SIZEMISMATCH);
}

/** \brief Converts polar co-ordinates to Cartesian co-ordinates
 *
 * \param[in] A Input matrix
 * \param[in] dim Data order ROWS/COLS
 * \return Cartesian co-ordinate matrix
 *
 */

MATRIX mat_pol2cart(MATRIX A, int dim, MATRIX result)
{
    int m, n, i;
    m = MatCol(A);
    n = MatRow(A);
    if(dim==0 && n>1)
    {
        if(result==NULL)if((result = mat_creat(2, m, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<m; ++i)
            __pol2cart(A[0][i], A[1][i], result[0][i], result[1][i]);
        return result;
    }
    if(dim==1 && m>1)
    {
        if(result==NULL)if((result = mat_creat(n, 2, ZERO_MATRIX))==NULL) return NULL;
        for(i=0; i<n; ++i)
            __pol2cart(A[i][0], A[i][1], result[i][0], result[i][1]);
        return result;
    }
    return mat_error (MAT_SIZEMISMATCH);
}

/**< \cond HIDDEN_SYMBOLS */
int __mat_powerof2(int width, int *m, int *twopm)
{
    *m = (int)ceil(log((mtype)width)/log((mtype)2.0));
    *twopm = (int)pow((mtype)2.0, (mtype)*m);

    return 1;
}
/**< \endcond */

/** \brief Extract only the unique integers from an integer vector
 *
 * \param[in] a Input vector
 * \return Unique vector
 *
 */

INT_VECTOR int_vec_unique(INT_VECTOR a)
{
    int i, l, uni_l = 0;
    mtype *ordered = NULL, **p_ordered;
    INT_VECTOR int_ordered = NULL;
    MAT_TREE s;
    p_ordered = &ordered;
    l = Int_VecLen(a);
    s = mat_bs_make_null();
    for (i = 0; i<l; ++i)
    {
        s = mat_bs_insert((mtype)a[i],s);
    }
    uni_l = mat_bs_inorder(s, uni_l, p_ordered);
    int_ordered = int_vec_creat(uni_l, UNDEFINED);
    for(i = 0; i< uni_l; ++i) int_ordered[i] = (int)ordered[i];
    free(ordered);
    s = mat_bs_free(s);
    return int_ordered;
}

