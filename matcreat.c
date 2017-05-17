#include "matrix.h"
#include <stdlib.h>


/** \cond HIDDEN_SYMBOLS */

MATRIX __mat_creat(int row, int col)
{
    MATRIX mat = NULL;
    int i;
    if((mat = (MATRIX)malloc(sizeof(mtype *)*(row+2)))==NULL)
        return (mat_error( MAT_MALLOC ));
    for(i=2; i<row+2; ++i)
    {
        if((*((mtype **)(mat + i)) = (mtype *)malloc(sizeof(mtype)*col))==NULL)
            return (mat_error(MAT_MALLOC));
    }
    *((int*)(mat)) = row;
    *((int*)(mat+1)) = col;
    return (mat+2);
}

/** \endcond */

/** \brief Creates a matrix
 *
 * \param[in] row Number of rows
 * \param[in] col Number of columns
 * \param[in] type Definition type (UNDEFINED/ZERO_MATRIX/UNIT_MATRIX/ONES_MATRIX)
 * \return Output matrix
 *
 */

MATRIX mat_creat(int row, int col, int type)
{
    MATRIX A;
    if((A =__mat_creat(row, col))!=NULL)
    {
        return (mat_fill_type(A, type));
    }
    else
        return NULL;
}

/** \cond HIDDEN_SYMBOLS */

MATSTACK __matstack_creat(int len)
{
    MATSTACK matrixstack;
    int i;
    if((matrixstack = (MATRIX *)malloc(sizeof(MATRIX)*(1+len)))==NULL)
        return (matstack_error( MATSTACK_MALLOC));

    for(i=0; i<=len; ++i)
    {
        matrixstack[i]= NULL;
    }
    matrixstack[0] = mat_creat(1,1,0);
    matrixstack[0][0][0] = len;
    return (matrixstack+1);
}

/** \endcond */

/** \brief Creates a matrix stack
 *
 * \param[in] len Length of the stack
 * \return Output matrix stack
 *
 */

MATSTACK matstack_creat(int len)
{
    MATSTACK A;
    if((A =__matstack_creat(len))!=NULL)
    {
        return A;
    }
    else
        return NULL;
}

/** \brief Appends a matrix to a matrix stack
 *
 * \param[in] s Input matrix stack
 * \param[in] A Input matrix to append
 * \return Output matrix stack
 *
 */

MATSTACK matstack_append(MATSTACK s, MATRIX A)
{
    MATSTACK m;
    int i, n;
    if(A==NULL) return s;
    n = MatStacklength(s);
    m = matstack_creat(n+1);
    for(i=0; i<n; ++i) m[i] = s[i];
    m[n] = A;
    return m;
}

/** \brief Frees a matrix stack
 *
 * \param[in] A Input matrix stack
 * \return Success
 *
 */

int matstack_free(MATSTACK A)
{
    int i, n;
    if(A==NULL) return 0;
    n = MatStacklength(A);
    A = A-1;
    for(i=0; i<n; ++i) mat_free(A[i]);
    A = NULL;
    return 1;
}

/** \brief Fills a matrix with a value
 *
 * \param[in] A Input matrix
 * \param[in] val Value to fill with
 * \return Filled matrix
 *
 */

MATRIX mat_fill(MATRIX A, mtype val)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
            A[i][j] = val;
    return A;
}

/** \brief Fills a matrix to a type
 *
 * \param[in] A Input matrix
 * \param[in] type Fill type (UNDEFINED/ZERO_MATRIX/UNIT_MATRIX/ONES_MATRIX)
 * \return Filled matrix
 *
 */

MATRIX mat_fill_type(MATRIX A, int type)
{
    int i, j, k, m, n;
    m = MatCol(A);
    n = MatRow(A);
    switch(type)
    {
    case UNDEFINED:
        break;
    case ZERO_MATRIX:
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j)
            {
                A[i][j] = 0.0;
            }
        }
        break;
    case UNIT_MATRIX:
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j)
            {
                A[i][j] = 0.0;
            }
        }
        k = (m<n)? m:n;
        for(i=0; i<k; ++i)
        {
            A[i][i] = 1.0;
        }
        break;
    case ONES_MATRIX:
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
        {
            for(j=0; j<m; ++j) A[i][j] = 1.0;
        }
        break;
    }
    return A;
}

/** \brief Frees a matrix
 *
 * \param[in] A Input matrix
 * \return Success
 *
 */

int mat_free(MATRIX A)
{
    int i, n;
    if(A==NULL)
        return (0);
    n = MatRow(A);
    for(i=0; i<n; ++i)
    {
        free(A[i]);
    }
    A = A-2;
    free(A);
    return 1;
}

/** \cond HIDDEN_SYMBOLS */

INT_VECTOR __int_vec_creat(int len)
{
    INT_VECTOR int_vector;
    if((int_vector = (int*)malloc(sizeof(int)*(len+1)))==NULL)
        return (int_vec_error(INT_VEC_MALLOC));
    int_vector[0] = len;
    return (int_vector+1);
}

/** \endcond */

/** \brief Creates an integer vector
 *
 * \param[in] len Length of the vector
 * \param[in] type Definition type (UNDEFINED/ZERO_INT_VECTOR/ONES_INT_VECTOR/SERIES_INT_VECTOR)
 * \return Output vector
 *
 */

INT_VECTOR int_vec_creat(int len, int type)
{
    INT_VECTOR A;
    if((A =__int_vec_creat(len))!=NULL)
    {
        return (int_vec_fill_type(A, type));
    }
    else
        return NULL;
}

/** \brief Fills an integer vector with a value
 *
 * \param[in] A Input vector
 * \param[in] val Value to fill with
 * \return Filled vector
 *
 */

INT_VECTOR int_vec_fill(INT_VECTOR A, int val)
{
    int i, n;
    n = Int_VecLen(A);
    for(i=0; i<n; ++i) A[i] = val;
    return A;
}

/** \brief Fills an integer vector to a type
 *
 * \param[in] A Input vector
 * \param[in] type Definition type (UNDEFINED/ZERO_INT_VECTOR/ONES_INT_VECTOR/SERIES_INT_VECTOR)
 * \return Filled vector
 *
 */

INT_VECTOR int_vec_fill_type(INT_VECTOR A, int type)
{
    int i, n;
    n = Int_VecLen(A);
    switch(type)
    {
    case UNDEFINED:
        break;
    case ZERO_INT_VECTOR:
        for(i=0; i<n; ++i) A[i] = 0;
        break;
    case ONES_INT_VECTOR:
        for(i=0; i<n; ++i) A[i] = 1;
        break;
    case  SERIES_INT_VECTOR:
        for(i=0; i<n; ++i) A[i] = i;
        break;
    }
    return A;
}

/** \brief Frees an integer vector
 *
 * \param[in] A Input vector
 * \return Success
 *
 */

int int_vec_free(INT_VECTOR A)
{
    if(A==NULL) return 0;
    Int_VecLen(A) = 0;
    free(A-1);
    A = NULL;
    return 1;
}

/** \cond HIDDEN_SYMBOLS */

INT_VECSTACK __int_vecstack_creat(int len)
{
    INT_VECSTACK int_vecstack;
    int i;
    if((int_vecstack = (INT_VECTOR *)malloc( sizeof(INT_VECTOR)*(1+len)))==NULL)
        return (int_vecstack_error(INT_VECSTACK_MALLOC));

    for(i=0; i<=len; ++i)
    {
        int_vecstack[i]= NULL;
    }
    int_vecstack[0] = int_vec_creat(1,0);
    int_vecstack[0][0] = len;
    return (int_vecstack+1);
}

/** \endcond */

/** \brief Creates an integer vector stack
 *
 * \param[in] len Length of the stack
 * \return Output vector stack
 *
 */

INT_VECSTACK int_vecstack_creat(int len)
{
    INT_VECSTACK A;
    if((A =__int_vecstack_creat(len))!=NULL)
    {
        return A;
    }
    else
        return NULL;
}

/** \brief Frees an integer vector stack
 *
 * \param[in] A Input vector stack
 * \return Success
 *
 */

int int_vecstack_free(INT_VECSTACK A)
{
    int i, n;
    if(A==NULL)
        return 0;
    n = Int_VecStackLength(A);
    A = A-1;
    for(i=0; i<=n; ++i)
        int_vec_free(A[i]);
    A = NULL;
    return 1;
}

/** \brief Creates a Bayes model
 *
 * \return Output Bayes model
 *
 */

MAT_BAYES_MODEL mat_bayes_model_creat(void)
{
    MAT_BAYES_MODEL a;
    a = (MAT_BAYES_MODEL)malloc(sizeof(mat_bayes_model));
    a->class_priors = NULL;
    a->class_covars = NULL;
    a->class_labels = NULL;
    a->class_means = 0;
    a->num_of_classes = 0;
    a->num_of_features = 0;
    return a;
}

/** \brief Frees a Bayes model
 *
 * \param[in] a Input Bayes model
 * \return Success
 *
 */

int mat_bayes_model_free(MAT_BAYES_MODEL a)
{
    matstack_free(a->class_means);
    matstack_free(a->class_covars);
    int_vec_free(a->class_labels);
    free(a);
    a = NULL;
    return 1;
}

/** \brief Creates a perceptron
 *
 * \return Output perceptron
 *
 */

MAT_PERCEPTRON mat_perceptron_creat(void)
{
    MAT_PERCEPTRON a;
    a = (MAT_PERCEPTRON)malloc(sizeof(mat_perceptron));
    a->class_weights = NULL;
    a->class_labels = NULL;
    a->istrained = 0;
    a->num_of_classes = 0;
    a->num_of_features = 0;
    a->num_of_iterations = 20;
    return a;
}

/** \brief Frees a perceptron
 *
 * \param[in] a Input perceptron
 * \return Success
 *
 */

int mat_perceptron_free(MAT_PERCEPTRON a)
{
    if(a->istrained!=0) mat_free(a->class_weights);
    int_vec_free(a->class_labels);
    free(a);
    a = NULL;
    return 1;
}

/** \brief Creates a matrix-vector pair
 *
 * \return Output matrix-vector pair
 *
 */

MATVEC_DPOINTER matvec_creat(void)
{
    MATVEC_DPOINTER a;
    a = (MATVEC_DPOINTER)malloc(sizeof(void *)*2);
    a[0] = NULL;
    a[1] = NULL;
    return a;
}

/** \brief Frees a matrix-vector pair
 *
 * \param[in] a Input matrix-vector pair
 * \return Success
 *
 *
 */

int matvec_free(MATVEC_DPOINTER a)
{
    mat_free((MATRIX)a[0]);
    int_vec_free((INT_VECTOR)a[1]);
    free(a);/* to be edited later */
    return 1;
}

/** \brief Appends an integer to an integer vector
 *
 * \param[in] a Input vector
 * \param[in] i Integer to append
 * \return Appended vector
 *
 */

INT_VECTOR int_vec_append(INT_VECTOR a, int i)
{
    int l;
    l = Int_VecLen(a);
    a--;
    if((a = (int *)realloc(a, sizeof(int)*(l+2)))==NULL) return int_vec_error(INT_VEC_MALLOC);
    a[0]= l+1;
    ++a;
    a[l] = i;
    return a;
}

/** \brief Copies an integer vector
 *
 * \param[in] a Input vector
 * \param[in] result Vector to store the result
 * \return Output vector
 *
 */

INT_VECTOR int_vec_copy(INT_VECTOR a, INT_VECTOR result)
{
    int i, m;
    m = Int_VecLen(a);
    if(result==NULL) if((result = int_vec_creat(m, UNDEFINED))==NULL) int_vec_error(INT_VEC_MALLOC);
    for(i=0; i<m; ++i) result[i] = a[i];
    return result;
}

/** \brief Copies a matrix
 *
 * \param[in] A Input matrix
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_copy(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = A[i][j];
        }
    }
    return result;
}

/** \brief Copies a sub-matrix
 *
 * \param[in] A Input matrix
 * \param[in] si Start of first index, \f$ s_i \f$
 * \param[in] ei End of first index, \f$ e_i \f$
 * \param[in] sj Start of second index, \f$ s_j \f$
 * \param[in] ej End of second index, \f$ e_j \f$
 * \param[in] result Matrix to store the result
 * \return Extracted matrix \f$ A_{s_i:e_i,s_j:e_j} \f$
 *
 */

MATRIX mat_xcopy(MATRIX A, int si, int ei, int sj, int ej, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(si<0 || sj<0 || ei>n || ej>m) mat_error(MAT_SIZEMISMATCH);
    if(result== NULL) if((result = mat_creat(ei-si, ej-sj, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=si; i<ei; ++i)
    {
        for(j=sj; j<ej; ++j)
        {
            result[i-si][j-sj] = A[i][j];
        }
    }
    return result;
}

/** \brief Copies a sub-matrix
 *
 * \param[in] A11 Input matrix, \f$ A_{11} \f$
 * \param[in] A12 Input matrix, \f$ A_{12} \f$
 * \param[in] A21 Input matrix, \f$ A_{21} \f$
 * \param[in] A22 Input matrix, \f$ A_{22} \f$
 * \param[in] result Matrix to store the result
 * \return Block matrix \f$ \left[\begin{array}{cc} A_{11} & A_{12}\\ A_{21} & A_{22} \end{array}\right] \f$
 *
 */

MATRIX mat_xjoin(MATRIX A11, MATRIX A12, MATRIX A21, MATRIX A22, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A11)+MatCol(A12);
    n = MatRow(A11)+MatRow(A21);
    if(result== NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A11); ++i)
    {
        for(j=0; j<MatCol(A11); ++j)
        {
            result[i][j] = A11[i][j];
        }
    }
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A12); ++i)
    {
        for(j=0; j<MatCol(A12); ++j)
        {
            result[i][j+MatCol(A11)] = A12[i][j];
        }
    }
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A21); ++i)
    {
        for(j=0; j<MatCol(A21); ++j)
        {
            result[i+MatRow(A21)][j] = A21[i][j];
        }
    }
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A22); ++i)
    {
        for(j=0; j<MatCol(A22); ++j)
        {
            result[i+MatRow(A11)][j+MatCol(A22)] = A22[i][j];
        }
    }
    return result;

}

/** \brief Copies a row from a matrix
 *
 * \param[in] A Input matrix
 * \param[in] rowa Source row
 * \param[in] rowb Destination row
 * \param[in] result Matrix to store the result
 * \return Copied matrix
 *
 */

MATRIX mat_rowcopy(MATRIX A, int rowa, int rowb, MATRIX result)
{
    int i, n;
    n = MatCol(result);
    for(i=0; i<n; ++i)
    {
        result[rowb][i] = A[rowa][i];
    }
    return result;
}

/** \brief Copies a column from a matrix
 *
 * \param[in] A Input matrix
 * \param[in] cola Source column
 * \param[in] colb Destination column
 * \param[in] result Matrix to store the result
 * \return Copied matrix
 *
 */

MATRIX mat_colcopy(MATRIX A, int cola, int colb, MATRIX result)
{
    int i, n;
    n = MatRow(result);
    for(i=0; i<n; ++i)
    {
        result[i][colb] = A[i][cola];
    }
    return result;
}

/** \brief Gets matrix data from opened file
 *
 * \param[in] A Matrix to store the data
 * \param[in] fp Pointer to opened file
 * \return Number of elements copied
 *
 */

int mat_fgetmat(MATRIX A, MAT_FILEPOINTER fp)
{
    int i, j, k=0, m, n;
    m = MatCol(A);
    n = MatRow(A);
    for(i=0; i<n; ++i)
#if mtype_n == 0
        for(j=0; j<m; ++j) k += fscanf(fp, "%f", &A[i][j]);
#elif mtype_n == 1
        for(j=0; j<m; j++) k += fscanf(fp, "%lf", &A[i][j]);
#endif
    return k;
}

/** \brief Creates a diagonal matrix from a 1-d matrix
 *
 * \param[in] diag_vals Input 1-d diagonal value matrix
 * \param[in] result Matrix to store the result
 * \return Diagonal matrix
 *
 */

MATRIX mat_creat_diag(MATRIX diag_vals, MATRIX result)
{
    int i, sz;
    if(MatCol(diag_vals)==1) sz = MatRow(diag_vals);
    else sz = MatCol(diag_vals);
    if(result==NULL) if((result = mat_creat(sz, sz, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    {
        if(MatCol(diag_vals)==1)
            for(i=0; i<sz; ++i) result[i][i] = diag_vals[i][0];
        else
            for(i=0; i<sz; ++i) result[i][i] = diag_vals[0][i];
        return result;
    }
}

