#include <malloc.h>
#include "matrix.h"


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

MATSTACK __matstack_creat(int length)
{
    MATSTACK matrixstack;
    int i;
    if((matrixstack = (MATRIX *)malloc(sizeof(MATRIX)*(1+length)))==NULL)
        return (matstack_error( MATSTACK_MALLOC));

    for(i=0; i<=length; ++i)
    {
        matrixstack[i]= NULL;
    }
    matrixstack[0] = mat_creat(1,1,0);
    matrixstack[0][0][0] = length;
    return (matrixstack+1);
}

MATSTACK matstack_creat(int length)
{
    MATSTACK A;
    if((A =__matstack_creat(length))!=NULL)
    {
        return A;
    }
    else
        return NULL;
}

MATSTACK matstack_append(MATSTACK s, MATRIX a)
{
    MATSTACK m;
    int i, n;
    if(a==NULL) return s;
    n = MatStacklength(s);
    m = matstack_creat(n+1);
    for(i=0; i<n; ++i) m[i] = s[i];
    m[n] = a;
    return m;
}

int matstack_free(MATSTACK A)
{
    int i, n;
    if(A==NULL) return 0;
    n = MatStacklength(A);
    A = A-1;
    for(i=0; i<n; ++i)
        mat_free(A[i]);
    A = NULL;
    return 1;
}

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


MATRIX mat_fill_type(MATRIX A, int type)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    switch(type)
    {
    case UNDEFINED:
        break;
    case ZERO_MATRIX:
    case UNIT_MATRIX:
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
            for(j=0; j<m; ++j)
            {
                if(type==UNIT_MATRIX)
                {
                    if(i==j)
                    {
                        A[i][j] = 1.0;
                        continue;
                    }
                }
                A[i][j] = 0.0;
            }
        break;
    case ONES_MATRIX:
        #pragma omp parallel for private(j)
        for(i=0; i<n; ++i)
            for(j=0; j<m; ++j) A[i][j] = 1.0;
        break;
    }
    return A;
}

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

INT_VECTOR __int_vec_creat(int length)
{
    INT_VECTOR int_vector;
    if((int_vector = (int*)malloc(sizeof(int)*(length+1)))==NULL)
        return (int_vec_error(INT_VEC_MALLOC));
    int_vector[0] = length;
    return (int_vector+1);
}

INT_VECTOR int_vec_creat(int length, int type)
{
    INT_VECTOR A;
    if((A =__int_vec_creat(length))!=NULL)
    {
        return (int_vec_fill_type(A, type));
    }
    else
        return NULL;
}

INT_VECTOR int_vec_fill(INT_VECTOR A, int val)
{
    int i, n;
    n = Int_VecLen(A);
    for(i=0; i<n; ++i) A[i] = val;
    return A;
}

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

int int_vec_free(INT_VECTOR A)
{
    if(A==NULL) return 0;
    Int_VecLen(A) = 0;
    free(A-1);
    A = NULL;
    return 1;
}

INT_VECSTACK __int_vecstack_creat(int length)
{
    INT_VECSTACK int_vecstack;
    int i;
    if((int_vecstack = (INT_VECTOR *)malloc( sizeof(INT_VECTOR)*(1+length)))==NULL)
        return (int_vecstack_error(INT_VECSTACK_MALLOC));

    for(i=0; i<=length; ++i)
    {
        int_vecstack[i]= NULL;
    }
    int_vecstack[0] = int_vec_creat(1,0);
    int_vecstack[0][0] = length;
    return (int_vecstack+1);
}

INT_VECSTACK int_vecstack_creat(int length)
{
    INT_VECSTACK A;
    if((A =__int_vecstack_creat(length))!=NULL)
    {
        return A;
    }
    else
        return NULL;
}

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

BAYES_MODEL pat_bayes_model_creat(void)
{
    BAYES_MODEL a;
    a = (BAYES_MODEL)malloc(sizeof(bayes_model));
    a->class_priors = NULL;
    a->class_covars = NULL;
    a->class_labels = NULL;
    a->class_means = 0;
    a->num_of_classes = 0;
    a->num_of_features = 0;
    return a;
}

int pat_bayes_model_free(BAYES_MODEL a)
{
    matstack_free(a->class_means);
    matstack_free(a->class_covars);
    int_vec_free(a->class_labels);
    free(a);
    a = NULL;
    return 1;
}

PERCEPTRON pat_perceptron_creat(void)
{
    PERCEPTRON a;
    a = (PERCEPTRON)malloc(sizeof(perceptron));
    a->class_weights = NULL;
    a->class_labels = NULL;
    a->istrained = 0;
    a->num_of_classes = 0;
    a->num_of_features = 0;
    a->num_of_iterations = 20;
    return a;
}

int pat_perceptron_free(PERCEPTRON a)
{
    if(a->istrained!=0) mat_free(a->class_weights);
    int_vec_free(a->class_labels);
    free(a);
    a = NULL;
    return 1;
}

MATVEC_DPOINTER matvec_creat(void)
{
    MATVEC_DPOINTER a;
    a = (MATVEC_DPOINTER)malloc(sizeof(void *)*2);
    a[0] = NULL;
    a[1] = NULL;
    return a;
}

int matvec_free(MATVEC_DPOINTER a)
{
    mat_free((MATRIX)a[0]);
    int_vec_free((INT_VECTOR)a[1]);
    free(a);/* to be edited later */
    return 1;
}

INT_VECTOR int_vec_append(INT_VECTOR A, int i)
{
    int l;
    l = Int_VecLen(A);
    A--;
    if((A = (int *)realloc(A, sizeof(int)*(l+2)))==NULL) return int_vec_error(INT_VEC_MALLOC);
    A[0]= l+1;
    ++A;
    A[l] = i;
    return A;
}

INT_VECTOR int_vec_copy(INT_VECTOR a)
{
    int i, m;
    INT_VECTOR b;
    m = Int_VecLen(a);
    if((b = int_vec_creat(m, UNDEFINED))==NULL)
        return NULL;
    for(i=0; i<m; ++i) b[i] = a[i];
    return b;
}

MATRIX mat_copy(MATRIX A, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
        for(j=0; j<m; ++j)
        {
            result[i][j] = A[i][j];
        }
    return result;
}

MATRIX mat_xcopy(MATRIX A, int si, int ei, int sj, int ej, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(si<0 || sj<0 || ei>n || ei>m) mat_error(MAT_SIZEMISMATCH);
    if(result== NULL) if((result = mat_creat(ei-si, ej-sj, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=si; i<ei; ++i)
        for(j=sj; j<ej; ++j)
        {
            result[i-si][j-sj] = A[i][j];
        }
    return result;
}

MATRIX mat_xjoin(MATRIX A11, MATRIX A12, MATRIX A21, MATRIX A22, MATRIX result)
{
    int i, j, m, n;
    m = MatCol(A11)+MatCol(A12);
    n = MatRow(A11)+MatRow(A21);
    if(result== NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A11); ++i)
        for(j=0; j<MatCol(A11); ++j)
        {
            result[i][j] = A11[i][j];
        }
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A12); ++i)
        for(j=0; j<MatCol(A12); ++j)
        {
            result[i][j+MatCol(A11)] = A12[i][j];
        }
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A21); ++i)
        for(j=0; j<MatCol(A21); ++j)
        {
            result[i+MatRow(A21)][j] = A21[i][j];
        }
    #pragma omp parallel for private(j)
    for(i=0; i<MatRow(A22); ++i)
        for(j=0; j<MatCol(A22); ++j)
        {
            result[i+MatRow(A11)][j+MatCol(A22)] = A22[i][j];
        }
    return result;

}

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

MATRIX mat_extractrows(MATRIX A, INT_VECTOR rows, MATRIX result)
{
    int i, j, m, n;
    n = MatCol(A);
    m = Int_VecLen(rows);
    if(result==NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n; ++j)
        {
            result[i][j] = A[rows[i]][j];
        }
    }
    return result;
}

MATRIX mat_extractcols(MATRIX A, INT_VECTOR cols, MATRIX result)
{
    int i, j, m, n;
    m = MatRow(A);
    n = Int_VecLen(cols);
    if(result==NULL) if((result = mat_creat(m, n, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n; ++j)
        {
            result[i][jj] = A[i][cols[j]];
        }
    }
    return result;
}

int mat_fgetmat(MATRIX A, FILEPOINTER fp)
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

