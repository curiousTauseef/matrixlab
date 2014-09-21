#include "matrix.h"


MATSTACK mat_cheby_series_table, mat_legendre_series_table, mat_binom_series_table;

/** \brief Evaluates polynomial at a point
 *
 * \param[in] A Input polynomial matrix
 * \param[in] x Value at which to evaluate
 * \param[in] dir Direction (ROWS/COLS)
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_eval(MATRIX A, mtype x, int dir, MATRIX result)
{
    int m, n, i, j;
    mtype r;
    m = MatRow(A);
    n = MatCol(A);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for(j=m-1; j>=0; --j) r=r*x+A[j][i];
            result[0][i] = r;
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for(j=n-1; j>=0; --j) r=r*x+A[i][j];
            result[i][0] = r;
        }
    }
    return result;
}

/** \brief Computes derivative polynomial of a polynomial
 *
 * \param[in] A Input polynomial matrix
 * \param[in] dir Direction (ROWS/COLS)
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_diff(MATRIX A, int dir, MATRIX result)
{
    int m, n, i, j;
    m = MatRow(A);
    n = MatCol(A);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(m-1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            for(j=m-2; j>=0; --j) result[j][i] = (j+1)*A[j+1][i];
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, n-1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            for(j=n-2; j>=0; --j) result[i][j]= (j+1)*A[i][j+1];
        }
    }
    return result;
}

/** \brief Evaluates derivative polynomial at a point
 *
 * \param[in] A Input polynomial matrix
 * \param[in] x Value at which to evaluate the derivative
 * \param[in] dir Direction (ROWS/COLS)
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_diff_eval(MATRIX A, mtype x, int dir, MATRIX result)
{
    int m, n, i, j;
    mtype r;
    m = MatRow(A);
    n = MatCol(A);
    if(dir==0)
    {
        if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<n; ++i)
        {
            r = 0.0;
            for(j=m-2; j>=0; --j) r=r*x+A[j+1][i]*(j+1);
            result[0][i] = r;
        }
    }
    else
    {
        if(result==NULL) if((result = mat_creat(m, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<m; ++i)
        {
            r = 0.0;
            for(j=n-2; j>=0; --j) r=r*x+A[i][j+1]*(j+1);
            result[i][0] = r;
        }
    }
    return result;
}

/** \brief Adds two polynomials
 *
 * \param[in] A First input polynomial matrix
 * \param[in] B Second input polynomial matrix
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_add(MATRIX A, MATRIX B, MATRIX result)
{
    int i, na, nb;
    na = MatCol(A);
    nb = MatCol(B);
    if(na==nb) return mat_add(A, B, result);
    if(na<nb)
    {
        if(result==NULL) if((result = mat_creat(1, nb, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<na; ++i) result[0][i] = A[0][i]+B[0][i];
        for(i=na; i<nb; ++i) result[0][i] = B[0][i];
    }
    else
    {
        if(result==NULL) if((result = mat_creat(1, na, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        for(i=0; i<nb; ++i) result[0][i] = A[0][i]+B[0][i];
        for(i=nb; i<na; ++i) result[0][i] = A[0][i];
    }
    return result;
}

/** \brief Multiplies two polynomials
 *
 * \param[in] a First input polynomial matrix
 * \param[in] b Second input polynomial matrix
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_mul(MATRIX A, MATRIX B, MATRIX result)
{
    int m, n, i, j, s;
    m = MatCol(A);
    n = MatCol(B);
    s = m+n-1;
    if(result==NULL) if((result = mat_creat(1, s, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    result = mat_fill(result, 0.0);
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n; ++j)
        {
            result[0][i+j] += A[0][i]*B[0][j];
        }
    }
    return result;
}

/** \brief Divides two polynomials
 *
 * \param[in] A First input polynomial matrix
 * \param[in] B Second input polynomial matrix
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATSTACK mat_poly_div(MATRIX A, MATRIX B, MATSTACK result)
{
    MATRIX tmp = NULL;
    int m, n, i, j;
    mtype s;
    m = MatCol(A);
    n = MatCol(B);
    if(m>=n)
    {
        if(result==NULL)
        {
            if((result = matstack_creat(2))==NULL) matstack_error(MATSTACK_MALLOC);
            if((result[0] = mat_creat(1, m-n+1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
            if((result[1] = mat_creat(1, n-1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
        }
        tmp = mat_copy(A, NULL);
        i = m;
        s = B[0][n-1];
        while((i)>=n)
        {
            result[0][0][i-n] = tmp[0][i-1]/s;
            for(j=i-1; j>=(i-n); --j) tmp[0][j] -= B[0][j-i+n]*result[0][0][i-n];
            --i;
        }
        for(i=0; i<(n-1); ++i) result[1][0][i] = tmp[0][i];
        mat_free(tmp);
    }
    else
    {
        if(result==NULL)
        {
            if((result = matstack_creat(2))==NULL) matstack_error(MATSTACK_MALLOC);
            if((result[0] = mat_creat(1, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
            if((result[1] = mat_creat(1, n-1, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
        }
        result[0] = mat_fill(result[0], 0.0);
        for(i=0; i<m; ++i) result[1][0][i] = A[0][i];
    }
    return result;
}

/** \brief Multiplies a polynomial with a scalar
 *
 * \param[in] A Input polynomial matrix
 * \param[in] s Scalar
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_scale(MATRIX A, mtype s, MATRIX result)
{
    return mat_muls(A, s, result);
}

/** \brief Shifts a polynomial
 *
 * \param[in] A Input polynomial matrix
 * \param[in] s Scalar shift
 * \param[in] result Matrix to store the result
 * \return Output matrix
 *
 */

MATRIX mat_poly_shift(MATRIX A, int s, MATRIX result)
{
    int i, n;
    n = MatCol(A)+s;
    if(result==NULL) if((result = mat_creat(1, (n>0)?n:0, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    mat_fill(result, 0.0);
    for(i=(s>0)?s:0; i<n; ++i) result[0][i] = A[0][i-s];
    return result;
}

/** \brief Initializes the Chebyshev polynomial series
 *
 *
 */

void mat_cheby_init()
{
    mat_cheby_series_table = matstack_creat(512);
    mat_cheby_series_table[0] = mat_creat(1, 1, ONES_MATRIX);
    mat_cheby_series_table[1] = mat_creat(1, 2, ZERO_MATRIX);
    mat_cheby_series_table[1][0][1] = 1;
}

/** \brief Initializes the Legendre polynomial series
 *
 *
 */

void mat_legendre_init()
{
    mat_legendre_series_table = matstack_creat(512);
    mat_legendre_series_table[0] = mat_creat(1, 1, ONES_MATRIX);
    mat_legendre_series_table[1] = mat_creat(1, 2, ZERO_MATRIX);
    mat_legendre_series_table[1][0][1] = 1;
}

/** \brief Initializes the binomial series
 *
 *
 */

void mat_binom_init()
{
    mat_binom_series_table = matstack_creat(512);
    mat_binom_series_table[0] = mat_creat(1, 1, ONES_MATRIX);
    mat_binom_series_table[1] = mat_creat(1, 2, ONES_MATRIX);
}

/** \brief Computes the \f$ n^{th} \f$ Chebyshev polynomial
 *
 * \param[in] n Polynomial series index
 * \return Output polynomial matrix
 *
 */

MATRIX mat_cheby(int n)
{
    if(mat_cheby_series_table[n]==NULL)
    {
        MATRIX p = mat_cheby(n-1);
        MATRIX pp = mat_cheby(n-2);
        MATRIX npp = mat_poly_scale(pp, -1.0, NULL);
        MATRIX px2 = mat_poly_shift(p, 1, NULL);
        px2 = mat_poly_scale(px2, 2, px2);
        mat_cheby_series_table[n] = mat_poly_add(px2, npp, px2);
        mat_free(npp);
    }
    return mat_cheby_series_table[n];
}

/** \brief Computes the \f$ n^{th} \f$ Legendre polynomial
 *
 * \param[in] n Polynomial series index
 * \return Output polynomial matrix
 *
 */

MATRIX mat_legendre(int n)
{
    if(mat_legendre_series_table[n]==NULL)
    {
        MATRIX p = mat_legendre(n-1);
        MATRIX pp = mat_legendre(n-2);
        MATRIX npp = mat_poly_scale(pp, ((1.0/n)-1.0), NULL);
        MATRIX px2 = mat_poly_shift(p, 1, NULL);
        px2 = mat_poly_scale(px2, (2.0-(1.0/n)), px2);
        mat_legendre_series_table[n] = mat_poly_add(px2, npp, px2);
        mat_free(npp);
    }
    return mat_legendre_series_table[n];
}

/** \brief Computes a binomial co-efficient
 *
 * \param[in] n \f$ 1^{st} \f$ argument
 * \param[in] k \f$ 2^{nd} \f$ argument
 * \return \f$ \binom{n}{k} \f$
 *
 */

mtype mat_binom(int n, int k)
{
    int i;
    if(mat_binom_series_table[n]==NULL)
    {
        mat_binom_series_table[n] = mat_creat(1, n+1, UNDEFINED);
        mat_binom_series_table[n][0][0] = 1;
        mat_binom_series_table[n][0][n] = 1;
        for(i=1; i<n; ++i)
        {
            mat_binom_series_table[n][0][i] = mat_binom(n-1, i-1)+mat_binom(n-1, i);
        }
    }
    return mat_binom_series_table[n][0][k];
}

/** \brief Converts Chebyshev co-efficients to a single polynomial
 *
 * \param[in] coeffs Chebyshev polynomial co-efficient matrix
 * \param[in] result Matrix to store the result
 * \return Polynomial matrix
 *
 */

MATRIX mat_cheby_coeffs_to_poly(MATRIX coeffs, MATRIX result)
{
    int i, n;
    n = MatCol(coeffs);
    if(result==NULL) if((result = mat_creat(1, n, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    mat_fill(result, 0.0);
    result[0][0] = coeffs[0][0]/2;
    for(i=1; i<n; ++i)
    {
        MATRIX sp = mat_poly_scale(mat_cheby(i), coeffs[0][i], NULL);
        result = mat_poly_add(result, sp, result);
        mat_free(sp);
    }
    return result;
}

/** \brief Approximates a function using Chebyshev polynomials
 *
 * \param[in] f Function to approximate
 * \param[in] a Lower limit of domain of the function
 * \param[in] b Upper limit of domain of the function
 * \param[in] n Degree of the approximate polynomial
 * \param[in] result Matrix to store the result
 * \return Approximate polynomial matrix
 *
 */

MATRIX mat_cheby_approx(mtype (*f)(mtype), mtype a, mtype b, int n, MATRIX result)
{
    int i, j;
    mtype r, c, xs, sum, m;
    MATRIX ys = NULL, coeffs = NULL, tmp = NULL, ctmp = NULL;
    r = (b-a)/2.0;
    c = -(a+b)/2.0;
    ys = mat_creat(1, n+1, UNDEFINED);
    coeffs = mat_creat(1, n+1, UNDEFINED);
    for(i=0; i<=n; ++i)
    {
        xs = r*cos(MAT_PI*(i+0.5)/(n+1))-c;
        ys[0][i] = f(xs);
    }
    for(j=0; j<=n; ++j)
    {
        sum = 0.0;
        for(i=0; i<=n; ++i)
        {
            sum += cos(MAT_PI*j*(i+0.5)/(n+1))*ys[0][i];
        }
        coeffs[0][j] = sum*2.0/(n+1);
    }
    mat_free(ys);
    result = mat_cheby_coeffs_to_poly(coeffs, result);
    mat_free(coeffs);
    m = r;
    for(i=1; i<=n; ++i)
    {
        result[0][i] /= m;
        m *= r;
    }
    if(fabs(c)>eps)
    {
        tmp = result;
        result = mat_creat(1, n+1, ZERO_MATRIX);
        ctmp = mat_creat(1, n+1, UNDEFINED);
        ctmp[0][0] = 1;
        for(i=1; i<=n; ++i) ctmp[0][i] = ctmp[0][i-1]*c;
        mat_binom_init();
        for(j=0; j<=n; ++j)
        {
            for(i=j; i<=n; ++i)
            {
                result[0][j] += tmp[0][i]*mat_binom(i,j)*ctmp[0][i];
            }
            result[0][j] /= ctmp[0][j];
        }
        mat_free(tmp);
        mat_free(ctmp);
    }
    return result;
}

