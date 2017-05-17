#include "matrix.h"


/** \brief Computes addition function
 *
 * \param[in] x
 * \param[in] y
 * \return \f$ x+y \f$
 *
 */

mtype __mat_addfunc(mtype x, mtype y)
{
    return (x+y);
}

/** \brief Computes subtraction function
 *
 * \param[in] x
 * \param[in] y
 * \return \f$ x-y \f$
 *
 */

mtype __mat_subfunc(mtype x, mtype y)
{
    return (x-y);
}

/** \brief Computes multiplication function
 *
 * \param[in] x
 * \param[in] y
 * \return \f$ xy \f$
 *
 */

mtype __mat_mulfunc(mtype x, mtype y)
{
    return (x*y);
}

/** \brief Computes division function
 *
 * \param[in] x
 * \param[in] y
 * \return \f$ \frac{x}{y} \f$
 *
 */

mtype __mat_divfunc(mtype x, mtype y)
{
    return (x/y);
}

/** \brief Computes square function
 *
 * \param[in] x
 * \return \f$ x^{2} \f$
 *
 */

mtype __mat_sqrfunc(mtype x)
{
    return (x*x);
}

/** \brief Computes square root function
 *
 * \param[in] x
 * \return \f$ \sqrt{x} \f$
 *
 */

mtype __mat_sqrtfunc(mtype x)
{
    return sqrt(x);
}

/** \brief Computes Huber weight function
 *
 * \param[in] x
 * \param[in] k
 * \return \f$ \begin{cases} 1, & \text{for } |x| \le k, \\ \frac{k}{|x|}, & \text{otherwise.}\end{cases} \f$
 *
 */

mtype __mat_huber_wt(mtype x, mtype k)
{
    if(fabs(x)<= k) return 1.0;
    else return (mtype)(k/fabs(x));
}

/** \brief Computes bisquare weight function
 *
 * \param[in] x
 * \param[in] k
 * \return \f$ \begin{cases} \left ( 1-\left ( \frac{x}{k} \right )^{2} \right )^{2}, & \text{for } |x| \le k, \\ 0, & \text{otherwise.}\end{cases} \f$
 *
 */

mtype __mat_bisquare_wt(mtype x, mtype k)
{
    mtype a;
    if(fabs(x)<= k)
    {
        a = x/k;
        a = 1 - (a*a);
        return (a*a);
    }
    else return 0.0;
}

/** \cond HIDDEN_SYMBOLS */
__inline mtype __huber_wt(mtype x, mtype k)
{
    if(fabs(x)<= k) return 1.0;
    else return (mtype)(k/fabs(x));
}

__inline mtype __bisquare_wt(mtype x, mtype k)
{
    mtype a;
    if(fabs(x)<= k)
    {
        a = x/k;
        a = 1 - (a*a);
        return (a*a);
    }
    else return 0.0;
}
/** \endcond */

/** \brief Computes inverse hyperbolic sine function
 *
 * \param[in] x
 * \return \f$ \sinh^{-1}\left(x\right) \f$
 *
 */

mtype __mat_arcsinh(mtype x)
{
    mtype y;
    if (fabs(x) > 1.0e10)
        return (mtype)((x > 0.0) ? 0.69314718055995+log(fabs(x)) :-0.69314718055995+log(fabs(x)));
    else
    {
        y=x*x;
        return (mtype)((x == 0.0f) ? 0.0f : ((x > 0.0f) ?
                                    __mat_logplusone((float)(fabs(x)+y/(1.0f+sqrt(1.0f+y)))) :
                                    -__mat_logplusone((float)(fabs(x)+y/(1.0f+sqrt(1.0f+y))))));
    }
}

/** \brief Computes inverse hyperbolic cosine function
 *
 * \param[in] x
 * \return \f$ \cosh^{-1}\left(x\right) \f$
 *
 */

mtype __mat_arccosh(mtype x)
{
    return (mtype)((x <= 1.0) ? 0.0 : ((x > 1.0e10) ?
                                0.69314718055995+log(x) :
                                log(x+sqrt((x-1.0)*(x+1.0)))));
}

/** \brief Computes inverse hyperbolic tangent function
 *
 * \param[in] x
 * \return \f$ \tanh^{-1}\left(x\right) \f$
 *
 */

mtype __mat_arctanh(mtype x)
{
    mtype ax;
    if (fabs(x) >= 1.0)return (mtype)((x > 0.0) ? DOUBLE_MAX : -DOUBLE_MAX);
    else
    {
        ax=(mtype)fabs(x);
        return (mtype)((x == 0.0) ? 0.0 : ((x > 0.0) ? 0.5*__mat_logplusone(2.0f*ax/(1.0f-ax)) :
                                    -0.5*__mat_logplusone(2.0f*ax/(1.0f-ax))));
    }
}

/** \brief Computes logarithm plus one function
 *
 * \param[in] x
 * \return \f$ \log\left(1+x\right) \f$
 *
 */

mtype __mat_logplusone(mtype x)
{
    mtype y,z;
    if(x==0.0) return 0.0;
    else if(x<-0.2928 || x>0.4142) return (mtype)log(1.0f+x);
    else
    {
        z=x/(x+2.0f);
        y=z*z;
        return (mtype)(z*(2.0+y*(0.66666666663366+y*(0.400000001206045+y*(0.285714091590488+y*(0.22223823332791+y*(0.1811136267967+y*0.16948212488)))))));
    }
}

/** \brief Rounds a number away from zero
 *
 * \param[in] x Input value
 * \return \f$ \textrm{sgn}(x) \left\lfloor \left| x \right| + 0.5 \right\rfloor \f$
 *
 */

mtype __mat_round_away_zero(mtype x)
{
    if(x>0) return ceil(x);
    else return floor(x);
}

/** \brief Rounds a number towards zero
 *
 * \param[in] x Input value
 * \return \f$ \textrm{sgn}(x) \left\lceil \left| x \right| - 0.5 \right\rceil \f$
 *
 */

mtype __mat_round_towards_zero(mtype x)
{
    if(x>0) return floor(x);
    else return ceil(x);
}


/** \brief Computes Huber weight function element-wise on a matrix
 *
 * \param[in] A Input matrix
 * \param[in] k Huber parameter
 * \return \f$ \mathbf{B},\, b_{ij}=f_k\left(a_{ij}\right) \f$ where \f$ f_k \f$ is the Huber weight function
 *
 */

MATRIX mat_huber_wt(MATRIX A, mtype k, mtype sigma, MATRIX result)
{
    int	i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
        return mat_error(MAT_MALLOC);;

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = __huber_wt(A[i][j]/sigma, k);
        }
    }
    return(result);
}

/** \brief Computes bisquare weight function element-wise on a matrix
 *
 * \param[in] A Input matrix
 * \param[in] k Bisquare parameter
 * \return \f$ \mathbf{B},\, b_{ij}=f_k\left(a_{ij}\right) \f$ where \f$ f_k \f$ is the biquare weight function
 *
 */

MATRIX mat_bisquare_wt(MATRIX A, mtype k, mtype sigma, MATRIX result)
{
    int	i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);
    if(result==NULL) if((result = mat_creat(n, m, UNDEFINED))==NULL)
        return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = __bisquare_wt(A[i][j]/sigma, k);
        }
    }
    return(result);
}


/** \brief Computes a given function element-wise on a matrix
 *
 * \param[in] A Input matrix
 * \param[in] f Given function
 * \return \f$ \mathbf{B},\, b_{ij}=f\left(a_{ij}\right) \f$
 *
 */

MATRIX mat_gfunc(MATRIX A, mtype (*pt2func)(mtype), MATRIX result)
{
    int	i, j, m, n;
    m = MatCol(A);
    n = MatRow(A);

    if(result== NULL|| result==A) if((result = mat_creat(n, m, UNDEFINED))==NULL)
            return mat_error(MAT_MALLOC);

    #pragma omp parallel for private(j)
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result[i][j] = (*pt2func)(A[i][j]);
        }
    }
    return(result);
}

