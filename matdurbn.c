#include "matrix.h"


/** \brief Runs Levinson-Durbin algorithm
 *
 * \param[in] R Input \f$ n^th \f$ correlation matrix \f$ (n+1) \times 1 \f$
 * \param[in] result Matrix to store the result
 * \return \f$ X \f$ where \f$ \tilde{R}X = B \f$ , \f$ \tilde{R} = \left[ \begin{array}{cccc} R[0][0] & R[1][0] & \cdots & R[n-1][0]\\  R[1][0] & R[0][0] & \cdots & R[n-2][0]\\  \vdots & \vdots & \ddots & \vdots \\  R[n-1][0] & R[n-2][0] & \cdots & R[0][0] \end{array} \right] \f$ and \f$ B = \left[ \begin{array}{cccc} R[1][0] & R[2][0] & \cdots & R[n][0] \end{array} \right] \f$
 *
 */

MATRIX mat_durbin(MATRIX R, MATRIX result)
{
    int i, i1, j, ji, p;
    MATRIX W, E, K, A;
    p = MatRow(R) - 1;
    W = mat_creat(p+2, 1, UNDEFINED);
    E = mat_creat(p+2, 1, UNDEFINED);
    K = mat_creat(p+2, 1, UNDEFINED);
    A = mat_creat(p+2, p+2, UNDEFINED);
    W[0][0] = R[1][0];
    E[0][0] = R[0][0];
    for(i=1; i<=p; ++i)
    {
        K[i][0] = W[i-1][0]/E[i-1][0];
        E[i][0] = E[i-1][0]*(1.0f - K[i][0]*K[i][0]);
        A[i][i] = -K[i][0];
        i1 = i-1;
        if(i1>=1)
        {
            for(j=1; j<=i1; ++j)
            {
                ji = i-j;
                A[j][i] = A[j][i1]-K[i][0]*A[ji][i1];
            }
        }
        if(i!=p)
        {
            W[i][0] = R[i+1][0];
            for(j=1; j<=i; ++j)
                W[i][0] += A[j][i]*R[i-j+1][0];
        }
    }
    if(result==NULL) if((result = mat_creat(p, 1, UNDEFINED))==NULL) mat_error(MAT_MALLOC);
    for(i=0; i<p; ++i) result[i][0] = -A[i+1][p];

    mat_free(A);
    mat_free(W);
    mat_free(K);
    mat_free(E);
    return result;
}

/** \brief Runs Levinson-Durbin algorithm
 *
 * \param[in] A Input correlation matrix \f$ A = \left[ \begin{array}{cccc} r_{0} & r_{1} & \cdots & r_{n-1}\\  r_{1} & r_{0} & \cdots & r_{n-2}\\  \vdots & \vdots & \ddots & \vdots \\  r_{n-1} & r_{n-2} & \cdots & r_{0} \end{array} \right] \f$
 * \param[in] B Input correlation matrix\f$ B = \left[ \begin{array}{c} r_{1} \\ r_{2} \\ \cdots \\ r_{n} \end{array} \right] \f$
 * \param[in] result Matrix to store the result
 * \return \f$ X \f$ where \f$ RX = B \f$
 *
 */

MATRIX mat_lsolve_durbin(MATRIX A, MATRIX B, MATRIX result)
{
    MATRIX R;
    int i, n;
    n = MatRow(A);
    R = mat_creat(n+1, 1, UNDEFINED);
    for(i=0; i<n; ++i)
    {
        R[i][0] = A[i][0];
    }
    R[n][0] = B[n-1][0];
    result = mat_durbin(R, result);
    mat_free(R);
    return result;
}

/** \brief Computes QR decomposition
 *
 * \param[in] A Input matrix
 * \param[in] qr Matrix stack to store result
 * \return Output QR Matrix stack
 *
 */

MATSTACK mat_qr(MATRIX A, MATSTACK qr)
{
    mtype mag, alpha;
    MATRIX u = NULL, v = NULL, P = NULL, I = NULL, tmp = NULL, tmp1 = NULL, tmp2 = NULL, tmp3 = NULL;
    int m, n, i, j, f = 1, s;
    m = MatRow(A);
    n = MatCol(A);
    if(qr==NULL)
    {
        if((qr = matstack_creat(2))==NULL) matstack_error(MATSTACK_MALLOC);
        qr[0] = mat_creat(m, m, UNDEFINED);
        qr[1] = mat_creat(m, n, UNDEFINED);
    }
    u = mat_creat(m, 1, UNDEFINED);
    v = mat_creat(m, 1, UNDEFINED);
    I = mat_creat(m, m, UNIT_MATRIX);
    qr[0] = mat_copy(I, qr[0]);
    qr[1] = mat_copy(A, qr[1]);
    tmp3 = qr[1];
    s = (m>n)?n:m;
    for(i=0; i<s; ++i)
    {
        for(j=0; j<m; ++j)
        {
            u[j][0] = 0.0;
            v[j][0] = 0.0;
        }
        mag = 0.0;
        for(j=i; j<m; ++j)
        {
            u[j][0] = tmp3[j][i];
            mag += u[j][0]*u[j][0];
        }
        mag = sqrt(mag);
        alpha = u[i][0]<0?mag:-mag;
        mag = 0.0;
        for(j=i; j<m; ++j)
        {
            v[j][0] = ((j==i)?(u[j][0]+alpha):u[j][0]);
            mag += v[j][0]*v[j][0];
        }
        mag = __mat_sqrtfunc(mag);
        if(mag<eps) continue;
        for(j=i; j<m; ++j) v[j][0]/= mag;
        tmp = mat_tran(v, tmp);
        P = mat_mul(v, tmp, P);
        P = mat_muls(P, -2.0, P);
        for(j=0; j<m; ++j) P[j][j] += 1.0;
        if(f>0)
        {
            tmp1 = mat_mul(P, qr[1], tmp1);
            tmp2 = mat_mul(qr[0], P, tmp2);
            tmp3 = tmp1;
        }
        else
        {
            qr[1] = mat_mul(P, tmp1, qr[1]);
            qr[0] = mat_mul(tmp2, P, qr[0]);
            tmp3 = qr[1];
        }
        f = f*-1;
    }
    if(f<0)
    {
        mat_free(qr[0]);
        mat_free(qr[1]);
        qr[0] = tmp2;
        qr[1] = tmp1;
    }
    else
    {
        mat_free(tmp1);
        mat_free(tmp2);
    }
    mat_free(u);
    mat_free(v);
    mat_free(P);
    mat_free(tmp);
    mat_free(I);
    return qr;
}

