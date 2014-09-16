#include "matrix.h"


/** \brief Computes fast Fourier transform
 *
 * \param[in] C Complex data matrix stack
 * \param[in] dir FFT direction (ROWS/COLS)
 * \param[in] result Matrix stack to store the result
 * \return Transformed matrix stack
 *
 */

MATSTACK mat_fft2(MATSTACK c, int dir, MATSTACK result)
{
    int i, j, m, n, mm, nn, twopm;
    MATRIX real = NULL, imag = NULL, realt = NULL, imagt = NULL;
    m = MatCol(c[0]);
    n = MatRow(c[0]);
    if(result==NULL) if((result = matstack_creat(2))==NULL) matstack_error(MATSTACK_MALLOC);

    real = mat_copy(c[0], NULL);
    imag = mat_copy(c[1], NULL);

    if(real==NULL || imag==NULL)
        return(matstack_error(MATSTACK_MALLOC));
    if(!__mat_powerof2(m, &mm, &twopm) || twopm!=m)
        gen_error(GEN_MATH_ERROR);
    if(!__mat_powerof2(n, &nn, &twopm) || twopm!=n)
        gen_error(GEN_MATH_ERROR);

    for(j=0; j<n; ++j) __mat_fft(dir, mm, real[j], imag[j]);
    realt = mat_tran(real, NULL);
    imagt = mat_tran(imag, NULL);

    for(i=0; i<m; ++i) __mat_fft(dir, nn,realt[i],imagt[i]);
    mat_free(real);
    mat_free(imag);

    result[0] = mat_tran(realt, result[0]);
    result[1] = mat_tran(imagt, result[1]);
    mat_free(realt);
    mat_free(imagt);
    return result;
}

/**< \cond HIDDEN_SYMBOLS */

int __mat_fft(int dir, int m, mtype *x,mtype *y)
{
    long nn, i, i1, j, k, i2, l, l1, l2;
    mtype c1, c2, tx, ty, t1, t2, u1, u2, z;

    nn = 1;
    for(i=0; i<m; ++i) nn *= 2;
    i2 = nn >> 1;
    j = 0;
    for(i=0; i<nn-1; ++i)
    {
        if(i<j)
        {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while(k<=j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for(l=0; l<m; ++l)
    {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for(j=0; j<l1; ++j)
        {
            for(i=j; i<nn; i+=l2)
            {
                i1 = i+l1;
                t1 = u1*x[i1]-u2*y[i1];
                t2 = u1*y[i1]+u2*x[i1];
                x[i1] = x[i]-t1;
                y[i1] = y[i]-t2;
                x[i] += t1;
                y[i] += t2;
            }
            z = u1*c1-u2*c2;
            u2 = u1*c2+u2*c1;
            u1 = z;
        }
        c2 = (mtype)sqrt((1.0f-c1)/2.0f);
        if(dir == MAT_FFT2_FORWARD) c2 = -c2;
        c1 = (mtype)sqrt((1.0f+c1)/2.0f);
    }

    if(dir==MAT_FFT2_BACKWARD) for(i=0; i<nn; i++)
    {
        x[i] /= (mtype)nn;
        y[i] /= (mtype)nn;
    }
    return 0;
}

/**< \endcond */
