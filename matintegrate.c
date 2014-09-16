#include "matrix.h"


/** \brief Computes trapezoid integration
 *
 * \param[in] func Function \f$ f \left (\cdot \right) \f$ to integrate
 * \param[in] n Number of subdivisions
 * \param[in] lower Lower Limit
 * \param[in] upper Upper Limit
 * \return \f$ \int _a^bf\left ( x \right )dx \f$
 *
 */

mtype mat_int_trapezoid(mtype (*func)(mtype), int n, mtype lower, mtype upper)
{
    mtype h = (upper-lower)/((mtype)n);
    mtype sum = 0.0;
    int i;
    for(i=1; i<n; ++i)
        sum += func(lower+h*i);
    return h/2.0*(func(lower)+func(upper)+2.0*sum);
}

/** \brief Computes Simpson's integration
 *
 * \param[in] func Function \f$ f \left (\cdot \right) \f$ to integrate
 * \param[in] n Number of subdivisions
 * \param[in] lower Lower Limit
 * \param[in] upper Upper Limit
 * \return \f$ \int _a^bf\left ( x \right )dx \f$
 *
 */

mtype mat_int_simpson(mtype (*func)(mtype), int n, mtype lower, mtype upper)
{
    mtype h = (upper-lower)/((mtype)n);
    mtype sum1 = 0.0;
    mtype sum2 = 0.0;
    int i;
    for(i=0; i<n; ++i)
        sum1 += func(lower+h*i+h/2.0);
    for(i=1; i<n; ++i)
        sum2 += func(lower+h*i);
    return h/6.0*(func(lower)+func(upper)+4.0*sum1+2.0*sum2);
}


/** \brief Computes Gauss quadrature integration
 *
 * \param[in] func Function \f$ f \left (\cdot \right) \f$ to integrate
 * \param[in] n Number of subdivisions
 * \param[in] lower Lower Limit
 * \param[in] upper Upper Limit
 * \return \f$ \int _a^bf\left ( x \right )dx \f$
 *
 */

mtype mat_int_qadrat(mtype(*func)(mtype), mtype lower, mtype upper)
{
    mtype f0, f2, f3, f5, f6, f7, f9, f14, hmin, hmax, re, ae, result, *x, tmp0 = 0.0f;
    x= &tmp0;
    hmax=(upper-lower)/16.0f;
    if(hmax==0.0f) return (mtype)0.0f;
    re=(mtype)Eps ;
    ae=(mtype)(2.0f*eps/fabs(upper-lower));
    hmin=(mtype)fabs(upper-lower)*re;
    *x=lower ;
    f0=(*func)(*x);
    *x=lower+hmax;
    f2=(*func)(*x) ;
    *x=lower+2.0f*hmax;
    f3=(*func)(*x) ;
    *x=lower+4.0f*hmax;
    f5=(*func)(*x);
    *x=lower+6.0f*hmax;
    f6=(*func)(*x);
    *x=lower+8.0f*hmax;
    f7=(*func)(*x) ;
    *x=upper-4.0f*hmax;
    f9=(*func)(*x);
    *x=upper ;
    f14=(*func)(*x);
    result = __mat_lint(x, func, lower, upper, f0, f2, f3, f5, f6, f7, f9, f14, hmin, hmax, re, ae)*16.0f;
    return result;
}

/**< \cond HIDDEN_SYMBOLS */
mtype __mat_lint(mtype *x, mtype (*func)(mtype), mtype x0, mtype xn, mtype f0, mtype f2, mtype f3, mtype f5, mtype f6, mtype f7, mtype f9, mtype f14, mtype hmin, mtype hmax, mtype re, mtype ae)
{
    mtype v,w,h,xm,f1,f4,f8,f10,f11,f12, f13;
    xm=(x0+xn)/2.0f;
    h=(xn-x0)/32.0f;
    *x=xm+4.0f*h;
    f8=(*func)(*x);
    *x=xn-4.0f*h;
    f11=(*func)(*x);
    *x=xn-2.0f*h;
    f12=(*func)(*x);
    v=(mtype)(0.330580178199226*f7+0.173485115707338*(f6+f8)+
              0.321105426559972*(f5+f9)+0.135007708341042*(f3+f11)+
              0.165714514228223*(f2+f12)+0.393971460638127e-1*(f0+f14));
    *x=x0+h;
    f1=(*func)(*x);
    *x=xn-h;
    f13=(*func)(*x);
    w=(mtype)(0.260652434656970*f7+0.239063286684765*(f6 + f8) +
              0.263062635477467*(f5+f9)+0.218681931383057*(f3+f11)+
              0.275789764664284e-1*(f2+f12)+0.105575010053846*(f1+f13)+
              0.157119426059518e-1*(f0+f14));
    if(fabs(v-w)<fabs(w)*re+ae || fabs(h)<hmin)
        return h*w;
    else
    {
        *x=x0+6.0f*h;
        f4=(*func)(*x);
        *x=xn-6.0f*h;
        f10=(*func)(*x);
        v=(mtype)(0.245673430093324*f7+0.255786258286921*(f6+f8)+
                  0.228526063690406*(f5+f9)+0.500557131525460e-1*(f4+f10)+
                  0.177946487736780*(f3+f11)+0.584014599347449e-1*(f2+f12)+
                  0.874830942871331e-1*(f1+f13)+0.189642078648079e-1*(f0+f14));
        return ((fabs(v-w)<fabs(v)*re+ae) ? h*v :
                (__mat_lint(x,func,x0,xm,f0,f1,f2,f3,f4,f5,f6,f7,hmin,hmax,re,ae)-
                 __mat_lint(x,func,xn,xm,f14,f13,f12,f11,f10,f9, f8, f7,hmin,hmax,re,ae)));
    }
}
/**< \endcond */
