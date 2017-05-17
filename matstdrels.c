#include "matrix.h"

/** \brief Checks if greater than zero
 *
 * \param[in] a Input value
 * \return int \f$ a>0\f$
 *
 */

int gen_gt(mtype a)
{
    if(a>0)
        return 1;
    else
        return 0;
}

/** \brief Checks if less than zero
 *
 * \param[in] a Input value
 * \return int \f$ a<0\f$
 *
 */

int gen_lt(mtype a)
{
    if(a<0)
        return 1;
    else
        return 0;
}

/** \brief Checks if equals to zero
 *
 * \param[in] a Input value
 * \return int \f$ a==0\f$
 *
 */

int gen_eq(mtype a)
{
    if(a==0)
        return 1;
    else
        return 0;
}

