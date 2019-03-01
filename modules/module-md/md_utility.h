/* 
 * File:   MDUtility.h
 * Author: tingnan
 *
 * Created on February 8, 2011, 10:24 PM
 */

#include "md_header.h"
#include "md_vec.h"
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>

#ifndef _MDUTILITY_H
#define	_MDUTILITY_H

std::string ftoa (double);
std::string int2str (int);

inline double
heaviside (double number)
{
    if (number >= 0)
        return 1.0;
    else
        return 0.0;

}

inline double
min (double a, double b)
{
    if (a < b)
        return a;
    else
        return b;
}

inline double
rsquare (RealVec r)
{
    return r.x * r.x + r.z * r.z + r.y * r.y;
}


int sign (double);

#endif /* MDUTILITY_H */
