/* 
 * File:   MDVec.cc
 * Author: tingnan
 * 
 * Created on February 8, 2011, 8:30 PM
 */

#include "md_vec.h"

/* following are some vector operations used in this program, VecI is integer vector,
 * used for cell operations, and VecR is real vector, used for discrib ball's position,
 * velocity, and so on. Some operations are none-standard. */

IntVec::IntVec ()
{
}

RealVec::RealVec ()
{
}

RealVec::RealVec (double a, double b)
{
    x = a;
    y = b;
}

RealVec::RealVec (double a, double b, double c)
{
    x = a;
    y = b;
    z = c;
}

RealVec
CrossProduct (const RealVec & a, const RealVec & b)
{
    RealVec c;

    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}
