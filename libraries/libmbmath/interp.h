/* $Header$ */

#ifndef INTERP_H
#define INTERP_H

#include <vector>
#include "ac/f2c.h"
#include "myassert.h"

/*
 * Compute spline coefficients
 */
extern void
spline(const std::vector<doublereal>& x, 
	const std::vector<doublereal>& y, 
	std::vector<doublereal>& b, 
	std::vector<doublereal>& c, 
	std::vector<doublereal>& d);

/*
 * Evaluate spline
*/
extern doublereal
seval(const doublereal& u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const std::vector<doublereal>& b,
	const std::vector<doublereal>& c,
	const std::vector<doublereal>& d,
	const int diff = 0);

/*
 * Evaluate multilinear function
*/
extern doublereal
leval(const doublereal& u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const int diff = 0);

#endif // INTERP_H

