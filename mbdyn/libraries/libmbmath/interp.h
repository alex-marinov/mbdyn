#include<vector>
#include "ac/f2c.h"
#include "myassert.h"
/*
 * Compute spline coefficients
*/
void spline (	
	const std::vector<doublereal>& x, 
	const std::vector<doublereal>& y, 
	std::vector<doublereal>& b, 
	std::vector<doublereal>& c, 
	std::vector<doublereal>& d);
/*
 * Evaluate spline
*/
double seval(const double u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const std::vector<doublereal>& b,
	const std::vector<doublereal>& c,
	const std::vector<doublereal>& d,
	const int diff = 0);
/*
 * Evaluate multilinear function
*/
double leval(const double u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const int diff = 0);

