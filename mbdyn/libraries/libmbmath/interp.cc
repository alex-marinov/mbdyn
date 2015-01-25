/* $Header$ */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "interp.h"

void
spline(const std::vector<doublereal>& x, 
	const std::vector<doublereal>& y, 
	std::vector<doublereal>& b, 
	std::vector<doublereal>& c, 
	std::vector<doublereal>& d)
{
/*
*
*  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
*  for a cubic interpolating spline
*
*    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
*
*    for  x(i) .le. x .le. x(i+1)
*
*  input..
*
*    n = the number of data points or knots (n.ge.2)
*    x = the abscissas of the knots in strictly increasing order
*    y = the ordinates of the knots
*
*  output..
*
*    b, c, d  = arrays of spline coefficients as defined above.
*
*  using  p  to denote differentiation,
*
*    y(i) = s(x(i))
*    b(i) = sp(x(i))
*    c(i) = spp(x(i))/2
*    d(i) = sppp(x(i))/6  (derivative from the right)
*
*  the accompanying function subprogram  seval  can be used
*  to evaluate the spline.
*
*
*/
	int n, nm1;
	ASSERTMSGBREAK(x.size() == y.size(), "Error in spline, x.size()!=y.size()");

	n = x.size();
	b.resize(n);
	c.resize(n);
	d.resize(n);

	nm1 = n - 1;
	if (n < 2) {
		return;
	}
	if (n < 3) {
		b[0] = (y[1] - y[0])/(x[1] - x[0]);
		c[0] = 0.;
		d[0] = 0.;
		b[1] = b[0];
		c[1] = 0.;
		d[1] = 0.;
	} else {
/*
*  set up tridiagonal system
*
*  b = diagonal, d = offdiagonal, c = right hand side.
*/
		d[0] = x[1] - x[0];
		c[1] = (y[1] - y[0])/d[0];
		for (int i = 1; i < nm1; i++) {
			d[i] = x[i + 1] - x[i];
			b[i] = 2.*(d[i - 1] + d[i]);
			c[i + 1] = (y[i + 1] - y[i])/d[i];
			c[i] = c[i + 1] - c[i];
		}
/*
*  end conditions.  third derivatives at  x(1)  and  x(n)
*  obtained from divided differences
*/
		b[0] = -d[0];
		b[n - 1] = -d[n - 2];
		c[0] = 0.;
		c[n - 1] = 0.;
		if (n != 3) {
			c[0] = c[2]/(x[3] - x[1]) - c[1]/(x[2] - x[0]);
			c[n - 1] = c[n - 2]/(x[n - 1] - x[n - 3]) - c[n - 3]/(x[n - 2] - x[n - 4]);
			c[0] = c[0]*d[0]*d[0]/(x[3] - x[0]);
			c[n - 1] = -c[n - 1]*d[n - 2]*d[n - 2]/(x[n - 1] - x[n - 4]);
		}
/*
*  forward elimination
*/
		for (int i = 1; i < n; i++) {
		         doublereal t = d[i - 1]/b[i - 1];
			 b[i] = b[i] - t*d[i - 1];
			 c[i] = c[i] - t*c[i - 1];
		}

/*
*  back substitution
*/
		c[n - 1] = c[n - 1]/b[n - 1];
		for (int ib = 0; ib < nm1; ib++) {
			int i = n - ib - 2;
			c[i] = (c[i] - d[i]*c[i + 1])/b[i];
		}
/*
*  c(i) is now the sigma(i) of the text
*
*  compute polynomial coefficients
*/
		b[n-1] = (y[n - 1] - y[nm1 - 1])/d[nm1 - 1] + d[nm1 - 1]*(c[nm1 - 1] + 2.*c[n - 1]);
		for (int i=0; i<nm1; i++) {
			b[i] = (y[i + 1] - y[i])/d[i] - d[i]*(c[i + 1] + 2.*c[i]);
			d[i] = (c[i + 1] - c[i])/d[i];
			c[i] = 3.*c[i];
		}
		c[n - 1] = 3.*c[n - 1];
		d[n - 1] = d[n - 2];
	}

	return;
}

doublereal
seval(const doublereal& u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const std::vector<doublereal>& b,
	const std::vector<doublereal>& c,
	const std::vector<doublereal>& d,
	const int diff)
{
/*
*  this subroutine evaluates the cubic spline function
*
*    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
*
*    where  x(i) .lt. u .lt. x(i+1), using horner's rule
*
*  if  u .lt. x(1) then  i = 1  is used.
*  if  u .ge. x(n) then  i = n  is used.
*
*  input..
*
*    n = the number of data points
*    u = the abscissa at which the spline is to be evaluated
*    x,y = the arrays of data abscissas and ordinates
*    b,c,d = arrays of spline coefficients computed by spline
*/
	int n;
	doublereal dx;

	ASSERTMSGBREAK(x.size() == y.size(), "Error in seval, x.size()!=y.size()");
	ASSERTMSGBREAK(x.size() == b.size(), "Error in seval, x.size()!=b.size()");
	ASSERTMSGBREAK(x.size() == c.size(), "Error in seval, x.size()!=c.size()");
	ASSERTMSGBREAK(x.size() == d.size(), "Error in seval, x.size()!=d.size()");

	n = x.size();


/*
* outside the interval
*/

	if (u <= x[0]) {
		dx = u - x[0];
		switch (diff) {
		case 0:
			return y[0] + dx*b[0];
			break;
		case 1:
			return b[0];
			break;
		default:
			return 0.;
			break;
		}

	} else if (u >= x[n - 1]) {
		dx = u - x[n - 1];
		switch (diff) {
		case 0:
			return y[n - 1] + dx*b[n - 1];
			break;
		case 1:
			return b[n - 1];
			break;
		default:
			return 0.;
			break;
		}

	} else {
		int i, j;
/*
*  binary search
*/

		i = 0;
		j = n;
		do {
			int k = (i + j)/2;
			if (u < x[k]) {
				j = k;
			} else {
				i = k;
			}
		} while (j > i + 1);
/*
*  evaluate spline
*/
		dx = u - x[i];
	
		ASSERTMSGBREAK(diff >=0, "Error in spline evaluation, diff<0");
		
		switch (diff) {
		case 0:
			return y[i] + dx*(b[i] + dx*(c[i] + dx*d[i]));
			break;
		case 1:
			return b[i] + dx*(2.*c[i] + 3.*dx*d[i]);
			break;
		case 2:
			return 2.*c[i] + 6.*dx*d[i];
			break;
		case 3:
			return 6.*d[i];
			break;
		default:
			return 0.;
			break;
		}
	}

	return 0.;
}

doublereal
leval(const doublereal& u,
	const std::vector<doublereal>& x,
	const std::vector<doublereal>& y,
	const int diff)
{
/*
*  this subroutine evaluates the multilinear funcion
*
*    seval = y(i) + (y(i+1)-y(i))(x(i+1)-x(i))*(u-x(i))
*
*    where  x(i) .lt. u .lt. x(i+1), using horner's rule
*
*  if  u .lt. x(1) then  i = 1  is used.
*  if  u .ge. x(n) then  i = n  is used.
*
*  input..
*
*    n = the number of data points
*    u = the abscissa at which the spline is to be evaluated
*    x,y = the arrays of data abscissas and ordinates
*    b,c,d = arrays of spline coefficients computed by spline
*/
	int i, j, n;
	doublereal dx, m;

	ASSERTMSGBREAK(x.size() == y.size(), "Error in leval, x.size()!=y.size()");

/*
*  binary search
*/

	n = x.size();
	i = 0;
	j = n;
	do {
		int k = (i + j)/2;
		if (u < x[k]) {
			j = k;
		} else {
			i = k;
		}
	} while (j > i + 1);
/*
*  evaluate 
*/
	if (i == n-1) {
		i--;
	}
	m = (y[i+1] - y[i])/(x[i+1] - x[i]);
	dx = u - x[i];
	
	
	ASSERTMSGBREAK(diff >=0, "Error in multilinear evaluation, diff<0");
	
	switch (diff) {
	case 0:
		return y[i] + m*dx;
		break;
	case 1:
		return m;
		break;
	default:
		return 0.;
		break;
	}
	return 0.;
}

