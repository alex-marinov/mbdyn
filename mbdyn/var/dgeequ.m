%function [R, C, ROWCND, COLCND, AMAX] = dgeequ(A)
%
% computes diagonal scaling matrices that modify matrix A
% such that the magnitude of the largest coefficient in each
% row/column is 1.
%
% usage: given a matrix A, intended for the solution of the
% linear problem
%
%	x = A\b
%
% the solution with scaling is given by
%
%	[R, C] = dgeequ(A);
%	x = diag(C) * ( ( diag(R) * A * diag(C) ) \ (diag(R) * b) )
%
% Based on lapack's dgeequ routine.
%
%Copyright 2008-2014 Pierangelo Masarati <masarati@aero.polimi.it>

% $Header$
%
% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2014
% 
% Pierangelo Masarati	<masarati@aero.polimi.it>
% Paolo Mantegazza	<mantegazza@aero.polimi.it>
% 
% Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
% via La Masa, 34 - 20156 Milano, Italy
% http://www.aero.polimi.it
% 
% Changing this copyright notice is forbidden.
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (version 2 of the License).
% 
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 

function [R, C, ROWCND, COLCND, AMAX] = dgeequ(A);

SMLNUM = 1e-15;
BIGNUM = 1./SMLNUM;

[N, M] = size(A);

X = abs(A);
R = max(X')';

RCMIN = min(min(R), BIGNUM);
RCMAX = max(max(R), 0.);

AMAX = RCMAX;

if (RCMIN == 0.),
	error(sprintf('R(%d) == 0', find(R == 0.)));
end

for I = 1:N,
	R(I) = 1./min(max(R(I), SMLNUM), BIGNUM);
end

ROWCND = max(RCMIN, SMLNUM)/min(RCMAX, BIGNUM);

X = diag(R)*X;
C = max(X)';

RCMIN = min(min(C), 1/1e-16);
RCMAX = max(max(C), 0.);

if (RCMIN == 0.),
	error(sprintf('C(%d) == 0', find(R == 0.)));
end

for I = 1:M,
	C(I) = 1./min(max(C(I), SMLNUM), BIGNUM);
end

COLCND = max(RCMIN, SMLNUM)/min(RCMAX, BIGNUM);

