%function [v, l] = eig3sym(A)
%
%	matrix A must be 3x3 symmetric
%
%	if nargout < 2, returns a vector containing the eigenvalues
%	otherwise, v contains the eigenvectors and diag(l) the eigenvalues
%
%From: `A robust algorithm for finding the eigenvalues and eigenvectors of
%3x3 symmetric matrices' W.M. Scherzinger, C.R. Dohrmann, Comput. Methods
%Appl. Mech. Engrg. 2008 doi:10.1016/j.cma.2008.03.031
%
%Copyright 2008-2017 Pierangelo Masarati <masarati@aero.polimi.it>

% $Header$
%
% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2017
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

function [v, l] = eig3sym(A)

if nargin ~= 1,
	error('invalid number of input parameters');
end

if nargout > 2,
	error('invalid number of output parameters');
end

[nr, nc] = size(A);
if (nr ~= 3 || nc ~= 3),
	error('`A'' must be 3x3');
end

if norm(A - A') > 1.e-15,
	error('`A'' must be symmetric')
end

v = zeros(3);
l = zeros(3, 1);

trA_3 = trace(A)/3;
AA = A - eye(3)*trA_3;
J2 = trace(AA*AA)/2;

if (abs(J2) < 1e-15),
	v = eye(3);
	l = diag(A);

else
	J3 = det(AA);

	dmy = J3/2*(3/J2)^(3/2);

	% NOTE: we want real eigenvalues; this requires the matrix to be
	% positive definite or semi-definite
	if (dmy < -1),
		dmy = -1;
		alpha = pi;
	elseif (dmy > 1),
		dmy = 1;
		alpha = 0;
	else
		alpha = acos(dmy)/3;
	end

	if (alpha < pi/6),
		idx1 = 1;
	else
		idx1 = 3;
	end

	eta_i = 2*sqrt(J2/3)*cos(alpha + 2/3*pi*(idx1 - 1));

	l(idx1) = eta_i + trA_3;

	% NOTE: there's a typo in the original paper;
	% AA must be used instead of A
	r = AA - eta_i*eye(3);
	nr = [norm(r(:,1)), norm(r(:,2)), norm(r(:,3))];
	imax = find(nr == max(nr));
	imax = imax(1);

	nrmax = nr(imax);
	nr(imax) = nr(1);
	nr(1) = nrmax;

	rmax = r(:, imax);
	r(:, imax) = r(:, 1);
	r(:, 1) = rmax;

	s = zeros(3, 2);
	s(:, 1) = r(:, 1)/nr(1);
	t2 = r(:, 2) - (s(:, 1)'*r(:, 2))*s(:, 1);
	nt2 = norm(t2);
	t3 = r(:, 3) - (s(:, 1)'*r(:, 3))*s(:, 1);
	nt3 = norm(t3);

	if (nt2 > nt3),
		s(:, 2) = t2/nt2;
	else
		s(:, 2) = t3/nt3;
	end

	v(:, idx1) = cross(s(:, 1), s(:, 2));

	AAA = zeros(3);
	AAA(1, 1) = eta_i;
	AAA(2:3, 2:3) = s'*AA*s;

	idx2 = fmod(idx1, 3) + 1;
	idx3 = fmod(idx1 + 1, 3) + 1;

	AAA22p33 = AAA(2, 2) + AAA(3, 3);
	AAA22m33 = AAA(2, 2) - AAA(3, 3);
	eta2 = (AAA22p33 - sign(AAA22m33)*sqrt(AAA22m33^2 + 4*AAA(2, 3)*AAA(3, 2)))/2;
	eta3 = AAA22p33 - eta2;

	l(idx2) = eta2 + trA_3;
	l(idx3) = eta3 + trA_3;

	u1 = (AA - eta2*eye(3))*s(:, 1);
	nu1 = norm(u1);
	u2 = (AA - eta2*eye(3))*s(:, 2);
	nu2 = norm(u2);

	if (nu1 > nu2),
		w1 = u1/nu1;
	else
		w1 = u2/nu2;
	end

	v(:, idx2) = cross(w1, v(:, idx1));
	v(:, idx3) = cross(v(:, idx2), v(:, idx1));
end

if nargout < 2,
	v = l;

else
	% NOTE: diag(l) makes little sense; it is done for consistency
	% with the output of eig().
	l = diag(l);
end

