%[S, Aout, B, mn, scl, ee, vv, X, H, BB] = pod(A, ns, dt, dec)
%
%input:
%	A	data (number of frames x number of outputs)
%	ns	desired POMs
%	dt	time lag between two frames (optional, defaults to 1.0)
%       uu      input signals (number of frames x number of input signals; Matlab only) 
%	dec	decimation factor (Matlab only)
%
%output:
%	S	singular values
%	Aout	data after POD
%	B	reduced dataset POMs
%	mn	mean value of data
%	scl	scale factor of data minus mn
%	ee	eigenvalues of the transition matrix
%	vv	eigenvectors of the transition matrix
%	X	physical eigenvectors
%	H	transition matrix
%	BB	expanded POMs
%
function [S, Aout, B, mn, scl, ee, vv, X, H, BB] = pod(A, ns, dt, uu, dec)

% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2003
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
%
% Authors:	Giuseppe Quaranta	<quaranta@aero.polimi.it>
%		Pierangelo Masarati	<masarati@aero.polimi.it>
%

if (nargin < 3),
	dt = 1.;
end

[r, c] = size(A);

if (nargin < 2),
	ns = 0;
end

if ( ns > min(r,c)),
 	error('too many sv required');
end

% detrend and normalize
mn = mean(A);
A = A - ones(r, 1)*mn;
% scl = max(abs(A));
scl = std(A);
thr = 1.e-9;
lt = find(scl <= thr);
nlt = length(lt);
gt = find(scl > thr);
ngt = length(gt);
A = A(:, gt)./(ones(r, 1)*scl(gt));
for i = 1:nlt,
	disp(sprintf('dof %d: output is negligible', lt(i)));
end
disp(sprintf('using %d dofs (out of %d)', ngt, c));

if (ns == 0),
	S = [];
	Aout = [];
	B = [];
	return
end

if ((exist('dec') == 1)),
	if dec <= 1,
		error(sprintf('dec = %d is illegal', dec));
	end
	for i = 1:ngt,
		AA(:, i) = decimate(A(:, i), dec);
	end
        if (exist('uu') & ~isempty(uu)),
                [ru, cu] = size(uu);
                for i = 1:cu,
                        uuu(:, i) = decimate(uu(:, i), dec);
                end
                uu = uuu;
        end
	A = AA;
	dt = dt*dec;
	r = fix(r/dec);
end

[nt, nd] = size(A);
nn = min(nt, nd);
if (nn < ns),
	error(sprintf('number of requested modes %d is too high', ns));
end

if exist('OCTAVE_HOME'),
	%%% This is the big octave drawback: no eigs() ...
	if (nt > nd),
		[Btmp, Etmp] = eig(A'*A);
		Etmp = diag(Etmp);
		[Etmp2, I] = sort(Etmp);
		E = Etmp(I(ngt:-1:ngt-ns+1));
		B = Btmp(:, I(ngt:-1:ngt-ns+1))';
	else 
		[Utmp, Etmp] = eig(A*A');
		Etmp = diag(Etmp);
		[Etmp2, I] = sort(Etmp);
		E = Etmp(I(r:-1:r-ns+1));
		U = Utmp(:, I(r:-1:r-ns+1));
	end
else
	if (nt > nd),
		[B, E] = eigs(A'*A, ns);
		E = diag(E);
		B = B';
	else
		[U, E] = eigs(A*A', ns);
		E = diag(E);
	end
end

if (nt > nd),
	U = A*B';
	for i = 1:ns,
		s = U(:, i)'*U(:, i);
		S(i, 1) = sqrt(s);		% norm(B(i, :)), since S = sqrt(E)
		U(:, i) = U(:, i)/S(i, 1);
	end
else
	B = U'*A;
	for i = 1:ns,
		s = B(i, :)*B(i, :)';
		S(i, 1) = sqrt(s);		% norm(B(i, :)), since S = sqrt(E)
		B(i, :) = B(i, :)/S(i, 1);
	end
end

% Aout = B*A';			% = E^-1*U'*A*A' = E^-1*U'*U*E*U' = U'
Aout = U;

if exist('OCTAVE_HOME'),
	%%% This is a very rough estimate of the transition matrix ...
	H = (Aout(1:r-1, :)\Aout(2:r, :))';
else
       if (exist('uu') & ~isempty(uu)),
                [ru, cu] = size(uu);
                yu = iddata(Aout, uu);
		%% si potrebbe dare la struttura della B noto il nodo eccitato... 
                H = ssdata(arx(yu, [ones(ns), ones(ns, cu), zeros(ns, cu)], 'Covariance', 'None'));
    	else
		H = ssdata(arx(Aout, ones(ns), 'Covariance', 'None'));
	end
end

% physical eigenvalues and eigenvectors ...
[vv, eetmp] = eig(H);
ee = log(diag(eetmp))/dt;
B = diag(S)*B.*(ones(ns, 1)*scl(gt));
X = zeros(ns, c);
X(:, gt) = vv'*B;
BB = zeros(ns, c);
BB(:, gt) = B;

%
% node indices:
% awk '/struct node dofs: / {printf("l=[%s];\n",substr($0,19,length($0)))}' f.log
%
% eigenvectors (e.g. xy):
% k = 1.;
% x = mn(l+1)+k*[1 1]*X([6 7], l+1);
% y = mn(l+2)+k*[1 1]*X([6 7], l+2);
%
% eigenvalues:
% [x,I]=sort(imag(ee));I=I(length(I)/2+1:length(I));[I,-100*real(ee(I))./abs(ee(I)),imag(ee(I))/2/pi,(S'*abs(vv(:,I)))']

