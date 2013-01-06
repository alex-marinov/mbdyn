%[S, Aout, B, mn, scl, ee, vv, X, H, BB, fm] = pod(A, ns, dt, uu, dec, ord, thr, wgt)
%
%input:
%	A	data (n. of frames x n. of outputs)
%	ns	desired POMs (ns(1): POM, ns(2): identification)
%	dt	time lag between two frames (optional, defaults to 1.0)
%	uu      input signals (n. of frames x n. of inputs; MATLAB only) 
%	dec	decimation factor
%	ord	model order (default 1)
%	thr	threshold for discarding signals
%	wgt	'norm', 'none' or weight matrix
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
%	fm      figure of merit: evaluates the relative quality of POMs

% $Header$

function [S, Aout, B, mn, scl, ee, vv, X, H, BB, fm] = pod(A, ns, dt, uu, dec, ord, thr, wgt)

% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2013
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

if (nargin < 8 || isempty(wgt)),
	% backwards compatibility
	wgt = 'norm';
end

if (nargin < 7 || isempty(thr)),
	thr = 1.e-16;
end

if (nargin < 6),
	ord = 1;
end

if (nargin < 5),
	dec = 0;
end

if (nargin < 3),
	dt = 1.;
end

if (nargin < 2),
	ns = 0;
end

[r, c] = size(A);

% ns(1): POM; ns(2): identification
if (length(ns) == 1),
	ns2 = ns;

elseif (length(ns) >= 2),
	ns2 = ns(2);
	ns = ns(1);

	if (ns2 > ns),
		error(sprintf('illegal ns(2) = %d > ns(1) = %d', ns2, ns));
	end
end

if ( ns > min(r, c)),
 	error('too many sv required');
end

if (ord < 1),
	error(sprintf('illegal order %d', ord));
end

% detrend and normalize
mn = mean(A);
A = A - ones(r, 1)*mn;
if (strcmp(wgt, 'norm')),
	scl = std(A);
	le = find(scl <= thr);
	nle = length(le);
	gt = find(scl > thr);
	ngt = length(gt);
	A = A(:, gt)./(ones(r, 1)*scl(gt));
	for i = 1:nle,
		disp(sprintf('dof %d: output is negligible', le(i)));
	end
	msg = sprintf('using %d dofs', ngt);
	if (ngt < c),
		msg = [msg, sprintf(' (out of %d)', c)];
	end
	msg = [msg, sprintf(' and %d', r)];
	if (dec > 1),
		msg = [msg, sprintf('/%d', dec)];
	end
	msg = [msg, ' samples'];
	disp(msg);
elseif (length(wgt) == 4 && wgt == 'none'),
	[dmy, ngt] = size(A);
	scl = ones(1, ngt);
	gt = [1:ngt];
	le = [];
	nle = 0;
else
	[dmy, ngt] = size(A);
	nwgt = size(wgt);
	if (nwgt(1) ~= 1 || nwgt(2) ~= ngt),
		error(sprintf('illegal dimensions %dx%d of vector wgt', nwgt(1), nwgt(2)));
	end
	scl = wgt;
	gt = [1:ngt];
	le = [];
	nle = 0;
end

if (ns == 0),
	S = [];
	Aout = [];
	B = [];
	return
end

if (dec > 1),
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
	error(sprintf('number of requested modes %d is too high (%d max)', ns, nn));
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

MATLAB_IS_BRAINDEAD = 1;
if MATLAB_IS_BRAINDEAD, % exist('OCTAVE_HOME'),
	%%% This is a very rough estimate of the transition matrix ...
	%%% H = (Aout(1:r-1, :)\Aout(2:r, :))';

	hh = [];
	for i = 1:ord,
		hh(:, ns2*(i - 1) + 1:ns2*i) = Aout(ord - i + 1:r - i, 1:ns2);
	end
	H = (hh\Aout(ord + 1:r, 1:ns2))';
	if (ord > 1),
		H(ns2 + 1:ord*ns2, 1:(ord - 1)*ns2) = eye((ord - 1)*ns2);
		C = [eye(ns2), zeros(ns2, (ord - 1)*ns2)];
	end
else
       if (exist('uu') & ~isempty(uu)),
                [ru, cu] = size(uu);
                yu = iddata(Aout(:, 1:ns2), uu);
		%% si potrebbe dare la struttura della B noto il nodo eccitato... 
                [H, Btmp, C] = ssdata(arx(yu, [ord*ones(ns2), ones(ns2, cu), zeros(ns, cu)], 'Covariance', 'None'));
    	else
		[H, Btmp, C] = ssdata(arx(Aout, ord*ones(ns2), 'Covariance', 'None'));
	end
end

% physical eigenvalues and eigenvectors ...
[vv, eetmp] = eig(H);
ee = log(diag(eetmp))/dt;
B = diag(S(1:ns2))*B(1:ns2, :).*(ones(ns2, 1)*scl(gt));
if (ord > 1),
	vv = C*vv;
end
X = zeros(ord*ns2, c);
X(:, gt) = vv.'*B(1:ns2, :);
BB = zeros(ns2, c);
BB(:, gt) = B;

% figure of merit...
fm = (S(1:ns2)'*abs(vv))';
fm = fm./max(fm);

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

