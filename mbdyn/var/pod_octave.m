%[S, Aout, B, mn, scl, ee, vv, X] = pod_octave(A, ns, dt)
%
%input:
%	A	data (number of frames x number of outputs)
%	ns	desired POMs
%	dt	time lag between two frames (optional, defaults to 1.0)
%
%output:
%	S	singular values
%	Aout	data after POD
%	B	POMs
%	mn	mean value of data
%	scl	scale factor of data-mn
%	ee	eigenvalues
%	vv	eigenvectors
%	X	physical eigenvectors
%
function [S, Aout, B, mn, scl, ee, vv, X] = pod_octave(A, ns, dt)

if nargin < 3,
	dt = 1.;
end

[r, c] = size(A);

if nargin < 2,
	ns = c;
else 
	if ( ns > c),
  	  error('too many sv required');
	end
end

% detrend and normalize
for i = 1:c,
    mn(i) = mean(A(:, i));
    A(:, i) = A(:, i)-mn(i);
    scl(i) = max(abs(A(:, i)));
    if (scl(i) > 1.e-9),
        A(:, i) = A(:, i)/scl(i);
    else
	A(:, i) = zeros(r, 1);
	disp(sprintf("dof %d: output is negligible", i));
    end
end

%%% This is the big octave drawback: no eigs() ...
[Utmp, Etmp] = eig(A*A');
Etmp = diag(Etmp);
[Etmp2, I] = sort(Etmp);
E = Etmp(I(r:-1:r-ns+1));
U = Utmp(:, I(r:-1:r-ns+1));

B = U'*A;
for i = 1:ns,
    S(i) = norm(B(i, :));	% S = sqrt(E)
    B(i, :) = B(i, :)/S(i)^2;
end
% Aout = B*A';			% = E^-1*U'*A*A' = E^-1*U'*U*E*U' = U'
Aout = U';

%%% This is a very rough estimate of the transition matrix ...
H = ((Aout(:, 1:r-1)')\(Aout(:, 2:r)'))';

% physical eigenvalues and eigenvectors ...
[vv, eetmp] = eig(H);
ee = diag(eetmp)/dt;
X = vv'*B;

% eigenvectors:
% plot(([1 1;sqrt(-1) -sqrt(-1)]*X([6 7], l+1))')
