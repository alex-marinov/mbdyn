%Input:
%  mbeig('name');
%  mbeig(Aminus, Aplus, dCoef);
%  mbeig(Aminus, Aplus, dCoef, vv);
%  mbeig(Aminus, Aplus, dCoef, fmin, fmax);
%
%  'name': name of file output by MBDyn (using 'output sparse matrices, output geometry')
%  Aminus, Aplus, dCoef: matrices and coefficients as output by MBDyn in Matlab-compatible file
%  vv: vector of indexes of modes to be used to reduce matrices
%  fmin, fmax: minimum and maximum angular frequency (imaginary part of eigenvalues) to be used to reduce matrices
%
%Output:
%  l = mbeig()
%  [y, l] = mbeig();
%  [y, yc, l] = mbeig();
%  [y, yc, l, A, E] = mbeig();
%
%  l: eigenvalues (vector if only output; on diagonal of diagonal matrix otherwise, as in Matlab's eig())
%  y: "right" eigenvectors
%  yc: "left" eigenvectors
%  A, E: descriptor form matrices (E x' = A x; reduced if either 'vv' or 'fmin', 'fmax' are provided)

function [out1, out2, out3, out4, out5] = mbeig(in1, in2, in3, in4, in5);

% NOTE: outdated
%
% If "output matrices" is set, eigenvalues and eigenvectors are computed
% by octave/matlab, and the reduced model is generated.
%
% If "ouput eigenvectors" is set, eigenvalues and eigenvectors computed
% by MBDyn are used; only if "output matrices" is set, the reduced model
% is generated.
%
% NOTE: you may set "upper frequency limit", but you should not
% set "lower frequency limit"

% mbeig('name');
% mbeig(Aminus, Aplus, dCoef);
% mbeig(Aminus, Aplus, dCoef, vv);
% mbeig(Aminus, Aplus, dCoef, fmin, fmax);

% l = mbeig();
% [y, l] = mbeig();
% [y, yc, l] = mbeig();
% [y, yc, l, A, E] = mbeig();

if (nargin == 1),
	name = in1;
	eval(sprintf('%s;', name));

elseif (nargin >= 3),
	if (nargin > 5),
		error(sprintf('mbeig: invalid input arg number (%d)', nargin));

	elseif (nargin == 5),
		fmin = in4;
		fmax = in5;

	elseif (nargin == 4),
		vv = in4;
	end

	Aminus = in1;
	Aplus = in2;
	dCoef = in3;

else
	error(sprintf('mbeig: invalid input arg number (%d)', nargin));
end

if ((nargout == 0) || (nargout == 1) || (nargout == 2) || (nargout == 3) || (nargout == 5)),
	% ...

else
	error(sprintf('mbeig: invalid output arg number (%d)', nargout));
end

if exist('VR') == 1 & exist('alpha') == 1,
	if (isempty(find(abs(alpha(:, 3)) < 1e-12)) ~= 1),
		error('Houston, we''ve had a problem: (alpha_r + j*alpha_i)/beta: some beta are too small');
	end
	Y = VR;
	L = (alpha(:, 1) + sqrt(-1)*alpha(:, 2))./alpha(:, 3);

	if (exist('VL') == 1),
		Yc = VL;
	end

else
	% "right" and "left" eigenproblem
	if exist('OCTAVE_HOME'),
		[dmy, dmy, dmy, dmy, Y, Yc, L] = qz(Aminus, Aplus);
		[n, dmy] = size(Aminus);
		nL = length(L);
		if (nL < n),
			error(sprintf('Houston, we''ve had a problem: too few valid eigenvalues (%d instead of %d)...', nL, n));
		end
	else
		[AA, BB, dmy, dmy, Y, Yc] = qz(Aminus, Aplus);
		% NOTE: by default, matlab performs the complex QZ transform.
		% As a consequence, the eigenvalues are diag(AA)./diag(BB).
		% This behavior is indicated as legacy; it might change
		% in the future.
		L = diag(AA)./diag(BB);
	end
end

% reduce "right" problem
l = (L - 1)./(L + 1)/dCoef;
if (nargout <= 1),
	out1 = l;
	return;

elseif (nargout == 2),
	% for compatibility with matlab's eig()...
	out2 = diag(l);

elseif (nargout >= 3),
	% for compatibility with matlab's eig()...
	out3 = diag(l);
end


if (exist('vv') == 1),
	v = vv;
elseif (exist('fmin')),
	v = find((abs(imag(l)) >= fmin) & (abs(imag(l)) <= fmax));
else
	v = [1:size(Aminus, 1)];
end

m = length(v);

% combine complex modes; leave real alone
T = zeros(m);
i = 1;
while (i <= m),
	if (norm(imag(Y(:, v(i)))) == 0.),
		T(i, i) = 1;
		i = i + 1;
	else
		T(i:i + 1, i:i + 1) = [.5 -.5*sqrt(-1); .5 .5*sqrt(-1)];
		i = i + 2;
	end
end
y = Y(:, v)*T;
out1 = y;

if (nargout == 2),
	return;
end

% matrices Aplus, Aminus and dCoef from MBDyn's run must be available
if exist('Aplus') ~= 1 | exist('Aminus') ~= 1 | exist('dCoef') ~= 1,
	error('need (at least) Aplus, Aminus and dCoef; execute ''<outfile>.m'' first');
end

Tc = zeros(m);
i = 1;
while (i <= m),
	if (norm(imag(Yc(:, v(i)))) == 0.),
		Tc(i, i) = 1;
		i = i + 1;
	else
		Tc(i:i + 1, i:i + 1) = [.5 -.5*sqrt(-1); .5 .5*sqrt(-1)];
		i = i + 2;
	end
end
yc = Yc(:, v)*Tc;
out2 = yc;

if (nargout == 3),
	return;
end

% matrices
E = .5*(yc'*Aplus*y + yc'*Aminus*y);
A = -.5/dCoef*(yc'*Aplus*y - yc'*Aminus*y);

out4 = A;
out5 = E;

% TODO...
if exist('b') == 1,
	B = yc'*b;
end
if exist('c') == 1,
	C = c*y;
end
