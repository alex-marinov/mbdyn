%[W, cnt] = pod_read( filename [ , ndof ] )
%
%input:
%	filename	name of the input file
%	ndof		number of dofs (optional)
%
%output:
%	W		data (array if ndof is not defined, else a matrix
%			nstep x ndofs)
%	cnt		number of values read
%
function [W, cnt] = pod_read(filename, ndof)

[fid, msg] = fopen(filename, 'r', 'native');
if (fid == -1),
	error('pod_read: %s', msg)
end

if 0
	nrows = fread(fid, 1, 'ulong', 0, 'native');
	ncols = fread(fid, 1, 'ulong', 0, 'native');
else
	ncols = ndof;
end

if nargin < 2,
	[W, cnt] = fread(fid, Inf, 'double', 0, 'native');
else
	[W, cnt] = fread(fid, [ncols, Inf], 'double', 0, 'native');
end

fclose(fid);

W = W';

