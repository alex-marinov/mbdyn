%[W, cnt] = pod_read( filename [ , ndof ] )
%
%input:
%	filename	name of the input file; contains the output
%			in binary format, preceded by the matrix
%			dimensions:
%
%				ulong	1		nrows
%				ulong	1		ncols
%				double	nrows*ncols	output
%
%output:
%	W		data (array if ndof is not defined, else a matrix
%			nstep x ndofs)
%	cnt		number of values read
%
function [W, cnt] = pod_read(filename)

[fid, msg] = fopen(filename, 'r', 'native');
if (fid == -1),
	error('pod_read: %s', msg)
end

nrows = fread(fid, 1, 'ulong', 0, 'native');
ncols = fread(fid, 1, 'ulong', 0, 'native');

if nargin < 2,
	[W, cnt] = fread(fid, Inf, 'double', 0, 'native');
else
	[W, cnt] = fread(fid, [nrows, Inf], 'double', 0, 'native');
end

fclose(fid);

W = W';

