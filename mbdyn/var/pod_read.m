%[W, cnt] = pod_read( filename [ , nstep ] )
%
%input:
%	filename	name of the input file
%	nstep		number of time steps (optional)
%
%output:
%	W		data (array if nstep is not defined, else a matrix
%			nstep x ndofs)
%	cnt		number of values read
%
function [W, cnt] = pod_read(filename, nstep)

[fid, msg] = fopen(filename, "r", "native");
if (fid == -1),
	error("pod_read: %s", msg)
end

if nargin < 2,
	[W, cnt] = fread(fid, Inf, "double", 0, "native");
else
	[W, cnt] = fread(fid, [nstep, Inf], "double", 0, "native");
end

fclose(fid);

