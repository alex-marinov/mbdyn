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

% $Header$

function [W, cnt] = pod_read(filename)

% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2007
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
% Author:	Pierangelo Masarati <masarati@aero.polimi.it>
%

[fid, msg] = fopen(filename, 'r', 'native');
if (fid == -1),
	error('pod_read: %s', msg)
end

nrows = fread(fid, 1, 'ulong', 0, 'native');
ncols = fread(fid, 1, 'ulong', 0, 'native');

[W, cnt] = fread(fid, [nrows, Inf], 'double', 0, 'native');

fclose(fid);

W = W';

