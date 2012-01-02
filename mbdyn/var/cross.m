%function out = cross(v1, v2)
%
%if only v1 is given, out is the cross-product matrix of v1, `v1 times'
%if v1 and v2 are given, out is the cross-product `v1 times v2'
%
%Copyright 2008-2012 Pierangelo Masarati <masarati@aero.polimi.it>

% $Header$
%
% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2012
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

function out = cross(v1, v2);
if nargin == 1,
	out = [0 -v1(3) v1(2); v1(3) 0 -v1(1); -v1(2) v1(1) 0];
elseif nargin == 2,
	out = [v1(2)*v2(3) - v1(3)*v2(2); v1(3)*v2(1) - v1(1)*v2(3); v1(1)*v2(2) - v1(2)*v2(1)];
else
	error('args');
end
