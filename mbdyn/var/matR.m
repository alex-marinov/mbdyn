% $Header$
%
% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
% 
% Copyright (C) 1996-2011
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

function R = matR(theta)

t = norm(theta);
if (t > 0.),
	a = sin(t)/t;
	b = (1 - cos(t))/t^2;
	R = eye(3) + a*cross(theta) + b*cross(theta)*cross(theta);
else
	R = eye(3);
end
