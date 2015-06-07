## -*- texinfo -*-
## @deftypefn {Function File} {[@var{y}, @var{J}] =} mbdyn_derivative(@var{F}, @var{x}, @var{varargin})
## Evaluate @var{F} for a given input @var{x} and compute the jacobian, such that
## @example
##
##            d
## @var{J}(i,j) = ----- @var{y}(i) where @var{y} = @var{F}(@var{x}, @var{varargin}@{:@})
##          d@var{x}(j)
## 
## @end example
##
## If @var{x} is complex, the above holds for the directional derivatives 
## along the real axis
##
## Derivatives are computed analytically via Automatic Differentiation
## @end deftypefn
## @seealso{use_sparse_jacobians}
#################################################################################

## $Header$

#################################################################################
## Copyright (C) 2006, 2007 Thomas Kasper, <thomaskasper@gmx.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.
#################################################################################

#################################################################################
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2014
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
# 
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
# 
# Changing this copyright notice is forbidden.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 
#################################################################################

#################################################################################
##
## AUTHOR: Reinhard Resch <r.resch@secop.com>
##        Copyright (C) 2011(-2013) all rights reserved.
##
##        The copyright of this code is transferred
##        to Pierangelo Masarati and Paolo Mantegazza
##        for use in the software MBDyn as described
##        in the GNU Public License version 2.1
##
#################################################################################

## This function is part of the MBDyn - octave interface.
## It is needed for automatic forward differentiation of octave functions
## called from a octave drive caller in MBDyn

function [y, J] = mbdyn_derivative(F, x, varargin)
  if nargin > 1
    if (3 == exist("gradinit", "file"))
      f = feval(F, gradinit(x), varargin{:});
   
      if isgradient (f)
        y = f.x;
        J = f.J;
      else
        warning ("AD: function not differentiable or not dependent on input")
        y = f;
        m = numel (y);
        n = numel (x);
        if use_sparse_jacobians () != 0
          J = sparse (m, n);
        else
          J = zeros (m, n);
        endif
      endif

      # FIXME: reshape needed for TplDriveCaller<Mat3x3>
      if ( isscalar(x) && ismatrix(y) )
        J = reshape(J, rows(y), columns(y));
      endif
    else
      y = feval(F, x, varargin{:});

      if (~(isscalar(x) && ismatrix(y)))
        J = zeros(length(y), length(x));
      endif

      for i=1:length(x)
        dx = zeros(size(x));
        dx(i) = sqrt(eps) * norm(x) + sqrt(eps);
        dy_dx = (feval(F, x + dx, varargin{:}) - y)  / dx(i);

        # FIXME: reshape needed for TplDriveCaller<Mat3x3>
        if ( isscalar(x) && ismatrix(y) )
          J = dy_dx;
        else
          J(:, i) = dy_dx;
        endif
      endfor
    endif
  else 
    usage ("[y, J] = mbdyn_derivative(F, x, varargin)");
  endif
endfunction

%!function y=f2(x)
%! y(1,1) = (x - 1) / (x + 1);
%! y(1,2) = 1/(1+x)^2;
%! y(2,1) = (x - 1) * (x + 1);
%! y(2,2) = x + 1 / x;

%!function y=f1(x)
%!  y(1)=sin(x(1)) * tan(x(2)) + exp(x(3));
%!  y(2)=cos(x(2)) + sin(x(3))^2 * x(1);

%!test
%! x = [0.57; 0.68; -1.3];
%! [y,J]=mbdyn_derivative(@f1, x)

%!test
%! x = 20;
%! [y, J]=mbdyn_derivative(@f2, x)

