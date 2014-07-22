# $Header$
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
#################################################################
##
## AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
##        Copyright (C) 2011(-2014) all rights reserved.
##
##        The copyright of this code is transferred
##        to Pierangelo Masarati and Paolo Mantegazza
##        for use in the software MBDyn as described
##        in the GNU Public License version 2.1
##
##################################################################

function [f, ridx] = InitialAssRes(elem, XCurr)

    ridx = elem.pMbElem.iGetFirstIndex() + int32(1:2).';
%{
    for i=1:length(ridx)
        X(i) = XCurr.dGetCoef(ridx(i));
    endfor
%}
    X = XCurr.GetVec(ridx);

    f = [ 500*X(1) + X(2) + 1; 
          300*X(2) + X(1) - 3];
%{
    elem
    XCurr
    pDM
    pElem
%}

%{
    for i=1:XCurr.iGetSize()
        printf("XCurr(i)=%g\n", XCurr.dGetCoef(i));
    endfor
%}
endfunction
