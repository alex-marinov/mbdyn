# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2012
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
##        Copyright (C) 2011(-2012) all rights reserved.
##
##        The copyright of this code is transferred
##        to Pierangelo Masarati and Paolo Mantegazza
##        for use in the software MBDyn as described
##        in the GNU Public License version 2.1
##
##################################################################

function [f, ridx] = AssRes(elem, dCoef, XCurr, XPrimeCurr)
    iFirstIndex = elem.pMbElem.iGetFirstIndex() + int32(1);

%    label = elem.pNode1.GetLabel()

    ridx = [ elem.pNode1.iGetFirstMomentumIndex() + int32(1:3).';
             iFirstIndex ];

    X1 = elem.pNode1.GetXCurr();
    V1 = elem.pNode1.GetVCurr();

    X2 = XCurr(iFirstIndex);
    XP2 = XPrimeCurr(iFirstIndex);

    f = [ -elem.S * X1 - elem.D * V1;
            X2 + 0.01*XP2 ];


    V2 = elem.pNode1.dGetPrivData(elem.V2_idx);
    PrmX1 = elem.pPrmNode1.dGetX();

%    XCurr(1:end)
%    size(XCurr)
%    typeinfo(XCurr)

%{
    
    elem
    XCurr
    XPrimeCurr

    for i=1:XCurr.iGetSize()
        printf("XCurr(i)=%g\n", XCurr.dGetCoef(i));
    endfor
%}
endfunction
