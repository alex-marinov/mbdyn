% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
%
% Copyright (C) 1996-2013
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

% AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
%        Copyright (C) 2011(-2013) all rights reserved.
%
%        The copyright of this code is transferred
%        to Pierangelo Masarati and Paolo Mantegazza
%        for use in the software MBDyn as described
%        in the GNU Public License version 2.1

function [X, XP] = GetStateVector(elem, XCurr, XPrimeCurr)
    X1 = elem.pNode1.GetXCurr();
    g1 = elem.pNode1.GetgCurr();
    X1P = elem.pNode1.GetVCurr();
    g1P = elem.pNode1.GetgPCurr();

    X2 = elem.pNode2.GetXCurr();
    g2 = elem.pNode2.GetgCurr();
    X2P = elem.pNode2.GetVCurr();
    g2P = elem.pNode2.GetgPCurr();

    iFirstIndex = elem.pMbElem.iGetFirstIndex();

    lambda = XCurr(iFirstIndex + int32(1).');
    lambdaP = XPrimeCurr(iFirstIndex + int32(1).');

    z = XCurr(iFirstIndex + int32(2:3).');
    zP = XPrimeCurr(iFirstIndex + int32(2:3).');

    X = [X1; g1; X2; g2; lambda; z];
    XP = [X1P; g1P; X2P; g2P; lambdaP; zP];
endfunction
