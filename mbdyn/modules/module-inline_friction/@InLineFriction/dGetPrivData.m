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

% AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
%        Copyright (C) 2011(-2012) all rights reserved.
%
%        The copyright of this code is transferred
%        to Pierangelo Masarati and Paolo Mantegazza
%        for use in the software MBDyn as described
%        in the GNU Public License version 2.1


function dPrivData = dGetPrivData(elem, i)
    switch (i)
    case {1, 2}
        dPrivData = elem.lambda(i);
    case 3
        dPrivData = elem.z;
    case 4
        dPrivData = elem.zP;
    case 5
        X1 = elem.pNode1.GetXCurr();
        R1 = elem.pNode1.GetRCurr();
        X1P = elem.pNode1.GetVCurr();
        W1 = elem.pNode1.GetWCurr();

        X2 = elem.pNode2.GetXCurr();
        R2 = elem.pNode2.GetRCurr();
        X2P = elem.pNode2.GetVCurr();
        W2 = elem.pNode2.GetWCurr();

        DeltaXP = ComputeDeltaXP(elem, X1, R1, X1P, W1, X2, R2, X2P, W2);
        dPrivData = ComputeTau(elem, DeltaXP, norm(elem.lambda), elem.z);
    otherwise
        error("inline friction(%d): invalid index %d for private data", elem.pMbElem.GetLabel(), i);
    endswitch
endfunction
