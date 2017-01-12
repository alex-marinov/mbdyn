% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
%
% Copyright (C) 1996-2017
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

% AUTHOR: Reinhard Resch <r.resch@secop.com>
%        Copyright (C) 2011(-2013) all rights reserved.
%
%        The copyright of this code is transferred
%        to Pierangelo Masarati and Paolo Mantegazza
%        for use in the software MBDyn as described
%        in the GNU Public License version 2.1

function elem = InLineFriction(pMbElem, pDM, HP)
    elem.pMbElem = pMbElem;
    elem.pDM = pDM;

    if ( ~HP.IsKeyWord("node1") )
        error("inline friction(%d): keyword node1 expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.pNode1 = pDM.ReadNode(HP, "STRUCTURAL");

    if ( HP.IsKeyWord("offset") )
        elem.o1 = HP.GetPosRel(elem.pNode1);
    else
        elem.o1 = zeros(3,1);
    endif

    if ( HP.IsKeyWord("orientation") )
        R1 = HP.GetRotRel(elem.pNode1);
    else
        R1 = eye(3);
    endif

    elem.e1 = R1(:,1);
    elem.e2 = R1(:,2);
    elem.e3 = R1(:,3);

    if ( ~HP.IsKeyWord("node2") ) 
        error("inline friction(%d): keyword node2 expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.pNode2 = pDM.ReadNode(HP, "STRUCTURAL");

    if ( HP.IsKeyWord("offset") )
        elem.o2 = HP.GetPosRel(elem.pNode2);
    else
        elem.o2 = zeros(3,1);
    endif

    if ( HP.IsKeyWord("coulombfrictioncoefficient") )
        elem.muc = HP.GetReal();
    else
        elem.muc = 0;
    endif

    if ( elem.muc < 0 )
        error("inline friction(%d): friction coefficient must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("staticfrictioncoefficient") )
        elem.mus = HP.GetReal();
    else
        elem.mus = 0;
    endif

    if ( elem.mus < 0 )
        error("inline friction(%d): friction coefficient must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("slidingvelocitycoefficient") )
        elem.vs = HP.GetReal();
    else
        elem.vs = 1;
    endif    

    if ( elem.vs < 0 )
        error("inline friction(%d): sliding velocity coefficient must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("slidingvelocityexponent") )
        elem.i = HP.GetReal();
    else
        elem.i = 1;
    endif

    if ( HP.IsKeyWord("microslipdisplacement") )
        elem.delta = HP.GetReal();
    else
        elem.delta = 1;
    endif

    if ( elem.delta <= 0 )
        error("inline friction(%d): micro slip displacement must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("initialstictionstate") )
        s0 = HP.GetReal();
        if ( abs(s0) > 1 )
            error("inline friction(%d): initial stiction state must be between -1 and 1 at line %s", pMbElem.GetLabel(), HP.GetLineData());
        endif
        elem.z = s0 * elem.delta;
    else
        elem.z = 0;
    endif

    if ( HP.IsKeyWord("initialstictionderivative") )
        elem.zP = HP.GetReal();
    else
        elem.zP = 0;
    endif

    if ( HP.IsKeyWord("viscousfrictioncoefficient") )
        elem.kv = HP.GetReal();
    else
        elem.kv = 0;
    endif

    if ( elem.kv < 0 )
        error("inline friction(%d): viscous friction coefficient must be greater then zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("stictionstateequationscale") )
        elem.PhiScale = HP.GetReal();
    else
        elem.PhiScale = 1;
    endif

    if ( elem.PhiScale == 0 )
        error("inline friction(%s): stiction state equation scale must not be equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.lambda = zeros(2, 1);

    elem = class(elem,"InLineFriction");
endfunction
