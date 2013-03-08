# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2013
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
##        Copyright (C) 2011(-2013) all rights reserved.
##
##        The copyright of this code is transferred
##        to Pierangelo Masarati and Paolo Mantegazza
##        for use in the software MBDyn as described
##        in the GNU Public License version 2.1
##
##################################################################

function elem = MyElem1(pMbElem, pDM, HP)
    if ( ~HP.IsKeyWord("node1") )
        error("mboct: keyword node1 expected at line %s", HP.GetLineData());
    endif

    elem.pNode1 = pDM.ReadNode(HP, "STRUCTURAL");

    if ( ~HP.IsKeyWord("offset") )
        error("mboct: keyword offset expected at line %s", HP.GetLineData());
    endif

    elem.o1 = HP.GetPosRel(elem.pNode1);

    if ( ~HP.IsKeyWord("orientation") )
        error("mboct: keyword orientation expected at line %s", HP.GetLineData());
    endif

    elem.R1 = HP.GetRotRel(elem.pNode1);

    if ( ~HP.IsKeyWord("stiffnessx") )
        error("mboct: keyword stiffness x expected at line %s", HP.GetLineData());
    endif

    elem.sx = HP.GetReal();

    if ( ~HP.IsKeyWord("stiffnessy") )
        error("mboct: keyword stiffness y expected at line %s", HP.GetLineData());
    endif

    elem.sy = HP.GetReal();

    if ( ~HP.IsKeyWord("stiffnessz") )
        error("mboct: keyword stiffness z expected at line %s", HP.GetLineData());
    endif

    elem.sz = HP.GetReal();
    
    if ( ~HP.IsKeyWord("dampingcoefficient") )
        error("mboct: keyword damping coefficient expected at line %s", HP.GetLineData());
    endif

    elem.ks = HP.GetReal();

    if ( ~HP.IsKeyWord("forcevalue") )
        error("mboct: keyword force value expected at line %s", HP.GetLineData());
    endif

    elem.f1 = HP.GetDriveCaller();

    if ( ~HP.IsKeyWord("forcedirection") )
        error("mboct: keyword force direction expected at line %s", HP.GetLineData());
    endif

    elem.v1 = HP.GetPosRel(elem.pNode1);

    % Do some tests with the MBDyn parser ...
    if ( HP.IsKeyWord("stringvalue") && HP.IsArg() )
        fprintf(stderr, "mboct: string value = %s\n", HP.GetString(""));
    endif

    if ( HP.IsKeyWord("boolvalue") && HP.IsArg() )
        fprintf(stderr, "mboct: bool value = %d\n", HP.GetBool(true));
    endif

    if ( HP.IsKeyWord("intvalue") && HP.IsArg() )
        fprintf(stderr, "mboct: int value = %d\n", HP.GetInt(int32(0)));
    endif

    if ( HP.IsKeyWord("realvalue") && HP.IsArg() )
        fprintf(stderr, "mboct: real value = %g\n", HP.GetReal(0));
    endif

    if ( HP.IsKeyWord("stringwithdelimiter") && HP.IsStringWithDelims("SQUAREBRACKETS") )
        fprintf(stderr, "mboct: string with delimiter = %s\n", HP.GetStringWithDelims("SQUAREBRACKETS"));
    endif

    if (HP.IsKeyWord("drive1"))
        elem.drive1 = HP.GetDriveCaller();
    endif

    if (HP.IsKeyWord("drive2"))
        elem.drive2 = HP.GetDriveCaller();
    endif

%   elem.pNode1 = pDM.GetStructNode(int32(1));
%   elem.pNode1 = pDM.pFindNode("STRUCTURAL",int32(1));

    elem.pPrmNode1 = pDM.pFindNode("PARAMETER", int32(1));
    elem.V2_idx = elem.pNode1.iGetPrivDataIdx("XP[2]");

    elem.S = [ elem.sx, 0,       0;
               0,       elem.sy, 0;
               0,       0,       elem.sz ];

    elem.D = elem.ks * elem.S;

    elem.nIter = 0;
    elem.t = [pDM.dGetTime()];
    elem.X2 = [1];
    elem.XP2 = [-100 * elem.X2];
    elem.pDM = pDM;
    elem.pMbElem = pMbElem;

%   uncomment to see all elements ...
%    elem

    elem = class(elem,"MyElem1");

%   uncomment to call display ...
%    elem
endfunction
