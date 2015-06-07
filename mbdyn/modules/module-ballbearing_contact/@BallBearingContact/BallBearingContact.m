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

% AUTHOR: Reinhard Resch <r.resch@secop.com>
%        Copyright (C) 2011(-2013) all rights reserved.
%
%        The copyright of this code is transferred
%        to Pierangelo Masarati and Paolo Mantegazza
%        for use in the software MBDyn as described
%        in the GNU Public License version 2.1

function elem = BallBearingContact(pMbElem, pDM, HP)
    elem.pMbElem = pMbElem;
    elem.pDM = pDM;

	if ( ~HP.IsKeyWord("ballradius") )
		error("ball bearing contact(%d): keyword \"ball radius\" expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
	endif

	elem.R = HP.GetReal();

	if ( HP.IsKeyWord("elasticmodulus") )
		elem.E = HP.GetReal();
	else
		if ( ~HP.IsKeyWord("elasticmodulusball") )
			error("ball bearing contact(%d): keyword \"elastic modulus ball\" expected at line %s", pMbElem.GetLabel(), HP.GetLineData());	
		endif

		E1 = HP.GetReal();

		if ( ~HP.IsKeyWord("poissonratioball") )
			error("ball bearing contact(%d): keyword \"poisson ratio ball\" expected at line %s", pMbElem.GetLabel(), HP.GetLineData());				
		endif

		nu1 = HP.GetReal();

		if ( ~HP.IsKeyWord("elasticmoduluswasherdisk") )
			error("ball bearing contact(%d): keyword \"elastic modulus washer disk\" expected at line %s", pMbElem.GetLabel(), HP.GetLineData());				
		endif

		E2 = HP.GetReal();

		if ( ~HP.IsKeyWord("poissonratiowasherdisk") )
			error("ball bearing contact(%d): keyword \"poisson ratio washer disk\" expected at line %s", pMbElem.GetLabel(), HP.GetLineData());				
		endif

		nu2 = HP.GetReal();

		elem.E = 1 / ((1 - nu1^2) / E1 + (1 - nu2^2) / E2);
	endif		

    if ( ~HP.IsKeyWord("node1") )
        error("ball bearing contact(%d): keyword node1 expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.pNode1 = pDM.ReadNode(HP, "STRUCTURAL");

    if ( ~HP.IsKeyWord("node2") ) 
        error("ball bearing contact(%d): keyword node2 expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.pNode2 = pDM.ReadNode(HP, "STRUCTURAL");

    if ( HP.IsKeyWord("offset") )
        elem.o2 = HP.GetPosRel(elem.pNode2);
    else
        elem.o2 = zeros(3,1);
    endif

    if ( HP.IsKeyWord("orientation") )
        elem.Rt2 = HP.GetRotRel(elem.pNode2);
    else
        elem.Rt2 = eye(3);
    endif

    if ( ~(HP.IsKeyWord("coulombfrictioncoefficient") || HP.IsKeyWord("coulombfrictioncoefficientx")) )
        error("ball bearing contact(%d): keyword coulomb friction coefficient x expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif
 
    mukx = HP.GetReal();

    if ( mukx <= 0 )
        error("ball bearing contact(%d): coulomb friction coefficient x must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("coulombfrictioncoefficienty") )
        muky = HP.GetReal();
    else
        muky = mukx;
    endif

    if ( muky <= 0 )
        error("ball bearing contact(%d): coulomb friction coefficient y must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.Mk = diag([mukx, muky]);

    if ( ~(HP.IsKeyWord("staticfrictioncoefficient") || HP.IsKeyWord("staticfrictioncoefficientx")) )
        error("ball bearing contact(%d): keyword static friction coefficient x expected at line %d", pMbElem.GetLabel(), HP.GetLineData());
    endif

    musx = HP.GetReal();

    if ( musx <= 0 )
        error("ball bearing contact(%d): static friction coefficient x must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("staticfrictioncoefficienty") )
        musy = HP.GetReal();
    else
        musy = musx;
    endif

    if ( musy <= 0 )
        error("ball bearing contact(%d): static friction coefficient y must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.Ms = diag([musx, musy]);

    if ( ~HP.IsKeyWord("slidingvelocitycoefficient") )
        error("ball bearing contact(%d): keyword sliding velocity coefficient expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif    

    elem.vs = HP.GetReal();

    if ( elem.vs <= 0 )
        error("ball bearing contact(%d): sliding velocity coefficient must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("slidingvelocityexponent") )
        elem.gamma = HP.GetReal();
    else
        elem.gamma = 1;
    endif

    if ( ~(HP.IsKeyWord("microslipstiffness") || HP.IsKeyWord("microslipstiffnessx")) )
        error("ball bearing contact(%d): keyword micro slip stiffness x expected at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    sigma0x = HP.GetReal();

    if ( sigma0x <= 0 )
        error("ball bearing contact(%d): micro slip stiffness x must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    if ( HP.IsKeyWord("microslipstiffnessy") )
        sigma0y = HP.GetReal();
    else
        sigma0y = sigma0x;
    endif

    if ( sigma0y <= 0 )
        error("ball bearing contact(%d): micro slip stiffness y must be greater than zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
    endif

    elem.sigma0 = diag([sigma0x, sigma0y]);

    if ( HP.IsKeyWord("microslipdamping") || HP.IsKeyWord("microslipdampingx") )
        sigma1x = HP.GetReal();

        if ( sigma1x < 0 )
            error("ball bearing contact(%d): micro slip damping x must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
        endif

        if ( HP.IsKeyWord("microslipdampingy") )
            sigma1y = HP.GetReal();

            if ( sigma1y < 0 )
                error("ball bearing contact(%d): micro slip damping y must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
            endif
        else
            sigma1y = sigma1x;
        endif
    else
        sigma1x = sigma1y = 0;
    endif

    elem.sigma1 = diag([sigma1x, sigma1y]);

    if ( HP.IsKeyWord("viscousfrictioncoefficient") || HP.IsKeyWord("viscousfrictioncoefficientx") )
        sigma2x = HP.GetReal();

        if ( sigma2x < 0 )
            error("ball bearing contact(%d): viscous friction coefficient x must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
        endif

        if ( HP.IsKeyWord("viscousfrictioncoefficienty") )
            sigma2y = HP.GetReal();

            if ( sigma2y < 0 )
                error("ball bearing contact(%d): viscous friction coefficient y must be greater than or equal to zero at line %s", pMbElem.GetLabel(), HP.GetLineData());
            endif
        else
            sigma2y = sigma2x;
        endif
    else
        sigma2x = sigma2y = 0;
    endif

    elem.sigma2 = diag([sigma2x, sigma2y]);

    elem.z = zeros(2,1);
    elem.zP = zeros(2,1);

	if (HP.IsKeyWord("initialstictionstate"))
		for i=1:2
			elem.z(i) = HP.GetReal();
		endfor
	endif

	if (HP.IsKeyWord("initialstictionstatederivative"))
		for i=1:2
			elem.zP(i) = HP.GetReal();
		endfor
	endif

    elem.rolling_friction = "no";

    if (HP.IsKeyWord("rollingfrictionmodel"))
        if (HP.IsKeyWord("SKF"))
            elem.rolling_friction = "SKF";
            if (~HP.IsKeyWord("staticloadcapacity"))
                error("ball bearing contact(%d): keyword \"staticloadcapacity\" expected at line %s\n", pMbElem.GetLabel(), HP.GetLineData());
            endif
            elem.C0 = HP.GetReal();
        elseif (HP.IsKeyWord("rigidball"))
            elem.rolling_friction = "rigid ball";
        elseif (HP.IsKeyWord("no"))
            elem.rolling_friction = "no";
        else
            error("ball bearing contact(%d): friction model not supported at line %s\n", pMbElem.GetLabel(), HP.GetLineData());
        endif
    endif

    elem = class(elem, "BallBearingContact");
endfunction
