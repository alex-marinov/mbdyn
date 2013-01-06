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

function f = ComputeResidual(elem, dCoef, XCurr, XPrimeCurr)

    R1_0 = elem.pNode1.GetRRef();
    R2_0 = elem.pNode2.GetRRef();

    w1_0 = elem.pNode1.GetWRef();
    w2_0 = elem.pNode2.GetWRef();

    X1 = XCurr(1:3);
    g1 = XCurr(4:6);
    X2 = XCurr(7:9);
    g2 = XCurr(10:12);
    lambda = XCurr(13:14);
    z = XCurr(15);
%{
    if ( elem.pDM.dGetTime() == 0)
        lambda
    endif
%}
    X1P = XPrimeCurr(1:3);
    g1P = XPrimeCurr(4:6);
    X2P = XPrimeCurr(7:9);
    g2P = XPrimeCurr(10:12);
    zP = XPrimeCurr(15);

    RDelta1 = rotation_increment_gibbs_rodriguez(g1);
    RDelta2 = rotation_increment_gibbs_rodriguez(g2);

    G1 = rotation_perturbation_matrix_gibbs_rodriguez(g1);
    G2 = rotation_perturbation_matrix_gibbs_rodriguez(g2);

    R1 = RDelta1 * R1_0;
    R2 = RDelta2 * R2_0;

    W1 = G1 * g1P + RDelta1 * w1_0;
    W2 = G2 * g2P + RDelta2 * w2_0;
%{
    switch( typeinfo(R1) )
        case "matrix"
            assert(R1, elem.pNode1.GetRCurr(), sqrt(eps));
            assert(R2, elem.pNode2.GetRCurr(), sqrt(eps));
        case "gradient"
        otherwise
            error("unexpected type: %s", typeinfo(R1));
    endswitch
%}
%{
    elem.pNode2.GetgCurr()
    elem.pNode2.GetWRef()
    elem.pNode2.GetWCurr()
%}

    switch( typeinfo(W1) )
        case "matrix"
            try
                assert(W1, elem.pNode1.GetWCurr(), sqrt(eps));
                assert(W2, elem.pNode2.GetWCurr(), sqrt(eps));
            catch
                fprintf(stderr, "%g: %s\n", elem.pDM.dGetTime(), lasterror.message);
            end_try_catch
        case "gradient"
        otherwise
            error("unexpected type: %s", typeinfo(W1));
    endswitch
%{%}

    DeltaXP = ComputeDeltaXP(elem, X1, R1, X1P, W1, X2, R2, X2P, W2);

    if ( lambda(1) == 0 ) % avoid division by zero during automatic differentiation
        lambda_res = abs(lambda(2));
    elseif ( lambda(2) == 0 )
        lambda_res = abs(lambda(1));
    else
        lambda_res = sqrt(sum(lambda.^2));
    endif

    tau = ComputeTau(elem, DeltaXP, lambda_res, z);

    if ( DeltaXP > 0 )
        sign_DeltaXP = 1;
    elseif ( DeltaXP < 0 )
        sign_DeltaXP = -1;
    else
        sign_DeltaXP = 0;
    endif

    Phi = zP - DeltaXP * ( 1 - z / elem.delta * sign_DeltaXP );
    
    F1 = R1 * (elem.e1 * tau + [elem.e2, elem.e3] * lambda);
    M1 = cross(X2 + R2 * elem.o2 - X1, F1);
    F2 = -F1;
    M2 = cross(R2 * elem.o2, F2);

    c = [ elem.e2.'; 
          elem.e3.' ] * ( R1.' * ( X2 + R2 * elem.o2 - X1 ) - elem.o1 );

%{
    switch ( typeinfo(c) )
        case "matrix"
            disp("AssRes:");
            disp("X2 + R2 * elem.o2 - X1="); disp(X2 + R2 * elem.o2 - X1);
            X1
            X2
            R1
            R2
            c
    endswitch
%}

    f = [ F1;
          M1;
          F2;
          M2;
          c / dCoef;
          Phi * elem.PhiScale ];
endfunction
