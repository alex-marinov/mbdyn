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

% AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
%        Copyright (C) 2011(-2019) all rights reserved.
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
    lambda = XCurr(13);
    z = XCurr(14:15);

    X1P = XPrimeCurr(1:3);
    g1P = XPrimeCurr(4:6);
    X2P = XPrimeCurr(7:9);
    g2P = XPrimeCurr(10:12);
    zP = XPrimeCurr(14:15);

    RDelta1 = rotation_increment_gibbs_rodriguez(g1);
    RDelta2 = rotation_increment_gibbs_rodriguez(g2);

    G1 = rotation_perturbation_matrix_gibbs_rodriguez(g1);
    G2 = rotation_perturbation_matrix_gibbs_rodriguez(g2);

    R1 = RDelta1 * R1_0;
    R2 = RDelta2 * R2_0;

    omega1 = G1 * g1P + RDelta1 * w1_0;
    omega2 = G2 * g2P + RDelta2 * w2_0;

    DeltaX = X2 + R2 * elem.o2 - X1 - R1 * elem.o1;
    DeltaXP = X2P + cross(omega2, R2 * elem.o2) - X1P - cross(omega1, R1 * elem.o1);

    DeltaXPb = [ elem.e1, elem.e2 ].' * R1.' * ( DeltaXP + cross(DeltaX, omega1) );

    tau = ( elem.sigma0 * z + elem.sigma1 * zP ) * abs(lambda) + elem.sigma2 * DeltaXPb;

    if ( sum(DeltaXPb.^2) == 0 )
        kappa = 0;
    else
        g = norm(elem.Mk^2 * DeltaXPb) / norm(elem.Mk * DeltaXPb) + (norm(elem.Ms^2 * DeltaXPb) / norm(elem.Ms * DeltaXPb) - norm(elem.Mk^2 * DeltaXPb) / norm(elem.Mk * DeltaXPb)) * exp(-(norm(DeltaXPb)/elem.vs)^elem.gamma);
    
        kappa = norm(elem.Mk^2 * DeltaXPb) / g;
    endif

    Phi = DeltaXPb - kappa * elem.Mk^-2 * elem.sigma0 * z - zP;
    
    F1 = R1 * ( [ elem.e1, elem.e2 ] * tau + elem.e3 * lambda);
    M1 = cross(X2 + R2 * elem.o2 - X1, F1);
    F2 = -F1;
    M2 = cross(R2 * elem.o2, F2);

    c = elem.e3.' * ( R1.' * ( X2 + R2 * elem.o2 - X1 ) - elem.o1 );

    f = [ F1;
          M1;
          F2;
          M2;
          c / dCoef;
          Phi ];
endfunction
