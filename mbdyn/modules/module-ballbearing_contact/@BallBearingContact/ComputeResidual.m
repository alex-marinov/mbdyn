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

function f = ComputeResidual(elem, dCoef, XCurr, XPrimeCurr)

    R1_0 = elem.pNode1.GetRRef();
    R2_0 = elem.pNode2.GetRRef();

    w1_0 = elem.pNode1.GetWRef();
    w2_0 = elem.pNode2.GetWRef();

    X1 = XCurr(1:3);
    g1 = XCurr(4:6);
    X2 = XCurr(7:9);
    g2 = XCurr(10:12);
    z = XCurr(13:14);

    X1P = XPrimeCurr(1:3);
    g1P = XPrimeCurr(4:6);
    X2P = XPrimeCurr(7:9);
    g2P = XPrimeCurr(10:12);
    zP = XPrimeCurr(13:14);

    RDelta1 = rotation_increment_gibbs_rodriguez(g1);
    RDelta2 = rotation_increment_gibbs_rodriguez(g2);

    G1 = rotation_perturbation_matrix_gibbs_rodriguez(g1);
    G2 = rotation_perturbation_matrix_gibbs_rodriguez(g2);

    R1 = RDelta1 * R1_0;
    R2 = RDelta2 * R2_0;

    omega1 = G1 * g1P + RDelta1 * w1_0;
    omega2 = G2 * g2P + RDelta2 * w2_0;

    switch( typeinfo(R1) )
        case "matrix"
            assert(R1, elem.pNode1.GetRCurr(), sqrt(eps));
            assert(R2, elem.pNode2.GetRCurr(), sqrt(eps));
        case "gradient"
        otherwise
            error("unexpected type: %s", typeinfo(R1));
    endswitch
%{%}
%{
    elem.pNode2.GetgCurr()
    elem.pNode2.GetWRef()
    elem.pNode2.GetWCurr()
%}

    switch( typeinfo(omega1) )
        case "matrix"
            assert(omega1, elem.pNode1.GetWCurr(), sqrt(eps));
            assert(omega2, elem.pNode2.GetWCurr(), sqrt(eps));
        case "gradient"
        otherwise
            error("unexpected type: %s", typeinfo(W1));
    endswitch
%{%}
	b = R2.' * (X1 - X2) - elem.o2;
    n = elem.Rt2(:, 3).' * b;
	v = elem.Rt2.' * b - [0; 0; n];

	d = elem.R - n;

	if (d > 0)
		norm_Fn = 4/3 * elem.E * sqrt(elem.R) * d^(3/2);
	else
		norm_Fn = 0;
	endif

	c1P = X1P - cross(omega1, R2 * elem.Rt2(:, 3) * elem.R);
	c2P = X2P + cross(omega2, R2 * (elem.o2 + elem.Rt2 * v));

    uP = elem.Rt2(:, 1:2).' * R2.' * (c1P - c2P);

    tau = norm_Fn * (elem.sigma0 * z + elem.sigma1 * zP) + elem.sigma2 * uP;

    if ( sum(uP.^2) == 0 )
        kappa = 0;
    else
        g = norm(elem.Mk^2 * uP) / norm(elem.Mk * uP) + (norm(elem.Ms^2 * uP) / norm(elem.Ms * uP) - norm(elem.Mk^2 * uP) / norm(elem.Mk * uP)) * exp(-(norm(uP)/elem.vs)^elem.gamma);
    
        kappa = norm(elem.Mk^2 * uP) / g;
    endif

    Phi = uP - kappa * elem.Mk^-2 * elem.sigma0 * z - zP;
    
	F1 = R2 * elem.Rt2 * [-tau; norm_Fn];
	M1 = -cross(R2 * elem.Rt2(:, 3) * n, F1);
	F2 = -F1;
	M2 = cross(X1 - X2 - R2 * elem.Rt2(:, 3) * n, F2);

    switch(elem.rolling_friction)
    case "no"
    otherwise
        wP = elem.Rt2(:, 1:2).' * R2.' * (X1P - (X2P + cross(omega2, X1 - X2)));

        switch (typeinfo(wP))
        case "matrix"
            gamma = atan2(wP(2), wP(1));
        case "gradient"
            % FIXME: atan2 is not supported by the octave AD package
            gamma = 0 * wP(1);
            gamma.x = atan2(wP(2).x, wP(1).x);
        endswitch

        switch(elem.rolling_friction)
        case "SKF"
            fn = 0.0016 * (norm_Fn / elem.C0)^0.33 * norm_Fn;
        case "rigid ball"
            if (d > 0)
                fn = sqrt(d/(2*elem.R)) * norm_Fn;
            else
                fn = 0;
            endif
        case "no"
            fn = 0;
        otherwise
            error("friction model not implemented");
        endswitch
%{
        switch(typeinfo(fn))
        case "scalar"
            printf("fn/norm_Fn=%g\n", fn/norm_Fn);
        endswitch
%}
        if (sum(wP.^2) == 0)
            fv = 0;
        else
            fv = tanh(2 * pi * norm(wP) / elem.vs);
        endif

        Ff1 = -R2 * (elem.Rt2(:, 1) * cos(gamma) + elem.Rt2(:, 2) * sin(gamma)) * fn * fv;
        Mf1 = zeros(3, 1);
        Ff2 = -Ff1;
        Mf2 = cross(X1 - X2, Ff2);

        F1 += Ff1;
        M1 += Mf1;
        F2 += Ff2;
        M2 += Mf2;
    endswitch
%{    
    switch (typeinfo(c))
        case "matrix"
            c/dCoef
        case "gradient"
            c.x/dCoef
    endswitch
%}
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
          Phi ];
endfunction
