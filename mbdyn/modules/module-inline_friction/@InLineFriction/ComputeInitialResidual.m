% MBDyn (C) is a multibody analysis code. 
% http://www.mbdyn.org
%
% Copyright (C) 1996-2014
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

function f = ComputeInitialResidual(elem, XCurr)

    R1_0 = elem.pNode1.GetRRef();
    R2_0 = elem.pNode2.GetRRef();

%{
    X = [X1; 
         g1; 
         X1P;
         W1;
         X2;
         g2;
         X2P;
         W2;
         lambda;
         lambdaP];
%}

    X1      = XCurr( 1:3 );
    g1      = XCurr( 4:6 );
    X1P     = XCurr( 7:9 );
    W1      = XCurr(10:12);
    X2      = XCurr(13:15);
    g2      = XCurr(16:18);
    X2P     = XCurr(19:21);
    W2      = XCurr(22:24);
    lambda  = XCurr(25:26);
    lambdaP = XCurr(27:28);

    RDelta1 = rotation_increment_gibbs_rodriguez(g1);
    RDelta2 = rotation_increment_gibbs_rodriguez(g2);

    R1 = RDelta1 * R1_0;
    R2 = RDelta2 * R2_0;

    switch( typeinfo(R1) )
        case "matrix"
            assert(R1, elem.pNode1.GetRCurr(), sqrt(eps));
            assert(R2, elem.pNode2.GetRCurr(), sqrt(eps));
        case "gradient"
        otherwise
            error("unexpected type: %s", typeinfo(R1));
    endswitch

    F1 = R1 * ([elem.e2, elem.e3] * lambda);
    M1 = cross(X2 + R2 * elem.o2 - X1, F1);
    F2 = -F1;
    M2 = cross(R2 * elem.o2, F2);

    F1P = cross(W1, R1 * ([elem.e2, elem.e3 ] * lambda)) + R1 * ( [elem.e2, elem.e3] * lambdaP );
    M1P = cross(X2P + cross(W2, R2 * elem.o2) - X1P, F1) + cross(X2 + R2 * elem.o2 - X1, F1P);
    F2P = -F1P;
    M2P = cross(cross(W2, R2 * elem.o2), F2) + cross(R2 * elem.o2, F2P);

    c = [ elem.e2.'; 
          elem.e3.' ] * ( R1.' * ( X2 + R2 * elem.o2 - X1 ) - elem.o1 );
%{
    switch ( typeinfo(c) )
        case "matrix"
            disp("InitialAssRes:");
            disp("X2 + R2 * elem.o2 - X1="); disp(X2 + R2 * elem.o2 - X1);
            X1
            X2
            R1
            R2
            c
    endswitch
%}
    cP = [ elem.e2.';
           elem.e3.' ] * R1.' * ( cross(X2 + R2 * elem.o2 - X1, W1) + X2P + cross(W2, R2 * elem.o2) - X1P );

    f = [ F1;
          M1;
          F1P;
          M1P;
          F2;
          M2;
          F2P;
          M2P;
          c;
          cP ];
endfunction
