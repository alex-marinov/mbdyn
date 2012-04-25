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

function y = my_func4(t, pDM)
    persistent f = pDM.GetVariable("f");
    persistent pNode1 = pDM.GetStructNode(1);
    persistent MBDYN_VERSION = pDM.GetVersion();

    %printf("mboct: MBDyn version %s\n", MBDYN_VERSION);
    %disp(typeinfo(pDM));

    X = pNode1.GetXCurr();
    R = pNode1.GetRCurr();
    V = pNode1.GetVCurr();
    W = pNode1.GetWCurr();

    XPP = pNode1.GetXPPCurr();
    WP = pNode1.GetWPCurr();

    [phi] = rotation_matrix_to_euler123(R);

%{
    disp("mboct: X=");     disp(X);
    disp("mboct: phi=");   disp(phi*180/pi);
    disp("mboct: V=");     disp(V);
    disp("mboct: W=");     disp(W);
%}
    y = norm(X);
endfunction
