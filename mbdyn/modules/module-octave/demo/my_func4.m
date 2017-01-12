# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2017
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
## AUTHOR: Reinhard Resch <r.resch@secop.com>
##        Copyright (C) 2011(-2017) all rights reserved.
##
##        The copyright of this code is transferred
##        to Pierangelo Masarati and Paolo Mantegazza
##        for use in the software MBDyn as described
##        in the GNU Public License version 2.1
##
##################################################################

function y = my_func4(t, pDM)
    persistent f = pDM.GetVariable("f");
    persistent pNode1 = pDM.GetStructNode(int32(1));
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
