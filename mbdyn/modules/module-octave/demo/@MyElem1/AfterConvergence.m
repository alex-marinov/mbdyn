# $Header$
#
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2012
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
##        Copyright (C) 2011(-2012) all rights reserved.
##
##        The copyright of this code is transferred
##        to Pierangelo Masarati and Paolo Mantegazza
##        for use in the software MBDyn as described
##        in the GNU Public License version 2.1
##
##################################################################

function [elem] = AfterConvergence(elemin, XCurr, XPrimeCurr)
    elem = elemin;
    iFirstIndex = elem.pMbElem.iGetFirstIndex() + int32(1);

    elem.t(++elem.nIter) = elem.pDM.dGetTime();
    elem.X2(elem.nIter) = XCurr.dGetCoef(iFirstIndex);
    elem.XP2(elem.nIter) = XPrimeCurr.dGetCoef(iFirstIndex);

    if ( false )
        if ( 0 == mod(elem.nIter, 200) )
            figure(1);
            subplot(2,1,1);
            plot(elem.t, elem.X2,'-;X2(t);1');
            xlabel('t [s]');
            ylabel('X');
            grid on;
            grid minor on;
            subplot(2,1,2);
            plot(elem.t, elem.XP2, '-;XP2(t);1');
            xlabel('t [s]');
            ylabel('XP');
            grid on;
            grid minor on;
            drawnow();
        endif
    endif
endfunction
