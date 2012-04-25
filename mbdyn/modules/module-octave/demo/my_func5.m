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

function F = my_func5(t)
    global f;
    global F0;
    F = F0 * [ sin(2*pi*f*t); 
               cos(2*pi*f*t); 
               0 ];
endfunction
