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

function y = my_func(x)
    global f;
    y = sin(2*pi*f*x);
    %printf("mboct: in my_func (x=%g, y=%g) ...\n", x, y);
endfunction
