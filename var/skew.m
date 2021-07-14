 % AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        % Copyright (C) 2011(-2019) all rights reserved.

        % The copyright of this code is transferred
        % to Pierangelo Masarati and Paolo Mantegazza
        % for use in the software MBDyn as described
        % in the GNU Public License version 2.1

function A = skew(a)
  persistent x = 1;
  persistent y = 2;
  persistent z = 3;
  
  A = [  0,     -a(z),   a(y);
         a(z),   0,     -a(x);
        -a(y),   a(x),   0     ];
endfunction
