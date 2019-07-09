 % AUTHOR: Reinhard Resch <r.resch@a1.net>
        % Copyright (C) 2011(-2011) all rights reserved.

        % The copyright of this code is transferred
        % to Pierangelo Masarati and Paolo Mantegazza
        % for use in the software MBDyn as described
        % in the GNU Public License version 2.1
        
function y = gradnorm(x)
    y = sqrt(sum(sum(x.^2)));
endfunction
