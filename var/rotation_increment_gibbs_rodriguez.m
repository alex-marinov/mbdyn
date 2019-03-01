 % AUTHOR: Reinhard Resch <r.resch@secop.com>
        % Copyright (C) 2011(-2011) all rights reserved.

        % The copyright of this code is transferred
        % to Pierangelo Masarati and Paolo Mantegazza
        % for use in the software MBDyn as described
        % in the GNU Public License version 2.1
        
function R = rotation_increment_gibbs_rodriguez(g)
    R = eye(3) + 4/(4 + g.' * g) * ( skew(g) + 1/2 * skew(g) * skew(g) );
endfunction
