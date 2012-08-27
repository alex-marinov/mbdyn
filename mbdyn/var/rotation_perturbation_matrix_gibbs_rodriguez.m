 % AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        % Copyright (C) 2011(-2011) all rights reserved.

        % The copyright of this code is transferred
        % to Pierangelo Masarati and Paolo Mantegazza
        % for use in the software MBDyn as described
        % in the GNU Public License version 2.1
        
function G = rotation_perturbation_matrix_gibbs_rodriguez(g)
    G = 4 / (4 + g.' * g) * (eye(3) + 1/2 * skew(g));
endfunction
