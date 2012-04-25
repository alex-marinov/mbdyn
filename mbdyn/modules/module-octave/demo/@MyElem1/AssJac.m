function [Jac, ridx, cidx] = AssJac(elem, dCoef, pDM)
    Jac = dCoef * elem.S;
    ridx = elem.ridx;
    cidx = elem.cidx;
endfunction
