function [f, ridx] = AssRes(elem, dCoef, pDM, pElem)
    global Time;

    ridx = elem.ridx;

    X1 = elem.pNode1.GetXCurr();

    f = -elem.S * X1;
endfunction
