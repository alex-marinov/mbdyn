function elem = MyElem1(a, b, c, pDM)
    elem.a = a;
    elem.b = b;
    elem.c = c;
    elem.pNode1 = pDM.GetStructNode(1);
    elem.ridx = elem.pNode1.iGetFirstMomentumIndex() + int32(1:3);
    elem.cidx = elem.pNode1.iGetFirstPositionIndex() + int32(1:3);

    elem.S = [ elem.a, 0,      0;
               0,      elem.b, 0;
               0,      0,      elem.c ];

    elem = class(elem,"MyElem1");
endfunction
