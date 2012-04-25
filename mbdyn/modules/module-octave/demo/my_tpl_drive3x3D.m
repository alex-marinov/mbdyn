function J = my_tpl_drive3x3D(t)
    J11 = 10 - 0.01 * t;
    J22 = J11;
    J33 = J11;

    J = [ J11, 0, -0.3*J11;
          0,   J22, 0;
          -0.1*J11,   0, J33 ];
endfunction
