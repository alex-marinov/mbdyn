function y = my_func(x)
    global f;
    y = sin(2*pi*f*x);
    printf("in my_func (f=%g, x=%g, y=sin(2*pi*f*x)=sin(2*pi*%g*%g)=%g) ...\n", f, x, f, x, y);
endfunction
