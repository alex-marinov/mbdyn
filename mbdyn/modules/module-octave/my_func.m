function y = my_func(x)
    global f;
    y = sin(2*pi*f*x);
    printf("in my_func (x=%g, y=%g) ...\n", x, y);
endfunction
