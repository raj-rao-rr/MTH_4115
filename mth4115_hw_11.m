clear 

table_form()

function tb = table_form()
    exact = cos(0)+sin(0);
    for i = 1:12
       x = 10^(-i);
       h(i,1) = x;
       
       error(i,1) = abs(three_point(0,x) - exact);
    end    
    plot(h, error)
    tb = table(h,error);
end

function val = three_point(x, h)
    val=(f(x+h)-f(x-h))/(2*h);
end

function val = f(x)
    val = sin(x)-cos(x);
end