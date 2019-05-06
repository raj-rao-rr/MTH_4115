clear 

function V = ode_high_rk4(y0, t0, n, h, f)
    %     Initial refrence point value for starting the itteration
    w_initial = y0;
%     adjust for step size and count value from intial starting 
    for i = t0:h:((n-1)*h)+t0
%         refrence the RK4 method for itterative step
        s1 = f(i, w_initial);
        s2 = f(i + (h/2), w_initial + (h/2)*s1);
        s3 = f(i + (h/2), w_initial + (h/2)*s2);
        s4 = f(i + h, w_initial + h*s3);
        w_next = w_initial + (h/6)*(s1 + 2*s2 + 2*s3 + s4);
        w_initial = w_next;
    end
%     returns the desired value
    V = w_next;
end
