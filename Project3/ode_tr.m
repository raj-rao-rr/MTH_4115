clear 

function V = ode_tr(y0, t0, n, h, f)
%     Initial refrence point value for starting the itteration
    w_initial = y0;
%     adjust for step size and count value from intial starting 
    for i = t0:h:((n-1)*h)+t0
%         refrence the explicit trapezoid method for itterative step
        w_next = w_initial + (h/2)*(f(i,w_initial) + f(i+h, w_initial+h*f(i,w_initial)));
        w_initial = w_next;
    end
%     returns the desired value
    V = w_next;
end
