clear 

function V = ode_sys_be(y0, t0, n, h, A)
    %     Initial refrence point value for starting the itteration
    w_initial = y0;
%     adjust for step size and count value from intial starting 
    for i = t0:h:((n-1)*h)+t0
%         refrence the explicit trapezoid method for itterative step
        new_A = (eye(length(w_initial))-h*A);
        w_next = new_A\w_initial;
        w_initial = w_next;
    end
%     returns the desired value
    V = w_next;
end