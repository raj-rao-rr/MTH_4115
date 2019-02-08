clear 

table_form()

function tb = table_form()
    x_arr = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    org_arr = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    mod_arr = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];   
    
    for i = 1:14
       x = 10^(-i);
       x_arr(i,1) = x;
       org_arr(i,1) = equ(x);
       mod_arr(i,1) = mod(x);
    end    
    
    tb = table(x_arr,org_arr,mod_arr);
end

function val = equ(x)
%     val = (1- sec(x))/(tan(x))^2;
    val = (1-(1-x)^3)/x;
end

function val = mod(x)
%     val = -1/(1+sec(x));
    val = (x^2)-(3*x)+3;
end