clear 
% format long in commmand line prior to executing for precision
val = secant_method(1,2, 0.00000001);
disp(val)

func = @g;
disp(fzero(func,[-2,1]))

function calc = newton_method(initial_guess, epsilon)
    calc = initial_guess - equ(initial_guess)/deri(initial_guess);
    old_calc = 0;
    %   range bounded epsilon checking adajcent sets
    while abs(old_calc - calc) > epsilon
        old_calc = calc;
        calc = calc - equ(calc)/deri(calc); 
    end
end


function calc = secant_method(a, b, epsilon)
    calc = b - (equ(b)*(b-a))/(equ(b)-equ(a));
    old_calc = 0;
    %   range bounded epsilon checking adajcent sets
    while abs(old_calc - calc) > epsilon
        old_calc = calc;
        calc = b - (equ(b)*(b-calc))/(equ(b)-equ(calc)); 
    end
end


% refrence equations for root-solving (un-comment your desired equation)
function value = equ(x)
%     value = (x^3) - (2*x) - 2;
%     value = exp(x) + x - 7;
    value = exp(x) + sin(x) - 4;
%     value = (1-(3/(4*x)))^(1/3);
end


% derivatives for root-solving (un-comment your desired equation)
function value = deri(x)
%     value = (3*x^2) - 2;
%     value = exp(x) + 1;
    value = exp(x) + cos(x);
%     value = 1/((4)^(1/3)*x^(4/3)*(4*x-3)^(2/3));
end
