clear 


function calc = newton_method(initial_guess, epsilon)
    calc = initial_guess - equ(initial_guess)/deri(initial_guess);
    old_calc = 0;
    while abs(old_calc - calc) > epsilon
        old_calc = calc;
        calc = calc - equ(calc)/deri(calc); 
    end
end

function calc = secant_method(a, b, epsilon)
    calc = b - (equ(b)*(b-a))/(equ(b)-equ(a));
    old_calc = 0;
    while abs(old_calc - calc) > epsilon
        old_calc = calc;
        calc = = b - (equ(b)*(b-calc))/(equ(b)-equ(calc)) 
    end
end
