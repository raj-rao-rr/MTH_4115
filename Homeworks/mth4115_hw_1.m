clear
% format long in commmand line prior to executing for precision
func = @f;
disp(fzero(func,0.1))


% used for ploting purposes (un-comment your desired equation)
function p = modify_plot(range) 
    v = -range:range;
%     equ = (2*v.^3) - (6*v) - 1;
%     equ = exp(v-2) + (v.^3) - v;
%     equ = 1 + (5*v) - (6*v.^3) - exp(2*v);
    plot(v, equ)
    p = "Plotted Function";
end


function [calc, count] = bisect_method(a, b, epsilon)
    %   weight the appropaite values of upper and lower
    if equation(a) > equation(b)
        upper = a;
        lower = b;
    else
        upper = b;
        lower = a;
    end

    % iterates over to find a solution while calling the equation function
    c = (upper + lower)/2;
    count = 0;
    %   bounded range by arbitrary epsilon for root finder
    while (equation(c) > epsilon)||(-epsilon > equation(c))  
        count = count + 1;
        c = (upper + lower)/2;
        %   decision process for bisection process
        if equation(c) > 0
            upper = c;
        elseif equation(c) < 0
            lower = c;
        end
    end
    calc = c;
end


function [calc, count] = fixed_itter(start_val, epsilon)
    c = start_val;
    old = c;
    c = g(c);
    count = 0;
    
    %   itterative function for succesively dimisnhing the fixed-point
    while abs(old-c) > epsilon
        count = count + 1;
        old = c;
        c = g(c);
    end
    calc = c;
end


% used to approximate the value of squares provided a guess and epsilon
function [square, count] = sq_root(approx, guess, epsilon)
    count=0;
    square = (guess + (approx/guess))/2;
    %   range bounded epsilon checking adajcent sets
    while abs((sqrt(approx)-square)) > epsilon
        count = count + 1;
        square = (square +(approx/square))/2;
    end   
end


% reformated for fixed-point itteration (un-comment your desired equation)
function ret = g(v)
    ret = (2*v+2)^(1/3);
%     ret = log(7-v);
%     ret = log(4-sin(v));
end


% refrence equations for root-solving (un-comment your desired equation)
function equ = equation(v)
    equ = v^3-9;
%     equ = (3*v^3) + (v^2) - v - 5;
%     equ = (cos(v))^2 + 6 - v;
%     equ = (2*v^3) - (6*v) - 1;
%     equ = exp(v-2) + (v^3) - v;
%     equ = 1 + (5*v) - (6*v^3) - exp(2*v);
%     equ = v^3-2;
%     equ = v^3-3;
%     equ = v^3-5;
end




