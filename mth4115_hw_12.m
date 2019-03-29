clear 

disp(trapezoid_approx(@F, 0, 1, 32))

function approx = trapezoid_approx(f, a, b, m)
    h = (b-a)/m;
    sum = f(a)+f(b);
    
    for i = a+(b-a)/m:(b-a)/m:b-(b-a)/m
       sum = sum + 2*f(i); 
    end
    
    approx = (h/2)*sum;
    format long
end

function v = F(x)
%     v = x/(sqrt(x^2+9));
%     v = x^3/(x^2+1);
%     v = x*exp(x);
%     v = x^2*log(x);
%     v = x^2*sin(x);
%     v = x^3/(sqrt(x^4-1));
%     v = 1/(sqrt(x^2+4));
    v = x/(sqrt(x^4+1));
end
