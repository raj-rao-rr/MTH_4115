% secant method for approximating the value of a given root 
function x = secant(x0, x1, f, tol, maxn)
x(1)=x1;
x(2)=x(1)-f(x(1))*(x1-x0)/(f(x1)-f(x0));
i=2;
n=0;
%      range bounded tollerance for error between checking adjacent values
    while (abs(x(i)-x(i-1))>tol)&&(n<=maxn)
%      checking for the number of itterations provided bound   
        x(i+1)=x(i)-f(x(i))*(x(i)-x(i-1))/(f(x(i))-f(x(i-1)));
        i=i+1;
        n=n+1;
    end
x = x(i);
end


