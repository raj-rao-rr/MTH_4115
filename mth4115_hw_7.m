clear 


x = [0,2,3];
y = [1,2,4];
disp(polyinterp(x,y,3))


 function y0 = polyinterp(x,y,x0)
    m = newtondd(x,y);
    y0 = nest(m, x, x0);
 end

 function a = newtondd(x,y)
    n = length(x);
    a = transpose(repelem(0,length(x)));
    for j = 1:n
        a(j,1) = y(j);
    end
    
    for i = 2:n
        for j = 1:n+1-i
            a(j,i) = (a(j+1,i-1)-a(j,i-1))/(x(j+i-1)-x(j));
        end
    end
 end  
 
 function sum = nest(m, x, x0)
    coef = m(1,:);
    vars = x(1:length(x)-1); 
    step = 1;
    
    sum = coef(1);
    for i = coef(2:length(coef))
        selection = vars(1:step);
        value = i*prod(x0 - selection);
        sum = sum + value;
        step = step + 1;
    end
 end
   
