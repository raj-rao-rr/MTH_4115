clear 

x = (-2:0.1:2);
y1 = sin(x);
plot(x,y1)

hold on

y2 = chebyshev(x);
plot(x,y2)

hold off

function poly = chebyshev(x)
    % as instructed to modify the Program 3.3 for 3.3-C.1
    
    % calcuate based on Chebyshev theroem nodes [0, pi/2]
    b=(pi/4)+ (pi/4)*cos(2*(1:4)-1*pi/8); 
    yb=sin(b); 
    c=newtdd(b,yb,4);
    
    %For each input x, move x to the fundamental domain and evaluate
    % the interpolating polynomial
    s=1; % Correct the sign of sin
    x1=mod(x,2*pi);
    
    if x1>pi
        x1 = 2*pi-x1;
        s = -1;
    end
    
    if x1 > pi/2
        x1 = pi-x1;
    end
    
    poly = s*nest(3,c,x1,b);
end 


function c=newtdd(x,y,n)
    for j=1:n
        v(j,1)=y(j); % Fill in y column of Newton triangle
    end
    for i=2:n % For column i,
        for j=1:n+1-i % fill in column from top to bottom
            v(j,i)=(v(j+1,i-1)-v(j,i-1))/(x(j+i-1)-x(j));
        end
    end
    for i=1:n
        c(i)=v(1,i); % Read along top of triangle
    end % for output coefficients
end    
    

function y = nest(d,c,x,b)
    if nargin<4, b=zeros(d,1); end 
    y=c(d+1);
    for i=d:-1:1
      y = y.*(x-b(i))+c(i);
    end
end  
    