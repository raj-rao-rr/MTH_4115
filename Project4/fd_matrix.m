clear

function [Lf, A] = fd_matrix(a,b,c1,c2,c3,f)
    Lf = zeros(length(f), 1);
    A = zeros(length(f));
    h = (b-a)/(length(f)-1);
    
    % itterated to form the A matrix, given column vectors c1, c2, c3
    for i = 1:length(f)-2
        r3 = c1(i+1)/(h^2) + c2(i+1)/(2*h);
        r2 = (-2*c1(i+1))/(h^2) + c3(i+1);
        r1 = c1(i+1)/(h^2) - c2(i+1)/(2*h);
        A(i+1,i:i+2) = [r1,r2,r3];
    end 
%     define the A matrix as itterative above and multpily function vector
    Lf = A*f;
end