clear 

function [A, L, U, P] = bvp_matrix(a,b,c1, c2, c3, a1, b1, a2, b2)
    h = (b-a)/(length(c1)-1);
%      modified values to be used in re-defining the upper/lower rows
    u1 = a1 - (3*b1)/(2*h);
    u2 = (4*b1)/(2*h);
    u3 = -b1/(2*h);
    
    q3 = a2 + (3*b2)/(2*h);
    q2 = (-4*b2)/(2*h);
    q1 = b2/(2*h);
    
%     recreating the A matrix to be used with FDF2 and BDF2 
    A = zeros(length(c1));
    for i = 1:length(c1)-2
        r3 = c1(i+1)/(h^2) + c2(i+1)/(2*h);
        r2 = (-2*c1(i+1))/(h^2) + c3(i+1);
        r1 = c1(i+1)/(h^2) - c2(i+1)/(2*h);
        A(i+1,i:i+2) = [r1,r2,r3];
    end 
    
%     redefining both the first and last rows of the A-diff-matrix 
    A(1,1:3) = [u1, u2, u3];
    A(length(c1),length(c1)-2:length(c1)) = [q1, q2, q3];
    
    [L,U,P] = lu(A);
end