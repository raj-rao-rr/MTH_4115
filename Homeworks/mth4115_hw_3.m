clear 

m = [[6,1];[3,7]];
[L , U] = lu_decomp(m);
disp(L)
disp(U)

matrix = [[2,-2,-1,-2];[4,1,-2,1];[-2,1,-1,-3]];
% ep_form = naive_ge(matrix);
% disp(backward_sub(ep_form))


% perform LU decomposition on an [m X m] matrix 
function [L,U] = lu_decomp(m)
    L = eye(length(m));
    for j = 1 : length(m)-1
%       if zero pivot encountered end the script 
        if abs(m(j,j))<1e-6; error("zero pivot encountered"); end
%       for each row below the jth row being active to end of matrix 
        for i = j+1 : size(m,1)
%       multiplication/scalar factor for performing G.E. 
            mult = m(i,j)/m(j,j);
            L(i,j) = mult;
            wipe = m(j,:).*mult;
            m(i,:) = m(i,:) - wipe;
        end
    end
    U = m;
end


% perform Guassian elimination on an [m X m+1] matrix 
function arr = naive_ge(m)
%   for each column up to the n-1 ** must retain final nth pivot
    for j = 1 : length(m)-1
%       if zero pivot encountered end the script 
        if abs(m(j,j))<1e-6; error("zero pivot encountered"); end
%       for each row below the jth row being active to end of matrix 
        for i = j+1 : size(m,1)
%       multiplication/scalar factor for performing G.E. 
            mult = m(i,j)/m(j,j);
            wipe = m(j,:).*mult;
            m(i,:) = m(i,:) - wipe;
        end
    end
    arr = m;   
end


function vals = backward_sub(m)
%   starting from the last row and working upward to substitute 
    x = [];
    for i = size(m,1) : -1 : 1
%       since we solve for the (ith,ith) value we start after ith place
%       when eliminating to the length of the vector 
        for j = i+1 : length(m)-1
%           check if pivot is terminal then simply terminate loop
            if j == length(m)
                break
            else
                m(i, length(m)) = m(i, length(m)) - m(i,j)*x(j);
            end
        end
        x(i) =  m(i, length(m))/m(i,i);
    end
    vals = x;
end
