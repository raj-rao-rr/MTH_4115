clear 

% formulazing the backward substituion of a given upper triangular matrix
function v = back(U, b)
    m = [U b];
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
    v = x;
end