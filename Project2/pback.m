clear 

% solving projected backward substitution 
function X = pback(U, b, g)
%     back substituion for defining T 
    T = back(U, b);
    X = ones(length(T), 1);
    
    for i = length(T):-1:1
%         computes the backward solution to find the equavlent value of T
        size = U(i,i+1:length(U));
        t_base = (b(i,1) - dot(size,X(i+1:length(X),1))) / U(i,i);
        
%         determines the correct value of X provided a T and g index value
        if t_base > g(i, 1)
            X(i,1) = t_base;
        else
            X(i,1) = g(i, 1);
        end
    end
end
