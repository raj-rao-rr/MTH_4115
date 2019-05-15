clear 

function X = pback(U, b, g)
    X = zeros(length(U), 1);
    
    for i = length(U):-1:1
%         computes the backward solution to find the equivlent value of U
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