clear 

[fe, emf, cond] = stats(10);
disp(fe)
disp(emf)
disp(cond)

%  calculates the forward error, error magnification factor and conditional
function [fe, emf, cond] = stats(n)
    [m,b] = spec_matrix(n);
    x = m\b;
    x_a = transpose(repelem(1,n));
 
    fe = max(abs(x-x_a));
    emf = (norm(x-x_a)*norm(b)) / (norm(b-(m*x_a))*norm(x));
    cond = m_norm(m)*m_norm(inv(m));
end

%  finds the infinte norm of a given vector 
function n = norm(v)
    n = max(abs(v));
end 

%  finds the infinte norm of a given matrix 
function n = m_norm(m)
    x = transpose(repelem(1,length(m)));
    for i = 1:length(m)
        val = sum(abs(m(i,:)));
        x(i) = val;
    end
    n = max(x);
end

% creates the special matrix and the b defined in problem set 
function [m,b] = spec_matrix(n)
    m = zeros(n);
    x = transpose(repelem(1,n));
    for i = 1:length(m)
        for j = 1:length(m)
            m(i,j) = 5/(i + 2*j - 1);
        end
    end
    b = m*x;
end


