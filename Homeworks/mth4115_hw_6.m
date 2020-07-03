clear 

n=5;
[m, b] = spare_matrix(n);
x = transpose(repelem(0,n));
[new_x, steps] = jacobi(m,b,x);
disp(new_x)
disp(steps)

% jacobi method for solving systems 
function [x, steps] = jacobi(a,b,x)
    absVal = inf;
    steps = 0;
    
    while absVal > 1e-6
        for i = 1:length(a)
            selc = [a(i,1:i-1), a(i,i+1:length(a))];
            x_selc = vertcat(x(1:i-1,1), x(i+1:length(x),1));
            val = selc*x_selc;
            x(i,1) = (b(i,1) - val)/a(i,i);
        end
        steps = steps + 1;
        absVal = max(abs(b - x));
    end
end

% Sparse matrix setup
function [m,b] = spare_matrix(size)
    m = eye(size) * 3;
    for i = 1:size-1
        m(i,i+1) = -1;
        m(i+1,i) = -1;
    end
    
    b = transpose(repelem(1,size));
    for x = 1:size
        if mod(x,2) == 0
            b(x, 1) = 1;
        else
            b(x, 1) = 2;
        end
    end
end 