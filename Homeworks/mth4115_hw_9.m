clear 

rand_matrix = randn(30, 10);
rand_vector = randn(30, 1);


% iterates through to generate 3 random vectors to test residual condition 
for i = 1:3 
    d = randn(10, 1);
    [x_hat, residual] = least_square(rand_matrix, rand_vector);	
    % secondary residual based on ||A(x+d)-b||^2     
    sec_residual = norm((rand_matrix * (x_hat+d))-rand_vector)^2;
    if sec_residual > residual
        disp('Condition holds for given vector d')
    end
end

% formulazines the least squares system to minimize errror
function [x_hat, residual] = least_square(A, b)
	x_hat = (transpose(A)*A)\transpose(A)*b;
	residual = norm((A* x_hat)-b)^2;
end
