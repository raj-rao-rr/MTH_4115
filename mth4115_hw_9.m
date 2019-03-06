clear 

rand_matrix = randn(30, 30);
rand_vector = randn(30, 1);


for i = 1:3 % iterates through to generate 3 random vectors to test residual condition 
    di = randn(10, 1);
    [x_hat, residual] = least_square(rand_matrix, rand_vector);	
    sec_residual = norm((A* (x_hat+di))-b)^2;
    if sec_residual > residual
        disp('Condition holds for given vector d')
        disp(d)
    end
end


function [x_hat, residual] = least_square(A, b)
	x_hat = (transpose(A)*A)\transpose(A)*b;
	residual = norm((A* x_hat)-b)^2;
end
