clear 

A = randn(20, 10);
C = randn(5, 10);
b = randn(20, 1);
d = randn(5, 1);

%   constrained least-squares problems
xht = lsqlin(A,b,C,d);

%   Least norm via QR factorization
xln = transpose(C)*inv(C*transpose(C))*d;

if norm(A*xln - b)^2 > norm(A*xht - b)^2
    disp("Conditions are satisfied")
end