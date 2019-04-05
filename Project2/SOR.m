clear 

% Successive over-relaxation method 
function x = SOR(A, b, x0, omega, tol, maxn)
%  defining the L,U,D forms to be used in the SOR decomposiiton method 
   L = dia(tril(A));
   D = diag(diag(A));
   U = dia(triu(A));
  
   error = 1;
   
%  determining the optimal method based on itterative block and tolerance 
   for i = 1:maxn
%      conditional tolerance break defined by relative error between x-vals
       if error > tol
           x1 = (inv(omega*L+D))*((1-omega)*D*x0 - omega*U*x0)+omega*(inv(D+omega*L))*b; 
           error = abs(max(x1-x0))/abs(max(x1));
           x0 = x1;
       else
           break
       end       
   end
   x = x1;
end

% helper function to be used when extracting the diagonal values
function d = dia(U)
    d = U;
%   sets all diagonals entries to be zero from a matrix (used for SOR-LU)
    for i = 1:length(U)
        d(i,i) = 0;
    end
end
