clear 

% black-scholes solver using backward euler (European Style)
function v = bs_eu_be(N, M, K, T, r, q, sigma, Sm)
    h = Sm/N; k = T/M; s = ones(N+1,1);
    for i=0:N
        s(i+1,1) = i*h + 0;
    end
    
%     defining the standard C-column vectors
    c1 = (-k*(sigma^2)*(s.^2))/2;
    c2 = -k*(r-q)*s;
    c3 = ones(N+1,1)*(1+k*r);
    
    [A, L, U, P] = bvp_matrix(0,Sm, c1, c2, c3, 1, 0, 1, 0);
    v = max(K-s, 0);
    
%     itteratively solving using the backward-euler method 
    for i = 0:M-1
        g1 = exp(-r*(k+k*i))*K; g2 = 0;
        [v] = bvp_solve(v, g1, g2, L, U, P);      
    end
end