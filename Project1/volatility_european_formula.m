%   used to calucted the implied volatility from the european formula
function sigma = volatility_european_formula(id, V, K, t, S, q, r, s0, s1, tol, maxn)
    sigma = secant(s0, s1, @g, tol, maxn);
    
    function y = g(sigma)
        y = european_formula(id, K, t, S, sigma, q, r)-V;
    end
end


