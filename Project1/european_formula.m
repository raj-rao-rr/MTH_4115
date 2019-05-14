
%   closed form price for a european option 
function V = european_formula(id, K, t, S, sigma, q, r)
%   id = 1 for call and 0 for put
    d1 = (log(S/K)+t*(r-q+0.5*sigma^2))/(sigma*sqrt(t));
    d2 = d1 - sigma*sqrt(t);
%     call price
    if id == 1
        V = S*exp(-q*t)*cdf('Normal',d1,0,1)-K*exp(-r*t)*cdf('Normal',d2,0,1);
    else
%      put price
        V = K*exp(-r*t)*cdf('Normal',-d2,0,1)-S*exp(-q*t)*cdf('Normal',-d1,0,1);
    end
end
