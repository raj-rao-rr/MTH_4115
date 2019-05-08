clear 

function [u] = bvp_solve(r, g1, g2, L, U, P)
    r(1) = g1;
    r(length(r)) = g2;
%     solved by seeting Ax=b -> LUx=Pb, where b = r
    u = (L*U)\(P*r);  
end