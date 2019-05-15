clear 

function [u] = bvp_solve_pj(r, g, g1, g2, L, U, P)
    r(1) = g1;
    r(length(r)) = g2;
%     solved by seeting projected back-substiution -> LUx=Pb, where b = r
    c = (L\P*r);
    u = pback(U, c, g);
end