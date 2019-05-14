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


% black-scholes solver using Crank-Nicolson (European Style)
function v = bs_eu_cn(N, M, K, T, r, q, sigma, Sm)
    h = Sm/N; k = T/M; s = ones(N+1,1);
    for i=0:N
        s(i+1,1) = i*h + 0;
    end
    
%     defining the standard C-column vectors
    c1 = (-k*(sigma^2)*(s.^2))/2;
    c2 = -k*(r-q)*s;
    c3 = ones(N+1,1)*(2+k*r);
    
%     creating the differentiation matrix 
    d1 = (k*(sigma^2)*(s.^2))/2;
    d2 = k*(r-q)*s;
    d3 = ones(N+1,1)*(2-k*r);
    
    [A, L, U, P] = bvp_matrix(0,Sm, c1, c2, c3, 1, 0, 1, 0);
    [B] = bvp_matrix(0,Sm, d1, d2, d3, 0, 0, 0, 0);
    v = max(K-s, 0);
    
%     itteratively solving using the backward-euler method 
    for i = 0:M-1
        g1 = exp(-r*(k+k*i))*K; g2 = 0;
        [v] = bvp_solve(B*v, g1, g2, L, U, P);      
    end
end


% black-scholes solver using backward euler (American Style)
function v = bs_am_be(N, M, K, T, r, q, sigma, Sm)
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
        g1 = K; g2 = 0;
        [v] = bvp_solve_pj(v, g, g1, g2, L, U, P);      
    end
end


% % % % % 
% Secoondary Helper Functions 
% % % % % 

% Boundary Value Problem solvers - back substition using Ax=b 
function [u] = bvp_solve(r, g1, g2, L, U, P)
    r(1) = g1;
    r(length(r)) = g2;
%     solved by seeting Ax=b -> LUx=Pb, where b = r
    u = (L*U)\(P*r);  
end

% Boundary Value Problem solvers - back substition using projection 
function [u] = bvp_solve_pj(r, g, g1, g2, L, U, P)
    r(1) = g1;
    r(length(r)) = g2;
%     solved by seeting Ax=b -> LUx=Pb, where b = r
    u = pback(L*U, P*r, g);  
end


% solving for LU=PA decomposition 
function [A, L, U, P] = bvp_matrix(a,b,c1, c2, c3, a1, b1, a2, b2)
    h = (b-a)/(length(c1)-1);
%      modified values to be used in re-defining the upper/lower rows
    u1 = a1 - (3*b1)/(2*h);
    u2 = (4*b1)/(2*h);
    u3 = -b1/(2*h);
    
    q3 = a2 + (3*b2)/(2*h);
    q2 = (-4*b2)/(2*h);
    q1 = b2/(2*h);
    
%     recreating the A matrix to be used with FDF2 and BDF2 
    A = zeros(length(c1));
    for i = 1:length(c1)-2
        r3 = c1(i+1)/(h^2) + c2(i+1)/(2*h);
        r2 = (-2*c1(i+1))/(h^2) + c3(i+1);
        r1 = c1(i+1)/(h^2) - c2(i+1)/(2*h);
        A(i+1,i:i+2) = [r1,r2,r3];
    end 
    
%     redefining both the first and last rows of the A-diff-matrix 
    A(1,1:3) = [u1, u2, u3];
    A(length(c1),length(c1)-2:length(c1)) = [q1, q2, q3];
    
    [L,U,P] = lu(A);
end


% solving projected backward substitution 
function X = pback(U, b, g)
%     back substituion for defining T 
    T = back(U, b);
    X = ones(length(T), 1);
    
    for i = length(T):-1:1
%         computes the backward solution to find the equavlent value of T
        size = U(i,i+1:length(U));
        t_base = (b(i,1) - dot(size,X(i+1:length(X),1))) / U(i,i);
        
%         determines the correct value of X provided a T and g index value
        if t_base > g(i, 1)
            X(i,1) = t_base;
        else
            X(i,1) = g(i, 1);
        end
    end
end