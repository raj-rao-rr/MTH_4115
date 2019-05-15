
% Problem 1
clear
N = 2^8; M = 2^7; 
q = 0; r = 0.01; T = 1; K = 50; sigma = 0.1;
ve = 1.919766444277223;
v = bs_eu_be(N, M, K, T, r, q, sigma, 2*K);
error_1 = norm(v(N/2)-ve);

% Problem 2
clear
N = 2^8; M = 2^7; 
q = 0; r = 0.01; T = 1; K = 50; sigma = 0.1;
ve = 1.921673472665866;
v = bs_eu_cn(N, M, K, T, r, q, sigma, 2*K);
error_2 = norm(v(N/2)-ve);


% Problem 3
clear
N = 2^8; M = 2^7; 
q = 0; r = 0.01; T = 1; K = 50; sigma = 0.1;
ve = 1.964605919902253;
v = bs_am_be(N, M, K, T, r, q, sigma, 2*K);
error_3 = norm(v(N/2)-ve);
