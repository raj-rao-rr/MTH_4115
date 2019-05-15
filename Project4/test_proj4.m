% Problem 1
clear
a=-pi/3; b=pi/3;
a1=2; b1=3; a2=4; b2=-5; g1=-sqrt(3)+3/2; g2=2*sqrt(3)-5/2;
n = 100; m = n+1; h = (b-a)/n;
x = [a:h:b]';
c1 = -4*x.*x; c2 = 3*x; c3 = 2*ones(m,1);
f = sin(x);  Lfe = c1.*(-f) + c2.*cos(x) + c3.*f;
[Lf A] = fd_matrix(a,b,c1,c2,c3,f);
er = Lf-Lfe;
error_1 = norm(er(2:end-1),'inf')


% Problem 2
clear
a=-pi/3; b=pi/3;
a1=2; b1=3; a2=4; b2=-5; g1=-sqrt(3)+3/2; g2=2*sqrt(3)-5/2;
n = 100; m = n+1; h = (b-a)/n;
x = [a:h:b]';
c1 = -4*x.*x; c2 = 3*x; c3 = 2*ones(m,1);
[ A L U P] = bvp_matrix(a,b,c1,c2,c3,a1,b1,a2,b2);
f = sin(x);  Lfe =  c1.*(-f) + c2.*cos(x) + c3.*f;
Lfe(1) = g1; Lfe(end) = g2;
Lf = A*f;
error_2 = norm(Lf-Lfe,'inf')



% Problem 3
clear
a=-pi/3; b=pi/3;
a1=2; b1=3; a2=4; b2=-5; g1=-sqrt(3)+3/2; g2=2*sqrt(3)-5/2;
n = 100; m = n+1; h = (b-a)/n;
x = [a:h:b]';
c1 = -4*x.*x; c2 = 3*x; c3 = 2*ones(m,1);
ue = sin(x);  
r = c1.*(-ue) + c2.*cos(x) + c3.*ue;
[A L U P ] = bvp_matrix(a,b,c1,c2,c3,a1,b1,a2,b2);
[u] = bvp_solve(r,g1,g2,L,U,P);
error_3 = norm(u-ue,'inf')
