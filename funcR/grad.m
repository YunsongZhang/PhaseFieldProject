% Derivative of u = u(x,y) on domain [-Lx, Lx]*[-Ly, Ly]
function [ux, uy] = grad(u, m, n, Lx, Ly)

k1 = [0:m/2-1, 0,  -(m/2-1):-1];
k2 = [0:n/2-1, 0,  -(n/2-1):-1];
[k1, k2] = meshgrid(k1, k2);


u1  = fft2(u);
ux  = ifft2(u1.*k1*1i)*pi/Lx;
uy  = ifft2(u1.*k2*1i)*pi/Ly;
