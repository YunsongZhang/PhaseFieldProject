% Compute Laplace of u on domain [-Lx, Lx]*[-Ly, Ly]
function uu = lap(u, m, n, Lx, Ly)

k1          = [0:m/2, -(m/2-1):-1];
k2          = [0:n/2, -(n/2-1):-1];
[k1,k2] = meshgrid(k1,k2);
k1_sqr     = k1.^2;
k2_sqr     = k2.^2;

u1  = fft2(u);
uu1  = u1.*( -(pi/Lx)^2*k1_sqr - (pi/Ly)^2*k2_sqr );
uu = ifft2(uu1);
