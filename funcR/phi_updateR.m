function phi=phi_updateR(phi0,PF,dt,x,y,Lx,Ly,chem)
%This function updates the phase field over dt
%semi-explicit Fourier spectral method is applied here
[m,n]=size(x);
dx=2*Lx/m;
dy=2*Ly/n;

Gp  =  36*(phi0-phi0.^2).*(1-2*phi0); 
[phi0x, phi0y] = grad(phi0,m, n, Lx, Ly);
phi_abs = sqrt(phi0x.^2 + phi0y.^2);

f=PF.Kprot*chem.*phi_abs/PF.tau-PF.lambda*Gp/(PF.tau*PF.epsilon^2); 
  
k1 = [0:m/2,  -(m/2-1):-1];
k2 = [0:n/2,  -(n/2-1):-1];
[k1, k2] = meshgrid(k1, k2);
k12  = k1.^2; 
k22  = k2.^2;
hx = (pi/Lx)^2; 
hy = (pi/Ly)^2;
p=dt*PF.lambda/PF.tau;     

f1=fft2(f);
phi1_old=fft2(phi0);
phi1_new=(f1*dt+phi1_old)./(1+p*(hx*k12+hy*k22));
phi=ifft2(phi1_new);

end
