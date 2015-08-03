%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------- moose code----------------------
%-------------Yunsong Zhang---------------------
%------------2015-08--02-------------------------
%-----------------------------------------------
% This code reproduces Dr Wouter-Jan Rappel's code about a simple diffusion
% process on a 1D phase field interface
% I will compare the result with the simulation in general 1D case and also
% with analytic solution of this problem

% It is noticed from the code that the only difference between my previous code and Dr Rappel's code is the way we interpret the 2D v field as 1D distribution, I did sort of integral things, while Wouter only took 1 row of concentration!
% In later versions, we will see how such things will be done in a circular geometry or some even more complex geometry
clear
clc
close all

nx=200;
ny=20;
save_interval=500;
Du=1;
Dv=1;
x0=1.0;
dx=0.1;
dy=dx;
straal=5;
dt=1e-5;
Tmax=1.0;
gamma=20.0/straal;
thresh_phi=1e-4;

ufac=Du*dt/(dx*dx);
vfac=Dv*dt/(dx*dx);

x=[(-nx-1):1:(nx+1)]*dx;
y=[(-ny):1:ny]*dy;
[x,y]=meshgrid(x,y);

u=zeros(1,2*nx+3);
un=zeros(1,2*nx+3);
v=zeros(size(x));
vn=zeros(size(x));


% u=zeros(1,2*nx+3);
% un=zeros(1,2*nx+3);
% v=zeros(2*ny+3,2*nx+3);
% vn=zeros(2*ny+1,2*nx+1);

phitmp=tanh(-gamma*abs(y));
phi0=0.5+0.5*phitmp;
gphi0=0.25*gamma*gamma*(1-phitmp.*phitmp).^2;

box=(x==0).*(gphi0>thresh_phi);
v=v+x0*box;
u(nx+2)=x0;

gphi0x=0.5*(circshift(gphi0,[0,-1])+gphi0);
gphi0y=0.5*(circshift(gphi0,[-1,0])+gphi0);

% figure(1)
% imagesc(phi0)
% figure(2)
% imagesc(g_phi0)
% figure(3)
% plot(g_phi0(:,2))
% figure(4)
% plot(phi0(:,2))
% shg

u_rec=[];
v_rec=[];

for nIter=1:floor(Tmax/dt)
    timenow=nIter*dt;
    
    tmpright = gphi0x.*(circshift(v,[0,-1])-v);
    tmpleft  = circshift(gphi0x,[0,1]).*(v-circshift(v,[0,1]));
    tmpup    = gphi0y.*(circshift(v,[-1,0])-v);
    tmpdown  = circshift(gphi0y,[1,0]).*(v-circshift(v,[1,0]));
    vn=v+vfac*(tmpright-tmpleft+tmpup-tmpdown)./gphi0;
    
    
    un=u+ufac*(circshift(u,[0,-1])+circshift(u,[0,1])-2*u);
    
    
    v(gphi0>thresh_phi)=vn(gphi0>thresh_phi);
    u=un;
    
    if mod(nIter,save_interval)==0
       fprintf('t=%.3f\n',timenow)
       u_rec(end+1,:)=u;
       v_rec(end+1,:)=v(21,:);
    end
    
end

figure(1)
imagesc(u_rec')
colorbar
figure(2)
imagesc(v_rec')
colorbar
figure(3)
imagesc(u_rec'-v_rec')
colorbar
figure(4)
plot(u_rec(1,:),'o')
hold on
plot(v_rec(1,:),'-')
legend('u','v')
shg




