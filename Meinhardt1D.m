%%%%%%%%%%%%%%%%%
%---Moose code--------
%---Yunsong Zhang------
%----2015-7-30---------
%%%%%%%%%%%%%%%%%
%This code solves the Meinhardt Model in 1D 

clear
close all
clc

set(0,'DefaultFigureVisible','off')
%rng(1467)
rng('shuffle')

tic
% Parameters
dt=1;
maxTime=80000;
interval=400;
numPts=200;
Radius=0.1;
L0=2*pi*Radius;
saveevery=floor(interval/dt);
image_directory='./result';

Da=4e-7; %2e-4;
Db=4e-5;  %4
Dc=2.8*Da; %1.8*Da
ra=0.02;
rb=0.03;
ba=0.1;
bc=0.005;
rc=0.013;
sc=0.2;      
sa=0.0005;
dr=0.01; 

a_rec=[];
b_rec=[];
c_rec=[];

%Initialize the nonuniform meshgrid
x=linspace(0,L0,numPts+1); x=x(1:end-1); dx=x(2)-x(1);
x0=36; %0.5*L0;


% Initialize the chemical system
a0=zeros(numPts,1);
b0=0.1*ones(numPts,1);
c0=zeros(numPts,1);



%According to Crank-Nicolson scheme, the set of PDEs turn to 3 linear
%matrix equations:
%  A1*a(n+1)=B1*a(n)+C1
%  A2*b(n+1)=B2*b(n)+C2
%  A3*c(n+1)=B3*c(n)+C3

% Finally, all 3 equations go to one A*u(n+1)=B*u(n)+C

% construct A1,A2,A3,B1,B2,B3
% while C1,C2change over time, C3=0 forever
   w1=Da*dt/(2*dx*dx);
   w2=Db*dt/(2*dx*dx);
   w3=Dc*dt/(2*dx*dx);
   
   ivec=[1:numPts,1:numPts,1:numPts];
   jvec=[1:numPts,2:numPts,1,numPts,1:(numPts-1)];
   A1_vec=[(1+ra*dt/2+2*w1)*ones(1,numPts),-w1*ones(1,2*numPts)];
   B1_vec=[(1-ra*dt/2-2*w1)*ones(1,numPts),w1*ones(1,2*numPts)];
   A2_vec=[(1+rb*dt/2+2*w2)*ones(1,numPts),-w2*ones(1,2*numPts)];
   B2_vec=[(1-rb*dt/2-2*w2)*ones(1,numPts),w2*ones(1,2*numPts)];
   A3_vec=[(1+rc*dt/2+2*w3)*ones(1,numPts),-w3*ones(1,2*numPts)];
   B3_vec=[(1-rc*dt/2-2*w3)*ones(1,numPts),w3*ones(1,2*numPts)];
   
   A1=sparse(ivec,jvec,A1_vec);
   B1=sparse(ivec,jvec,B1_vec);
   A2=sparse(ivec,jvec,A2_vec);
   B2=sparse(ivec,jvec,B2_vec);
   A3=sparse(ivec,jvec,A3_vec);
   B3=sparse(ivec,jvec,B3_vec);
   
    %  B21=sparse(1:numPts,1:numPts,0.5*rb*dt*ones(1,numPts));
   B31=sparse(1:numPts,1:numPts,0.5*bc*dt*ones(1,numPts));
 %  A21=sparse(1:numPts,1:numPts,-0.5*rb*dt*ones(1,numPts));
   A31=sparse(1:numPts,1:numPts,-0.5*bc*dt*ones(1,numPts));
   
   C3=sparse(zeros(numPts,1));
   
   
   O=sparse(zeros(size(A1)));
   A=[A1,O,O;O,A2,O;A31,O,A3]; A=sparse(A);
   B=[B1,O,O;O,B2,O;B31,O,B3]; B=sparse(B);

for k=1:floor(maxTime/dt)
   
   time=k*dt;
   
   %construct the C vector in the Crank-Nicolson scheme
    signal=ra*(1+0.02*cos(2*pi*(x-x0)/L0)).*(1+dr*rand(size(x)));  
    C1=dt*signal'.*(a0.*a0./b0+ba)./((sc+c0).*(1+sa*a0.*a0)); C1=sparse(C1);
    tmp=dt*rb*sum(a0)*dx/L0;
    C2=tmp*ones(numPts,1);  C2=sparse(C2);
    C=sparse([C1;C2;C3]);
   
   
   u_old=[a0;b0;c0];
   u_new=A\(B*u_old+C);
   
   a=u_new(1:numPts);
   b=u_new((numPts+1):(2*numPts));
   c=u_new((2*numPts+1):(3*numPts));
   
   a0=a;
   b0=b;
   c0=c;  
   
   %---record & check module----------
    if mod(k,saveevery)==0
         fprintf('time=%.1f\n',k*dt);
         a_rec(:,end+1)=a0;
         b_rec(:,end+1)=b0;
         c_rec(:,end+1)=c0;
               
        
     end
    
     if max(max(isnan(u_new)))==1
         fprintf('t=%.3f',k*dt)
         error('Something is wrong!')
     end
    %---end of record & check module--------- 
    
    
end

toc

handle=imagesc((interval:interval:maxTime)/100,x,a_rec);
title('spatio-temporal pattern of a','FontSize',40);
colorbar
shading interp
saveas(handle,'./result_a','jpg');
shg

save('./1D.mat')