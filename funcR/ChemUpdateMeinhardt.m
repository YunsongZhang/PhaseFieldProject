function [a,b,c,signal]=ChemUpdateMeinhardt(x_bp,y_bp,x_bp0,y_bp0,Model,dt,a0,b0,c0)
%This function updates the Meinhardt model on the cell membrane
%A Crank-Nicolson finite difference scheme is applied here
%Since a typical time scale for one to see pseudopod pattern is much longer than the one for phase field evolution
%We scale the time here by a factor of timescale
%i.e. for each time interval "dt" for the phase field, we actully let the membrane chemical field go for "timescale*dt"

Da=Model.Da;
Db=Model.Db;
Dc=Model.Dc;
ra=Model.ra;
rb=Model.rb;
ba=Model.ba;
bc=Model.bc;
rc=Model.rc;
sc=Model.sc;   
sa=Model.sa;
dr=Model.dr;

m=length(x_bp0);
tmpx0=[x_bp0,x_bp0(1)];
tmpy0=[y_bp0,y_bp0(1)];
tmpx=[x_bp,x_bp(1)];
tmpy=[y_bp,y_bp(1)];
h_old=sqrt((tmpx0(2:end)-tmpx0(1:end-1)).^2+(tmpy0(2:end)-tmpy0(1:end-1)).^2);
h_new=sqrt((tmpx(2:end)-tmpx(1:end-1)).^2+(tmpy(2:end)-tmpy(1:end-1)).^2);
x_old=[0,cumsum(h_old)]; L_old=x_old(end);  x_old=x_old(1:end-1);
x_new=[0,cumsum(h_new)]; L_new=x_new(end);  x_new=x_new(1:end-1);

 %construct the Crank-Nicolson matrix for the reaction-diffusion
 %------------------------------
    %equation, applying the same method of "DiffNU.m", but we have three
    %components here:  A*u(n+1)=B*u(n)+C
    ivec=[1:m,1:m,1:m]; jvec=[1:m,2:m,1,m,1:(m-1)];
    hvec_new=[-1./(h_new.*circshift(h_new,[0,1]))];
    hvec_new=[hvec_new,1./(h_new.*(circshift(h_new,[0,1])+h_new))];
    hvec_new=[hvec_new,1./(circshift(h_new,[0,1]).*(circshift(h_new,[0,1])+h_new))];
    
    hvec_old=[-1./(h_old.*circshift(h_old,[0,1]))];
    hvec_old=[hvec_old,1./(h_old.*(circshift(h_old,[0,1])+h_old))];
    hvec_old=[hvec_old,1./(circshift(h_old,[0,1]).*(circshift(h_old,[0,1])+h_old))];
    
    delta_new=sparse(ivec,jvec,hvec_new);
    delta_old=sparse(ivec,jvec,hvec_old);
    
    
    
    Imatrix=sparse(1:m,1:m,ones(1,m));
    
    A1=(1+0.5*ra*dt)*Imatrix-Da*dt*delta_new;
    B1=(1-0.5*ra*dt)*Imatrix+Da*dt*delta_old;
    A2=(1+0.5*rb*dt)*Imatrix-Db*dt*delta_new;
    B2=(1-0.5*rb*dt)*Imatrix+Db*dt*delta_old;
    A3=(1+0.5*rc*dt)*Imatrix-Dc*dt*delta_new;
    B3=(1-0.5*rc*dt)*Imatrix+Dc*dt*delta_old;
   
 %  B21=sparse(1:m,1:m,0.5*rb*dt*ones(1,m));
   B31=sparse(1:m,1:m,0.5*bc*dt*ones(1,m));
 %  A21=sparse(1:m,1:m,-0.5*rb*dt*ones(1,m));
   A31=sparse(1:m,1:m,-0.5*bc*dt*ones(1,m));
   
   O=sparse(zeros(size(A1)));
   A=[A1,O,O;O,A2,O;A31,O,A3]; A=sparse(A);
   B=[B1,O,O;O,B2,O;B31,O,B3]; B=sparse(B);
   
   %Chemotaxis module
   %------------------------
   x_center=mean(x_bp0);
   y_center=mean(y_bp0);
   nx_gradient=cos(pi/2);
   ny_gradient=sin(pi/2);
   %tmp=(x_bp0-x_center)*nx_gradient+(y_bp0-y_center)*ny_gradient;
   Px=x_bp0-x_center; Py=y_bp0-y_center;
   tmp=Px*nx_gradient+Py*ny_gradient; tmp=tmp./sqrt(Px.*Px+Py.*Py);
   %tmp1=2*(tmp-min(tmp))/(max(tmp)-min(tmp))-1;
   %tmp1=0.02*tmp+ones(size(tmp));
   tmp1=0*tmp+ones(size(tmp));
  
   
    signal=ra*tmp1.*(1+dr*rand(size(x_bp0)));  
    C1=dt*signal.*(a0.*a0./b0+ba)./((sc+c0).*(1+sa*a0.*a0)); C1=sparse(C1');
    tmp=dt*rb*sum(0.5*h_old.*(a0+circshift(a0,[0,-1])))/L_old;
    C2=tmp*ones(m,1);  C2=sparse(C2);
    C=sparse([C1;C2;zeros(m,1)]);
    
   %end of construction of Crank-Nicolson matrix
   %-----------------------

 
    u_old=[a0';b0';c0'];
    
    u_new=A\(B*u_old+C);
    
    a=transpose(u_new(1:m));
    b=transpose(u_new(m+1:2*m));
    c=transpose(u_new(2*m+1:3*m));

end
