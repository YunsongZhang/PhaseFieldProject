%-----------------------------------
% 05-10-2015
% Yunsong Zhang
%----------------------------
% This code tries to reproduce Robert Insall's SIAM paper in our own framework
% We take parameters and technical details from Robert Insall's SIAM paper as much as possible
% However,it turns out impossible for us to follow every parameter of Insall's paper to get their result
%
% Notice here we put all functions in the folder ./funcR
% Attention! You need to compile the file NDirac2.c by the command "mex NDirac2.c" before running this code!
% And also many functions and other .m files have names ending in "R", which is short for "reproduce"
%-------------------------------------

% Notice: Each time when you run this script,please write down a very
% brief description of the simulation into the variable simu_descrip!



clear
close all
clc


tic
set(0,'DefaultFigureVisible','off')
addpath('./funcR');   % Here we put all the functions in the folder funcR

%------------
%random number generator setting
rng(1467)    % 1467 is the plate # of my car -_-
%rng('shuffle')
%------------

%% ----Parameter Setting---------------------------------

%--Phase Field Parameters--
PF=struct();
PF.lambda0=2e-6;     %cortical tension factor
PF.Kprot=0;  % 1e-5;
PF.beta=1e-6;
PF.r0=0.1;           % initial radius of the round Phase Field cell
PF.epsilon=0.01;     % the width coefficient of the initial Phase Field cell
PF.tau=1.0;
%---------

%--Meinhardt Model Parameters--
Model=struct();
Model.Da=4e-7;               %Diffusion constant 
Model.Db=4e-5;
Model.Dc=2.5*Model.Da;       % The ratio of Dc/Da determines the pattern: shall be larger than 1 for bifurcations, and the larger, the higher frequency for bifurcations
Model.ra=0.02;   
Model.rb=0.03;
Model.ba=0.1;
Model.bc=0.005;
Model.rc=0.013;
Model.sc=0.2;      
Model.sa=0.0005;
Model.dr=0.01;               % noise magnitude
%---------

%--Technical Parameters for simulating--
Lx=0.5;
Ly=0.5;
m=2^9;
n=2^9;
M=200; 
dt=5e-3;
timescale=200; %To scale the time of membrane reactions and the mechanics
time_interval=0.5;
saveevery=floor(time_interval/dt);
maxTime=400;
global eps
eps=1e-5;  %minimum error value below which any value is considered to be 0
pic=0;
%---------

%--Other non-independent parameters--
%PF.A0=pi*PF.r0^2;   %Initial area of the cell: should stay almost constant over the whole 2D simulation
%---------

%--Saving directory Setting--
simu_descrip='uniformly shrinking';
timeofstart=clock;
save_directory=['./',date,'/',simu_descrip];
savename=sprintf('%s/data.mat',save_directory);
saveevery=floor(time_interval/dt);
%---------



%-------End of Parameter Setting-------------------------------------
if isempty(simu_descrip)
    disp('Oops! You forgot to describe this simulation!');
    break;
elseif ~exist(save_directory)
    mkdir(save_directory)
end
% Do NOT forget to write some description for the simulation



%% Initialize the system

% ---Phase Field Initialization---
x=linspace(-Lx,Lx,m+1); x=x(1:end-1); dx=x(2)-x(1);
y=linspace(-Ly,Ly,n+1); y=y(1:end-1); dy=y(2)-y(1);
[x,y]=meshgrid(x,y);
r=sqrt(x.*x+y.*y);
phi0=0.5+0.5*tanh((PF.r0-r)/PF.epsilon);
PF.A0=sum(sum(phi0))*dx*dy;   Area_old=PF.A0;
PF.lambda=PF.lambda0;
d1_Area=0;                %first derivative over time of the cell area
%-------

%---Boundary point Initialization---
theta0=linspace(0,2*pi,M+1);  theta0=theta0(1:end-1);
x_bp0=PF.r0*cos(theta0);
y_bp0=PF.r0*sin(theta0);  % x_bp0 and y_bp0 record the positions of all boundary marker points
s0=sqrt((circshift(x_bp0,[0,-1])-x_bp0).^2+(circshift(y_bp0,[0,-1])-y_bp0).^2);
Perimeter0=sum(s0); PerimeterGrowthRate=0;
%------

%Initialize the chemical field on the membrane
a0=zeros(size(x_bp0));
b0=0.1*ones(size(x_bp0));  %Here the chemical field of b is assigned any very small value to avoid singularity
c0=zeros(size(x_bp0));
%------

% Initialize all tools for recording
a_rec=[];
b_rec=[];
c_rec=[];
contour_rec={};
chem_rec=[];
xbp_rec=x_bp0;
ybp_rec=y_bp0;
curvature_rec=[];
Perimeter_rec=[];
Velocity_rec=[];
CellArea_rec=[];
%------

%% Start moving the cell!
for kSteps=1:floor(maxTime/dt)
    
        time=kSteps*dt;
        
   % step 1     
   % Put membrane chemical field into 2D phase field framework
   %-------   
    chem=NDirac2(x_bp0,y_bp0,a0,x,y);   %immersed boundary condition applied
   % chem=zeros(size(x));
   %-------
   
   % step 2
   % Update the phase field
   %-------
    PF.lambda=PF.lambda0;
    phi=phi_updateR(phi0,PF,dt,x,y,Lx,Ly,chem);
    %phi=0.5+0.5*tanh((PF.r0+0.0001-r)/PF.epsilon);
   %-------
   
   % step 3
   % Update the cortial tension factor PF.lambda
   %-------
     Area_new=sum(sum(phi))*dx*dy;
     d1_Area=(Area_new-Area_old)/dt;
     d1_lambda=PF.lambda*PF.lambda0*(Area_new-Area_old+d1_Area)/((PF.lambda+PF.lambda)*PF.A0)-PF.beta*PF.lambda;
     PF.lambda=PF.lambda+d1_lambda*dt;
     Area_old=Area_new;
   %------
   
   % step 4
   % Figure out the contour information of the phase field
   %-------
    [Contour,~]=contour(x,y,phi,[0.5,0.5]); 
    membrane=Contour(:,2:end);
    x_cont=membrane(1,:);
    y_cont=membrane(2,:);
   %-------
    
   % step 5
   % Figure out the velocity,curvatures of all membrane marker points
   %------
     [Vx_bp,Vy_bp,curvature]=velocity_curvature(phi0,phi,x_bp0,y_bp0,x,y,dx,dy,dt,PerimeterGrowthRate);
     Velocity_bp=sqrt(Vx_bp.*Vx_bp+Vy_bp.*Vy_bp);
   %------
     
   % step 6  
   %------
   % move the boundary marker points along their velocity directions to their crossings with the contour of the phase field cell
    [x_bp,y_bp]=bp_update(x_cont,y_cont,x_bp0,y_bp0,Vx_bp,Vy_bp,dt);
    s=sqrt((circshift(x_bp,[0,-1])-x_bp).^2+(circshift(y_bp,[0,-1])-y_bp).^2);
    Perimeter=sum(s);
    PerimeterGrowthRate=(Perimeter-Perimeter0)/dt;
   %-------  
  
   % step 7
   %-------
   % solve the reaction-diffusion equation on the membrane
   % Crank-Nicolson finite difference scheme applied on a nonuniform and changing 1D meshgrid
   [a,b,c]=ChemUpdateMeinhardt(x_bp,y_bp,x_bp0,y_bp0,Model,dt*timescale,a0,b0,c0);
   %------- 

   %-------
   % Values updating!
    phi0=phi;
    Perimeter0=Perimeter;
    x_bp0=x_bp;
    y_bp0=y_bp;
    a0=a;
    b0=b;
    c0=c;
   %------
    

   %------
   % Module for saving and recording datas
    if mod(kSteps,saveevery)==0
            fprintf('time=%1.3f\n',kSteps*dt);
            handle=imagesc([min(x(:)),max(x(:))],[min(y(:)),max(y(:))],phi0,[0,1]);
           % handle=imagesc([min(x(:)),max(x(:))],[min(y(:)),max(y(:))],chem);
            colorbar
            pic=pic+1;
            picname=sprintf('%s/figure%d',save_directory,pic);
            saveas(handle,picname,'jpg');
            
             a_rec(:,end+1)=a0';
             b_rec(:,end+1)=b0';
             c_rec(:,end+1)=c0';
             
             contour_rec{end+1}=membrane;
             chem_rec(:,:,end+1)=chem;
             xbp_rec(end+1,:)=x_bp0;
             ybp_rec(end+1,:)=y_bp0;
             curvature_rec(end+1,:)=curvature;
             Perimeter_rec(end+1)=Perimeter0;
             
             Velocity_bp=sqrt(Vx_bp.*Vx_bp+Vy_bp.*Vy_bp);
             Velocity_rec(end+1,:)=Velocity_bp;
             CellArea_rec(end+1)=sum(sum(phi0))*dx*dy;
             
             if mod(kSteps,10*saveevery)==0
             save(savename);
             end
            
            if max(phi0(:))>1.05 || min(phi0(:))<-0.08
                error('phi value is not OK');
            end
    end
     %-----end of module for saving and recording datas 
 end
% end of cell movement 
 


AreaChangeRatio=sum(sum(phi0))*dx*dy/PF.A0 



% plot and outputs starts:
handle=imagesc(a_rec);
colorbar
shading interp
picname=sprintf('%s/a_rec',save_directory);
saveas(handle,picname,'jpg');

handle=imagesc(b_rec);
colorbar
shading interp
picname=sprintf('%s/b_rec',save_directory);
saveas(handle,picname,'jpg');

handle=imagesc(c_rec);
colorbar
shading interp
picname=sprintf('%s/c_rec',save_directory);
saveas(handle,picname,'jpg')

moviename=sprintf('%s/movie.gif',save_directory);

%save(savename);
quiver(x_bp0,y_bp0,Vx_bp,Vy_bp)
for k=1:floor(maxTime/time_interval)
    picname=sprintf('%s/figure%d.jpg',save_directory,k);
    I=imread(picname);
    [I,map]=rgb2ind(I,256);
    if k==1
        imwrite(I,map,moviename,'gif','Loopcount',1,'DelayTime',0.05);
    else
        imwrite(I,map,moviename,'gif','WriteMode','append','DelayTime',0.05);
    end
end


total_time=toc/3600;
fprintf('The code has used %.3f hours\n',total_time);



Velocity_rec=sqrt(((xbp_rec(2:end,:)-xbp_rec(1:end-1,:)).^2+(xbp_rec(2:end,:)-xbp_rec(1:end-1,:)).^2))/time_interval;


 handle=imagesc(Velocity_rec');
colorbar
shading interp
picname=sprintf('%s/Velocity_rec',save_directory);
saveas(handle,picname,'jpg')
 
save(savename);

