function [Vx_bp,Vy_bp,curvature]=velocity_curvature(phi0,phi,x_bp0,y_bp0,x,y,dx,dy,dt,PerimeterGrowthRate)

% This function is used to calculate the velocities and curvatures of all boundary marker points
% Here a spline interpolation is applied 
% We calculate the 1st derivatives and hessian of phi and phi0 at all boundary marker points
% The curvature is calculated as \nabla \dot (-\nabla phi/|\nabla phi|)

    %global eps
    eps=1e-4;

    PHI0=griddedInterpolant(x',y',phi0','spline');
    PHI=griddedInterpolant(x',y',phi','spline');
    
    % calculate values and 1st & 2nd derivatives of phi and phi0
    %---------------------------
    phi0_bp=PHI0(x_bp0,y_bp0);
    phi_bp=PHI(x_bp0,y_bp0);
    dphi0x_bp=(PHI0(x_bp0+eps,y_bp0)-PHI0(x_bp0-eps,y_bp0))/(2*eps);
    dphi0y_bp=(PHI0(x_bp0,y_bp0+eps)-PHI0(x_bp0,y_bp0-eps))/(2*eps);
    d2phi0x2_bp=(PHI0(x_bp0+eps,y_bp0)+PHI0(x_bp0-eps,y_bp0)-2*PHI0(x_bp0,y_bp0))/(eps*eps);
    d2phi0y2_bp=(PHI0(x_bp0,y_bp0+eps)+PHI0(x_bp0,y_bp0-eps)-2*PHI0(x_bp0,y_bp0))/(eps*eps);
    d2phi0xy_bp=(PHI0(x_bp0+eps,y_bp0+eps)+PHI0(x_bp0-eps,y_bp0-eps)-PHI0(x_bp0+eps,y_bp0-eps)-PHI0(x_bp0-eps,y_bp0+eps))/(4*eps*eps);
    absdphi0_bp=sqrt(dphi0x_bp.*dphi0x_bp+dphi0y_bp.*dphi0y_bp);
    %---------------------------

    %calculate the normal components of the velocities of boundary marker points
    %--------------------------
    Vn_bp=(phi_bp-phi0_bp)./(absdphi0_bp*dt);  
    Vn_bp(absdphi0_bp<eps)=0;
    nx_bp=-dphi0x_bp./absdphi0_bp;
    nx_bp(absdphi0_bp<eps)=0;
    ny_bp=-dphi0y_bp./absdphi0_bp;
    ny_bp(absdphi0_bp<eps)=0;
    Vnx_bp=Vn_bp.*nx_bp;
    Vny_bp=Vn_bp.*ny_bp;
    %-------------------------

    %calculate the curvature values at each of the boundary marker points
    %-------------------------
    % curvatures are calculated by \nabla \dot (-\nabla phi/ |\nabla phi|)
    % The formula can be referred in the manual for level set toolbox Page 102 
    % curvatureSecond, the only difference is a "minus" sign
    tmp1=2*dphi0x_bp.*dphi0y_bp.*d2phi0xy_bp;
    tmp2=-d2phi0x2_bp.*dphi0y_bp.^2-d2phi0y2_bp.*dphi0x_bp.^2;
    tmp3=absdphi0_bp.^3;
    temp=(tmp1+tmp2)./tmp3;
    temp(absdphi0_bp<eps)=0;
    curvature=temp;
    %------------------------
   
    %calculate the tangential components of the velocities of boundary marker points
    %-------------------------
    s0=sqrt((circshift(x_bp0,[0,-1])-x_bp0).^2+(circshift(y_bp0,[0,-1])-y_bp0).^2);
    tmp=0.5*(s0+circshift(s0,[0,1]));
    tmp1=-tmp.*curvature.*Vn_bp;
    Vt_bp=cumsum(tmp1)+cumsum(s0)*PerimeterGrowthRate/sum(s0);
    tx_bp=-ny_bp;
    ty_bp=nx_bp;
    Vtx_bp=Vt_bp.*tx_bp;
    Vty_bp=Vt_bp.*ty_bp;
    %-------------------------

    %total velocities
    %--------------------------
    Vx_bp=Vnx_bp+Vtx_bp;
    Vy_bp=Vny_bp+Vtx_bp;

end




    

    





    

     
