 function [x_bp,y_bp]=bp_update(x_cont,y_cont,x_bp0,y_bp0,Vx_bp,Vy_bp,dt)
  % This function tries to apply the built in function 'polyxpoly' to move the marker points to their crossings with the contour of the updated phase field
  % I set a numerical minimum value for velocity here
  % Any membrane point with velocity below Velocity_num_minumun will stay where they are
  Velocity_num_minimum=1e-12;
  V_bp=sqrt(Vx_bp.*Vx_bp+Vy_bp.*Vy_bp);
%  if ~isempty(find(V_bp==0))
%      disp('Lazy point found! I can NOT go on -_-! \n');
%  end
  
  for iPoint=1:length(x_bp0)
      if V_bp(iPoint)>Velocity_num_minimum
     [x,y]=crossing_search(x_cont,y_cont,x_bp0(iPoint),y_bp0(iPoint),Vx_bp(iPoint),Vy_bp(iPoint),dt);
     x_bp(iPoint)=x;
     y_bp(iPoint)=y;
     %disp('Lazy point found!');
      else
     x_bp(iPoint)=x_bp0(iPoint);
     y_bp(iPoint)=y_bp0(iPoint);
      end
   end
 
 return
 end
