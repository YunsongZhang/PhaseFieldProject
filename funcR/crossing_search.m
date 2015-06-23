function [x_cross,y_cross]=crossing_search(x_cont,y_cont,x0,y0,vx,vy,dt)
% This function is to figure out the crossing of a given membrane marker point to the contour of the updated phase field in its velocity direction
% In case no crossing is found, the membrane marker point stays where it is 
% In case of more than 1 crossing,the point goes to the most likely one according to its velocity
 factor1=-50;
 factor2=50;
 xline=[x0+vx*dt*factor1,x0+vx*dt*factor2];
 yline=[y0+vy*dt*factor1,y0+vy*dt*factor2];
 [xc,yc]=polyxpoly(x_cont,y_cont,xline,yline);
 
 if length(xc)==1
     x_cross=xc;
     y_cross=yc;
 elseif isempty(xc)
     %disp('no crossing!\n');
     x_cross=x0;
     y_cross=y0;
 elseif length(xc)>1
     xp=x0+vx*dt; yp=y0+vy*dt;
     error_compare=sqrt((xc-xp).^2+(yc-yp).^2);
     %error_compare=abs(sqrt((xc-x0).^2+(yc-y0).^2)/sqrt(vx*vx+vy*vy)-dt);
     %error_compare=[xc-x0,yc-y0]./[vx,vy];
     min_error=min(error_compare);
     x_c=xc(find(error_compare==min_error)); x_cross=x_c(1);
     y_c=yc(find(error_compare==min_error)); y_cross=y_c(1);
 end
 
 return
    



end
