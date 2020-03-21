function [ resids ] = rotationobj( params,x_shift,y_shift,x_pos,y_pos)
%ROTATIONOBJ Summary of this function goes here
%   Detailed explanation goes here

xc=params(1);
yc=params(2);
theta=params(3);

%calculate a rotation matrix
rot_mat=[cosd(theta),-sind(theta);sind(theta),cosd(theta)];

%rotate the existing grid about (xc,yc)
new_points=[x_pos(:)-xc,y_pos(:)-yc]*rot_mat;

%calculate the shifts
shifts = [x_pos(:)-xc-new_points(:,1),y_pos(:)-yc-new_points(:,2)];

%find the residuals
resids=abs([shifts(:,1)-x_shift;shifts(:,2)-y_shift]);
 
end

