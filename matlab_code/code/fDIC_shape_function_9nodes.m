function [N] = fDIC_shape_function_9nodes(xi_X,xi_Y)
% Shpae function for 3x3 grid
N(1,1) =   0.25*xi_X*(1-xi_X)*xi_Y*(1-xi_Y);
N(1,2) =  -0.25*xi_X*(1+xi_X)*(1-xi_Y)*xi_Y;
N(1,3) =   0.25*xi_X*(1+xi_X)*(1+xi_Y)*xi_Y;
N(1,4) =  -0.25*xi_X*(1-xi_X)*(1+xi_Y)*xi_Y;

N(1,5) =  -0.5* (1+xi_X)*(1-xi_X)*(1-xi_Y)*xi_Y;
N(1,6) =   0.5* xi_X*(1+xi_X)*(1+xi_Y)*(1-xi_Y);
N(1,7) =   0.5* (1+xi_X)*(1-xi_X)*(1+xi_Y)*xi_Y;
N(1,8) =  -0.5*xi_X*(1-xi_X)*(1-xi_Y)*(1+xi_Y);
N(1,9) =  (1-xi_X)*(1+xi_X)*(1-xi_Y)*(1+xi_Y);





