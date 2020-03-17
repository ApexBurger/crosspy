function [dN_dXi] = fDIC_shape_function_deriviative_9nodes(xi_X,xi_Y)
% determine shape function derivative for displacement gradient
dN_dXi(1,1) = 0.25*(1-xi_X)*xi_Y*(1-xi_Y);
dN_dXi(2,1) = 0.25*xi_X*(1-xi_X)*(1-xi_Y);

dN_dXi(1,2) = -0.25*(1+xi_X)*(1-xi_Y)*xi_Y;
dN_dXi(2,2) = -0.25*xi_X*(1+xi_X)*(1-xi_Y);

dN_dXi(1,3) = 0.25*(1+xi_X)*(1+xi_Y)*xi_Y;
dN_dXi(2,3) = 0.25*xi_X*(1+xi_X)*(1+xi_Y);

dN_dXi(1,4) = -0.25*(1-xi_X)*(1+xi_Y)*xi_Y;
dN_dXi(2,4) = -0.25*xi_X*(1-xi_X)*(1+xi_Y);


dN_dXi(1,5) =  -0.5* (-2*xi_X)*(1-xi_Y)*xi_Y;
dN_dXi(2,5) = -0.5* (1+xi_X)*(1-xi_X)*(1-xi_Y);

dN_dXi(1,6) =  0.5* (1+xi_X)*(1+xi_Y)*(1-xi_Y);
dN_dXi(2,6) =  0.5* xi_X*(1+xi_X)*(-2*xi_Y);

dN_dXi(1,7) =  0.5* (-2*xi_X)*(1+xi_Y)*xi_Y;
dN_dXi(2,7) =  0.5* (1+xi_X)*(1-xi_X)*(1+xi_Y);

dN_dXi(1,8) = -0.5*(1-xi_X)*(1-xi_Y)*(1+xi_Y);
dN_dXi(2,8) = -0.5*xi_X*(1-xi_X)*(-2*xi_Y);

dN_dXi(1,9) = (-2*xi_X)*(1-xi_Y)*(1+xi_Y);
dN_dXi(2,9) = (1-xi_X)*(1+xi_X)*(-2*xi_Y);


