function [dN_dXi] = fDIC_shape_function_deriviative_4nodes(xi_X,xi_Y)
% determine shape function derivative for displacement gradient
dN_dXi(1,1) =  -0.25*(1-xi_Y);
dN_dXi(2,1) =  -0.25*(1-xi_X);

dN_dXi(1,2) =  0.25*(1-xi_Y);
dN_dXi(2,2) = -0.25*(1+xi_X);

dN_dXi(1,3) =  0.25*(1+xi_Y);
dN_dXi(2,3) =  0.25*(1+xi_X);

dN_dXi(1,4) = -0.25*(1+xi_Y);
dN_dXi(2,4) =  0.25*(1-xi_X);



