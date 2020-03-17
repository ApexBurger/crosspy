function [B] = fDIC_B_matrix_8nodes(dN_dxi)
% this function is designed to work out the B matrix       

B = [dN_dxi(1,1), 0, dN_dxi(1,2), 0, dN_dxi(1,3), 0, dN_dxi(1,4),0,dN_dxi(1,5), 0, dN_dxi(1,6), 0, dN_dxi(1,7), 0, dN_dxi(1,8),0;
    0, dN_dxi(2,1), 0, dN_dxi(2,2), 0, dN_dxi(2,3), 0, dN_dxi(2,4), 0, dN_dxi(2,5), 0, dN_dxi(2,6), 0, dN_dxi(2,7), 0, dN_dxi(2,8);
    dN_dxi(2,1), dN_dxi(1,1), dN_dxi(2,2), dN_dxi(1,2), dN_dxi(2,3),dN_dxi(1,3),dN_dxi(2,4), dN_dxi(1,4), dN_dxi(2,5), dN_dxi(1,5), dN_dxi(2,6), dN_dxi(1,6), dN_dxi(2,7),dN_dxi(1,7),...
    dN_dxi(2,8), dN_dxi(1,8)];
       