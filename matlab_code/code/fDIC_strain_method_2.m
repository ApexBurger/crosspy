function [strain_11, strain_22, strain_12] = fDIC_strain_method_2(nodeCoords_X,nodeCoords_Y,displacements_X,displacements_Y, num_row, num_col, Image_size)
% 8 nodes Finite element analysis for strain 

[local_elem] = fDIC_shape_positions_8nodes(nodeCoords_X,nodeCoords_Y);
[local_disp] = fDIC_shape_displacements_8nodes(displacements_X,displacements_Y);

strain_11  = zeros(size(displacements_X));
strain_22  = strain_11;
strain_12  = strain_11;
for I=1:(num_row-2)
    for J =1:(num_col-2) % the two rows and two colums from the edges of top right not be tested as they have been included in the 3x3 matrix
        % SHAPE FUNCTIONS AND DERIVATES WRT LOCAL COORDS
        [N] = fDIC_shape_function_8nodes(-1,-1);
        [dN_dXi] = fDIC_shape_function_deriviative_8nodes(-1,-1);
        % CALCULATE JACOBIAN MATRIX
        Jacobian = dN_dXi*local_elem{I,J};
        % CONVERT SFs TO GLOBAL COORDS
        dN_dxi   = Jacobian\dN_dXi;
        % Calculate strain matrix B (B=LN)
        [B] = fDIC_B_matrix_8nodes(dN_dxi);
        % CALCULATE  STRAIN
        strain = B*local_disp{I,J};
        strain_11(I,J) =strain(1);
        strain_22(I,J) =strain(2);
        strain_12(I,J) =strain(3);
    end
end 
