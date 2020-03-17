function [strain_11, strain_22, strain_12] = fDIC_strain_method_3(nodeCoords_X,nodeCoords_Y,displacements_X,displacements_Y, num_row, num_col, Image_size)
% % four nodes finite analysis anlysis-bilinear   
[local_elem,local_elem_x,local_elem_y] = fDIC_shape_positions_4nodes(nodeCoords_X,nodeCoords_Y);
    [local_disp]                           = fDIC_shape_displacements_4nodes(displacements_X,displacements_Y);
    
    % Determine the missing points within the 3x3 grid in natural axes
    Integ_sep     = local_elem_x{1}(3)-local_elem_x{1}(1);
    Integ_coord_x(1) = -1;
    for n=2:Integ_sep+1
        Integ_coord_x(n) = Integ_coord_x(n-1)+2/(Integ_sep);
    end
    Integ_coord_x = Integ_coord_x';
    
    for M= 1:Integ_sep+1
        xi_X(M,:) = Integ_coord_x;
        xi_Y(:,M) = Integ_coord_x';
    end
    
    for I=1:(num_col-1)*(num_row-1) % the two rows and two colums from the edges of top right not be tested as they have been included in the 3x3 matrix
        parfor K=1:numel(xi_X)
            % SHAPE FUNCTIONS AND DERIVATES WRT LOCAL COORDS
            [N] = fDIC_shape_function_4nodes(xi_X(K),xi_Y(K));
            [dN_dXi] = fDIC_shape_function_deriviative_4nodes(xi_X(K),xi_Y(K));
            % COORDINATES OF INTEGRATION POINT
            location_x(I,K)  = int32(N*local_elem_x{I});
            location_y(I,K)  = int32(N*local_elem_y{I});
            % CALCULATE JACOBIAN MATRIX
            Jacobian = dN_dXi*local_elem{I};
            % CONVERT SFs TO GLOBAL COORDS
            dN_dxi   = Jacobian\dN_dXi;
            % Calculate strain matrix B (B=LN)
            [B] = fDIC_B_matrix_4nodes(dN_dxi);
            % CALCULATE  STRAIN
            strain = B*local_disp{I};
            strain_11_temp(I,K) =strain(1);
            strain_22_temp(I,K) =strain(2);
            strain_12_temp(I,K) =strain(3);
        end
    end
    % allocate memeory space
    strain_11  = zeros(Image_size);
    strain_22  = zeros(Image_size);
    strain_12  = zeros(Image_size);
    %sort out the results
    for i =1:(num_col-1)*(num_row-1)
        for j=1:numel(xi_X)
            strain_11(location_y(i,j),location_x(i,j)) = strain_11_temp(i,j);
            strain_22(location_y(i,j),location_x(i,j)) = strain_22_temp(i,j);
            strain_12(location_y(i,j),location_x(i,j)) = strain_12_temp(i,j);
        end
    end
end
