function [strain, rotation, strain_effective, F_output] = fDIC_StrainCalc(Shift_X, Shift_Y, Position_X, Position_Y, Image_size, strain_method, Image_num, DIC_method, Image_interval)
% Created by Jun Jiang
% 02/2014 Imperial college
% ---2D plane strain using quartic finite element rectangle 9 nodes approach is adopted in this function
% ***all points within the 9 points grid can be calculated.
% ***9 node shape function is used here with the 9th one at the central.
% *** 9 points are transformed into a natural space and using Jacobian
% matrix to return to the physical space
%----------------------------------------------------------------------------------
%Input results
%node/ROI positions
nodeCoords_X      = Position_X; % careful on index and xy axis coordinate
nodeCoords_Y      = Position_Y;
if DIC_method ==1
    Image_num_DIC     = Image_num-1;
else
    Image_num_DIC     = Image_num-1;
end
strain = cell(1,Image_num_DIC);
rotation = cell(1,Image_num_DIC);

for I=1:Image_num_DIC
    
    %node/ROI displacements
    displacements_X   = Shift_X(:,:,I);
    displacements_Y   = Shift_Y(:,:,I);
    %determine the size of map
    [num_row, num_col] = size(Position_X);
    % -----------------------------------------------
    % LOCAL ELEMENT DISPLACMENT ARRAY
    % LOOP OVER ALL ELEMENTS
    % LOCAL ELEMENT CO-ORDINATES AND DISPLCAMENTS
    if  strain_method ==1 % 9nodes
        [strain_11(:,:,I), strain_22(:,:,I), strain_12(:,:,I)] = fDIC_strain_method_1(nodeCoords_X,nodeCoords_Y,displacements_X,displacements_Y, num_row, num_col, Image_size);
        
    elseif strain_method ==2 % 8nodes
        [strain_11(:,:,I), strain_22(:,:,I), strain_12(:,:,I)] = fDIC_strain_method_2(nodeCoords_X,nodeCoords_Y,displacements_X,displacements_Y, num_row, num_col, Image_size);
    elseif strain_method ==3 % 4nodes
        [strain_11(:,:,I), strain_22(:,:,I), strain_12(:,:,I)] = fDIC_strain_method_3(nodeCoords_X,nodeCoords_Y,displacements_X,displacements_Y, num_row, num_col, Image_size);
    elseif strain_method ==4 % least square quatartic approach
        % initalised outputs
        warning off
%         strain_temp = NaN(size(displacements_X,1),size(displacements_X,1),3);
%         rotation_temp = strain_temp;
%         strain_effective_temp = zeros(size(displacements_X));
       
        strain_11_temp = zeros(size(displacements_X));
        strain_22_temp = strain_11_temp;
        strain_12_temp = strain_11_temp;
        % determine strain at the edges and corners with a mix of first and second order polyfit
        [strain_11_temp, strain_22_temp, strain_12_temp] = fDIC_strain_method_4_corners(displacements_X,displacements_Y,Position_X,Position_Y, strain_11_temp, strain_22_temp, strain_12_temp);
        [strain_11_temp, strain_22_temp, strain_12_temp] = fDIC_strain_method_4_edges(displacements_X,displacements_Y,Position_X,Position_Y, strain_11_temp, strain_22_temp, strain_12_temp);
        
        % deal with the bulk part using "plus" which may have difficulty to work out the shear part
        for i=2: num_row-1
            for j =2:num_col-1
                displacement_hor_3points_X  = [displacements_X(i,j-1),displacements_X(i,j),displacements_X(i,j+1)];
                displacement_hor_3points_Y  = [displacements_Y(i,j-1),displacements_Y(i,j),displacements_Y(i,j+1)];
                displacement_ver_3points_X  = [displacements_X(i-1,j),displacements_X(i,j),displacements_X(i+1,j)];
                displacement_ver_3points_Y  = [displacements_Y(i-1,j),displacements_Y(i,j),displacements_Y(i+1,j)];
                position_X_3points      = [Position_X(i,j-1), Position_X(i,j), Position_X(i,j+1)];
                position_Y_3points      = [Position_Y(i-1,j), Position_Y(i,j), Position_Y(i+1,j)];
                
                % determine the second order of polynomial fit and its derivative which is equavalent to
                % strain
                polycoef_x  = polyder(polyfit(position_X_3points,displacement_hor_3points_X,2));
                polycoef_y  = polyder(polyfit(position_Y_3points,displacement_ver_3points_Y,2));
                polycoef_xy = polyder(polyfit(position_X_3points,displacement_hor_3points_Y,2));
                polycoef_yx = polyder(polyfit(position_Y_3points,displacement_ver_3points_X,2));
                % generate deformation gradient tesnor "F" as the third axis is not accessible assuming dz is 0
                % F = du/dX +I                 
                du_dx = polyval(polycoef_x,Position_X(i,j));
                dv_dy = polyval(polycoef_y,Position_Y(i,j));
                du_dy = polyval(polycoef_xy,Position_X(i,j));
                dv_dx = polyval(polycoef_yx,Position_Y(i,j));
                
                F = [du_dx, du_dy 0; dv_dx, dv_dy 0; 0 0 -(du_dx+dv_dy)]+eye(3); % compatibility du/dx + dv/dy _dw/dz =0 and assume shear strain is small 
                F_temp(i,j,:,:) = F; 
                C = F'*F; %Green-Lagranage tensor 
                % determine the eigen vector V and eigen value Q 
                [Q, V] = eig(C);
                Q = Q'; % Matlab have to transpose this Q to be conventional one 
                Ut     = sqrt(V); % square root of transformed C
                U      = Q'*Ut*Q; % transform back to orginal axes 
                U_1    = inv(U); % inverse of U 
                R      = F*U_1;   
                % determine Green strain tesnor and rotation tensor from F by symmetrical and anti symmetrical parts 
                strain_temp(i,j,:,:)   = 1/2*(F'*F-eye(3));
                rotation_temp(i,j,:,:) = R;

                % determine the effective strain to convert strain tensor
                % as a scalar p = (2/3 x strain:strain)^1/2 
                strain_effective_temp(i,j) = sqrt(2/3*trace(strain_temp(i,j)*strain_temp(i,j))); 
            end
        end
        strain{1,I}             = strain_temp;
        rotation{1,I}           = rotation_temp; 
        strain_effective(:,:,I) = strain_effective_temp; 
        F_output{1,I}           = F_temp; 
    end
end


