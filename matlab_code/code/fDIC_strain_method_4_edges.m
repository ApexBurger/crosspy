
function  [strain_11, strain_22, strain_12] = fDIC_strain_method_4_edges(displacements_X,displacements_Y,Position_X,Position_Y, strain_11, strain_22, strain_12)
% treat four edges differently from the bulk part.
% Second order polynomial is used for one axis on the edge while the first order is used for the second axis
[num_row, num_col] = size(displacements_X);
% top edge 
for i=1
    for j=2:num_col-1
        displacement_X_3points_X  = [displacements_X(i,j-1),displacements_X(i,j),displacements_X(i,j+1)];
        displacement_Y_2points_Y  = [displacements_Y(i,j),displacements_Y(i+1,j)];
        displacement_X_3points_Y  = [displacements_Y(i,j-1),displacements_Y(i,j),displacements_Y(i,j+1)];
        displacement_Y_2points_X  = [displacements_X(i,j),displacements_X(i+1,j)];
        position_X_3points      = [Position_X(i,j-1), Position_X(i,j), Position_X(i,j+1)];
        position_Y_2points      = [Position_Y(i,j), Position_Y(i+1,j)];
        
        % determine the first order of polynomial fit which is equavalent to
        % strain
        polycoef_x  = polyder(polyfit(position_X_3points,displacement_X_3points_X,2));
        polycoef_y  = polyfit(position_Y_2points,displacement_Y_2points_Y,1);
        polycoef_xy = polyder(polyfit(position_X_3points,displacement_X_3points_Y,2));
        polycoef_yx = polyfit(position_Y_2points,displacement_Y_2points_X,1);
        strain_11(i,j) = polyval(polycoef_x,Position_X(i,j));
        strain_22(i,j) = polycoef_y(1,1);
        strain_12(i,j) = 0.5*(polyval(polycoef_xy,Position_X(i,j))+polycoef_yx(1,1));
    end
end
% bottom edge 
for i=num_row
    for j=2:num_col-1
        displacement_X_3points_X  = [displacements_X(i,j-1),displacements_X(i,j),displacements_X(i,j+1)];
        displacement_Y_2points_Y  = [displacements_Y(i-1,j),displacements_Y(i,j)];
        displacement_X_3points_Y  = [displacements_Y(i,j-1),displacements_Y(i,j),displacements_Y(i,j+1)];
        displacement_Y_2points_X  = [displacements_X(i-1,j),displacements_X(i,j)];
        position_X_3points      = [Position_X(i,j-1), Position_X(i,j), Position_X(i,j+1)];
        position_Y_2points      = [Position_Y(i-1,j), Position_Y(i,j)];
        
        % determine the first order of polynomial fit which is equavalent to
        % strain
        polycoef_x = polyder(polyfit(position_X_3points,displacement_X_3points_X,2));
        polycoef_y = polyfit(position_Y_2points,displacement_Y_2points_Y,1);
        polycoef_xy = polyder(polyfit(position_X_3points,displacement_X_3points_Y,2));
        polycoef_yx = polyfit(position_Y_2points,displacement_Y_2points_X,1);
        strain_11(i,j) = polyval(polycoef_x,Position_X(i,j));
        strain_22(i,j) = polycoef_y(1);
        strain_12(i,j) = 0.5*(polyval(polycoef_xy,Position_X(i,j))+polycoef_yx(1));
    end
end
%left side edge 
for i=2:num_row-1
    for j=1
        displacement_X_2points_X  = [displacements_X(i,j),displacements_X(i,j+1)];
        displacement_Y_3points_Y  = [displacements_Y(i-1,j),displacements_Y(i,j),displacements_Y(i+1,j)];
        displacement_X_2points_Y  = [displacements_Y(i,j),displacements_Y(i,j+1)];
        displacement_Y_3points_X  = [displacements_X(i-1,j),displacements_X(i,j),displacements_X(i+1,j)];
        position_X_2points      = [Position_X(i,j), Position_X(i,j+1)];
        position_Y_3points      = [Position_Y(i-1,j), Position_Y(i,j), Position_Y(i+1,j)];
        
        % determine the first order of polynomial fit which is equavalent to
        % strain
        polycoef_x = polyfit(position_X_2points,displacement_X_2points_X,1);
        polycoef_y = polyder(polyfit(position_Y_3points,displacement_Y_3points_Y,2));
        polycoef_xy = polyfit(position_X_2points,displacement_X_2points_Y,1);
        polycoef_yx = polyder(polyfit(position_Y_3points,displacement_Y_3points_X,2));
        strain_11(i,j) = polycoef_x(1);
        strain_22(i,j) = polyval(polycoef_y,Position_Y(i,j));
        strain_12(i,j) = 0.5*( polycoef_xy(1)+polyval(polycoef_yx,Position_Y(i,j)));
    end
end
% right side edge
for i=2: num_row-1
    for j =num_col
        displacement_X_2points_X  = [displacements_X(i,j-1),displacements_X(i,j)];
        displacement_Y_3points_Y  = [displacements_Y(i-1,j),displacements_Y(i,j),displacements_Y(i+1,j)];
        displacement_X_2points_Y  = [displacements_Y(i,j-1),displacements_Y(i,j)];
        displacement_Y_3points_X  = [displacements_X(i-1,j),displacements_X(i,j),displacements_X(i+1,j)];
        position_X_2points      = [Position_X(i,j-1), Position_X(i,j)];
        position_Y_3points      = [Position_Y(i-1,j), Position_Y(i,j), Position_Y(i+1,j)];
        
        % determine the first order of polynomial fit which is equavalent to
        % strain
        polycoef_x = polyfit(position_X_2points,displacement_X_2points_X,1);
        polycoef_y = polyder(polyfit(position_Y_3points,displacement_Y_3points_Y,2));
        polycoef_xy = polyfit(position_X_2points,displacement_X_2points_Y,1);
        polycoef_yx = polyder(polyfit(position_Y_3points,displacement_Y_3points_X,2));
        strain_11(i,j) = polycoef_x(1);
        strain_22(i,j) = polyval(polycoef_y,Position_Y(i,j));
        strain_12(i,j) = 0.5*(polycoef_xy(1)+polyval(polycoef_yx,Position_Y(i,j)));
    end
end