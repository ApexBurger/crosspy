function  [strain_11, strain_22, strain_12] = fDIC_strain_method_4_corners(displacements_X,displacements_Y,Position_X,Position_Y, strain_11, strain_22, strain_12)
% use first order polynormial fitting for the 4 corners points
[num_row, num_col] = size(displacements_X);
% first corner
displacement_X_corner1_X  = [displacements_X(1,1),displacements_X(1,2)];
displacement_X_corner1_Y  = [displacements_Y(1,1),displacements_Y(1,2)];
displacement_Y_corner1_X  = [displacements_X(1,1),displacements_X(2,1)];
displacement_Y_corner1_Y  = [displacements_Y(1,1),displacements_Y(2,1)];
position_X_corner1      = [Position_X(1,1), Position_X(1,2)];
position_Y_corner1      = [Position_Y(1,1), Position_Y(2,1)];
% determine the first order of polynomial fit which is equavalent to
% strain
polycoef_x = polyfit(position_X_corner1,displacement_X_corner1_X,1);
polycoef_y = polyfit(position_Y_corner1,displacement_Y_corner1_Y,1);
polycoef_xy = polyfit(position_X_corner1,displacement_X_corner1_Y,1);
polycoef_yx = polyfit(position_Y_corner1,displacement_Y_corner1_X,1);
strain_11(1,1) = polycoef_x(1);
strain_22(1,1) = polycoef_y(1);
strain_12(1,1) = 0.5*(polycoef_xy(1,1)+polycoef_yx(1,1));


% second corner
displacement_X_corner2_X  = [displacements_X(1,num_col-1),displacements_X(1,num_col)];
displacement_Y_corner2_Y  = [displacements_Y(1,num_col),displacements_Y(2,num_col)];
displacement_X_corner2_Y  = [displacements_Y(1,num_col-1),displacements_Y(1,num_col)];
displacement_Y_corner2_X  = [displacements_X(1,num_col),displacements_X(2,num_col)];
position_X_corner2      = [Position_X(1,num_col-1), Position_X(1,num_col)];
position_Y_corner2      = [Position_Y(1,num_col), Position_Y(2,num_col)];
% determine the first order of polynomial fit which is equavalent to
% strain
polycoef_x = polyfit(position_X_corner2,displacement_X_corner2_X,1);
polycoef_y = polyfit(position_Y_corner2,displacement_Y_corner2_Y,1);
polycoef_xy = polyfit(position_X_corner2,displacement_X_corner2_Y,1);
polycoef_yx = polyfit(position_Y_corner2,displacement_Y_corner2_X,1);
strain_11(1,num_col) = polycoef_x(1);
strain_22(1,num_col) = polycoef_y(1);
strain_12(1,num_col) = 0.5*(polycoef_xy(1,1)+polycoef_yx(1,1));


% third corner
displacement_X_corner3_X  = [displacements_X(num_row,1),displacements_X(num_row,2)];
displacement_Y_corner3_Y  = [displacements_Y(num_row-1,1),displacements_Y(num_row,1)];
displacement_X_corner3_Y  = [displacements_Y(num_row,1),displacements_Y(num_row,2)];
displacement_Y_corner3_X  = [displacements_X(num_row-1,1),displacements_X(num_row,1)];
position_X_corner3      = [Position_X(num_row,1), Position_X(num_row,2)];
position_Y_corner3      = [Position_Y(num_row-1,1), Position_Y(num_row,1)];
% determine the first order of polynomial fit which is equavalent to
% strain
polycoef_x = polyfit(position_X_corner3,displacement_X_corner3_X,1);
polycoef_y = polyfit(position_Y_corner3,displacement_Y_corner3_Y,1);
polycoef_xy = polyfit(position_X_corner2,displacement_X_corner2_Y,1);
polycoef_yx = polyfit(position_Y_corner2,displacement_Y_corner2_X,1);
strain_11(num_row,1) = polycoef_x(1);
strain_22(num_row,1) = polycoef_y(1);
strain_12(num_row,1) = 0.5*(polycoef_xy(1,1)+polycoef_yx(1,1));

% fourth corner
displacement_X_corner4_X  = [displacements_X(num_row,num_col-1),displacements_X(num_row,num_col)];
displacement_Y_corner4_Y  = [displacements_Y(num_row-1,num_col),displacements_Y(num_row,num_col)];
displacement_X_corner4_Y  = [displacements_Y(num_row,num_col-1),displacements_Y(num_row,num_col)];
displacement_Y_corner4_X  = [displacements_X(num_row-1,num_col),displacements_X(num_row,num_col)];
position_X_corner4      = [Position_X(num_row,num_col-1), Position_X(num_row,num_col)];
position_Y_corner4      = [Position_Y(num_row-1,num_col), Position_Y(num_row,num_col)];
% determine the first order of polynomial fit which is equavalent to
% strain
polycoef_x = polyfit(position_X_corner4,displacement_X_corner4_X,1);
polycoef_y = polyfit(position_Y_corner4,displacement_Y_corner4_Y,1);
polycoef_xy = polyfit(position_X_corner2,displacement_X_corner2_Y,1);
polycoef_yx = polyfit(position_Y_corner2,displacement_Y_corner2_X,1);
strain_11(num_row,num_col) = polycoef_x(1);
strain_22(num_row,num_col) = polycoef_y(1);
strain_12(num_row,num_col) = 0.5*(polycoef_xy(1,1)+polycoef_yx(1,1));
