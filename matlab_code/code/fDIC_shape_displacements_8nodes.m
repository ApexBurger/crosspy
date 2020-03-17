 function [local_disp] = fDIC_shape_displacements_8nodes(displacements_X,displacements_Y) 

for i=1:size(displacements_X,1)-2
    for j=1:size(displacements_X,2)-2
         local_disp{i,j}   = [displacements_X(i,j);     displacements_Y(i,j);     displacements_X(i,j+2);   displacements_Y(i,j+2);...
                             displacements_X(i+2,j+2); displacements_Y(i+2,j+2); displacements_X(i+2,j);   displacements_Y(i+2,j);...
                             displacements_X(i,j+1);   displacements_Y(i,j+1);   displacements_X(i+1,j+2); displacements_Y(i+1,j+2);...
                             displacements_X(i+2,j+1); displacements_Y(i+2,j+1); displacements_X(i+1,j);   displacements_Y(i+1,j)];
    end 
end 
