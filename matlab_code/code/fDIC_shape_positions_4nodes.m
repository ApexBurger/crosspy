function [local_elem, local_elem_x, local_elem_y] = fDIC_shape_positions_4nodes(nodeCoords_X,nodeCoords_Y) 

for i=1:size(nodeCoords_X,1)-1
    for j=1:size(nodeCoords_X,2)-1
        local_elem{i,j}   = [nodeCoords_X(i,j)     nodeCoords_Y(i,j);     nodeCoords_X(i,j+1)   nodeCoords_Y(i,j+1);...
                             nodeCoords_X(i+1,j+1) nodeCoords_Y(i+1,j+1); nodeCoords_X(i+1,j)   nodeCoords_Y(i+1,j)];
        local_elem_x{i,j} = [nodeCoords_X(i,j);   nodeCoords_X(i,j+1);  nodeCoords_X(i+1,j+1); nodeCoords_X(i+1,j)];
                             
        local_elem_y{i,j} = [nodeCoords_Y(i,j);   nodeCoords_Y(i,j+1);   nodeCoords_Y(i+1,j+1); nodeCoords_Y(i+1,j)];
    end 
end 
