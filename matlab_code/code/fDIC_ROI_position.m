
function [ROI_position_X, ROI_position_Y,ROI_num_x,ROI_num_y,ROI_coordinator,ROI_num] = fDIC_ROI_position(ROI_size,ROI_overlap,boundary)
% POI position is calculated within the function

ROI_posi_x =  ROI_size/2+boundary(1):(floor(ROI_size*(1-ROI_overlap))):(boundary(3)+boundary(1)-ROI_size/2);
ROI_posi_y = ROI_size/2+boundary(2):(floor(ROI_size*(1-ROI_overlap))):(boundary(4)+boundary(2)-ROI_size/2);
for n=1:5
    if ROI_posi_x(end)+ROI_size/2 >boundary(3)+boundary(1)
        ROI_posi_x(end) =[];
    end
    if ROI_posi_y(end)+ROI_size/2 >boundary(4)+boundary(2)
        ROI_posi_y(end) =[];
    end
    if ROI_posi_y(end)+ROI_size/2 >boundary(4)+boundary(2)
        ROI_posi_y(end) =[];
    end
end
[ROI_position_X, ROI_position_Y] = meshgrid(ROI_posi_x, ROI_posi_y);
ROI_num_x = numel(ROI_posi_x);
ROI_num_y = numel(ROI_posi_y);
% set up the cooridnates of the centre of ROIs, convert x-y axis from image/xy to matrix index
for i=1: size(ROI_position_X,1)
    for j=1:size(ROI_position_Y,2)
        ROI_coordinator{i,j} =[ROI_position_X(i,j) ROI_position_Y(i,j)];
    end
end
ROI_num = numel(ROI_coordinator);