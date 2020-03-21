function [Image_re] = fDIC_Image_correction(Image_folder, Image_list, Shift_X, Shift_Y, Position_X, Position_Y)
% this funciton is designed by Ben Britton to correct the translation as
% well as rotation of the SEM images in order to allow smaller ROI to be used.
% ------05/2014
Image_1 = im2double(imread([Image_folder,Image_list{1}]));
[pix_xgrid,pix_ygrid]=meshgrid(1:size(Image_1,2),1:size(Image_1,1));
Image_re = zeros(size(Image_1,1),size(Image_1,2),numel(Image_list)-1);
Image_re(:,:,1) = Image_1; 
for I=1:numel(Image_list)-1
    
    Image = im2double(imread([Image_folder,Image_list{I+1}])); 
    X_Shifts=Shift_X(:,:,I);
    Y_Shifts=Shift_Y(:,:,I);
    % Correct the rigid translation in x and y axes
    X_Shifts_new = X_Shifts(:)-mean(X_Shifts(:));
    Y_Shifts_new = Y_Shifts(:)-mean(Y_Shifts(:));
    % first guess of the rotation origin and rotation angle
    params0=[size(Image,2)/8,size(Image,1)/1,-2];
    %define the bounds and starting guess for our rotation correction
    %params0 = initial guess
    %UB = upper bound
    %LB = lower bound
    %given as [xc,yc,angle] in [pixels,pixels,degrees]
    
    UB=[size(Image,2),size(Image,1),45];
    LB=[1,1,-45];
    [ params1 ] = rotationcalc(X_Shifts_new(:),Y_Shifts_new(:),Position_X(:),Position_Y(:),params0,UB,LB );
    
    %now fit a model to find out how to have the minimum global rotation for
    %the shifts
    %free parameters are the centre coordinate of the rotation (Xc,Yc) and the
    %rotation angle (theta)
    
    % correct the image rotation
    Centre_X = params1(1);
    Centre_Y = params1(2);
    Theta    = params1(3);
    Rotation = [cosd(Theta), -sind(Theta); sind(Theta), cosd(Theta)];
    
    % apply correction to the image
    
    Image_shifts    = [pix_xgrid(:)-mean(X_Shifts(:)), pix_ygrid(:)-mean(Y_Shifts(:))];
    Image_rotation  = [(Image_shifts(:,1)-Centre_X),(Image_shifts(:,2)-Centre_Y)]*Rotation;
    Image_correct = [Image_rotation(:,1)+Centre_X, Image_rotation(:,2)+Centre_Y];
    
    F=scatteredInterpolant(Image_correct(:,1), Image_correct(:,2),Image(:),'natural');
    
    Image_re(:,:,I+1) =F(pix_xgrid,pix_ygrid);
    
end


