
function [Shift_X,Shift_Y,CCmax] = fDIC_xcf(Image_folder,Image_list,Image_num,ROI,Filters_setting,XCF_mesh,hfilter,FFTfilter)
% cross correlation is performed and shift X, shift Y and peak height are
% determined.
Image_ref = im2double(imread([Image_folder Image_list{1}]));
Image_size = size(Image_ref);
for i=1:Image_num-1
    for j=1:ROI.num_pass_1
        ROI_ref_1 = Image_ref((ROI.coordinator_pass_1{j}(1,2)-ROI.size_pass_1/2):(ROI.coordinator_pass_1{j}(1,2)+ROI.size_pass_1/2-1),(ROI.coordinator_pass_1{j}(1,1)-ROI.size_pass_1/2):(ROI.coordinator_pass_1{j}(1,1)+ROI.size_pass_1/2-1));
        % zero mean and normalise standard deviation
        ROI_ref_1 = (ROI_ref_1- mean(ROI_ref_1(:)))./std(ROI_ref_1(:));
        ROI_ref_1 = ROI_ref_1.*hfilter; % han filtering
        ROI_ref_1 = fft2(ROI_ref_1); % 2D fast fourier transform
        data_fill =[1:(Filters_setting(3)+Filters_setting(4)),ROI.size_pass_1-(Filters_setting(3)+Filters_setting(4)-1):ROI.size_pass_1];
        data_fill = data_fill';
        ROI_ref = FFTfilter(data_fill,data_fill).*ROI_ref_1(data_fill,data_fill); % apply high and low frequence filter
        Image_test = im2double(imread([Image_folder Image_list{i+1}]));
        ROI_test_1 = Image_test((ROI.coordinator_pass_1{j}(1,2)-ROI.size_pass_1/2):(ROI.coordinator_pass_1{j}(1,2)+ROI.size_pass_1/2-1),(ROI.coordinator_pass_1{j}(1,1)-ROI.size_pass_1/2):(ROI.coordinator_pass_1{j}(1,1)+ROI.size_pass_1/2-1));
        % zero mean and normalise standard deviation
        ROI_test_1 = (ROI_test_1- mean(ROI_test_1(:)))./std(ROI_test_1(:));
        ROI_test_1 = ROI_test_1.*hfilter; % han filtering
        ROI_test_1 = fft2(ROI_test_1); % 2D fast fourier transform
        ROI_test = FFTfilter(data_fill,data_fill).*ROI_test_1(data_fill,data_fill);
        [col_shift(j),row_shift(j),ccmax(j)] = fReg(ROI_test,ROI_ref, ROI.size_pass_1,XCF_mesh, data_fill);
    end
    Shift_X_temp(:,i) = col_shift(:);
    Shift_Y_temp(:,i) = row_shift(:);
    CCmax(:,i)   = ccmax(:);
    Shift_X(:,:,i) = reshape(Shift_X_temp(:,i),size(ROI.position_X_pass_1));
    Shift_Y(:,:,i) = reshape(Shift_Y_temp(:,i),size(ROI.position_X_pass_1));
end
