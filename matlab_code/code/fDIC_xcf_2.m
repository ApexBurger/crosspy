
function [Shift_X_2,Shift_Y_2,CCmax_2] = fDIC_xcf_2(Image_re, Image_num, Image_Interval, method,ROI,Filters_setting,XCF_mesh,hfilter,FFTfilter)
% cross correlation is performed and shift X, shift Y and peak height are
% determined.
if method ==1
    Image_ref     = Image_re(:,:,1);
    Image_size = size(Image_ref);
    for i=1:ROI.num_pass_2
        ROI_ref_1 = Image_ref((ROI.coordinator_pass_2{i}(1,2)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{i}(1,2)+ROI.size_pass_2/2-1),(ROI.coordinator_pass_2{i}(1,1)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{i}(1,1)+ROI.size_pass_2/2-1));
        % zero mean and normalise standard deviation
        ROI_ref_1 = (ROI_ref_1- mean(ROI_ref_1(:)))./std(ROI_ref_1(:));
        ROI_ref_1 = ROI_ref_1.*hfilter; % han filtering
        ROI_ref_1 = fft2(ROI_ref_1); % 2D fast fourier transform
        data_fill =[1:(Filters_setting(3)+Filters_setting(4)),ROI.size_pass_2-(Filters_setting(3)+Filters_setting(4)-1):ROI.size_pass_2];
        ROI_ref(:,:,i) = FFTfilter(data_fill,data_fill).*ROI_ref_1(data_fill,data_fill); % apply high and low frequence filter
    end
    
    for i=1:Image_num-1
         parfor j=1:ROI.num_pass_2
            Image_test = Image_re(:,:,i);
            ROI_test_1 = Image_test((ROI.coordinator_pass_2{j}(1,2)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{j}(1,2)+ROI.size_pass_2/2-1),(ROI.coordinator_pass_2{j}(1,1)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{j}(1,1)+ROI.size_pass_2/2-1));
            % zero mean and normalise standard deviation
            ROI_test_1 = (ROI_test_1- mean(ROI_test_1(:)))./std(ROI_test_1(:));
            ROI_test_1 = ROI_test_1.*hfilter; % han filtering
            ROI_test_1 = fft2(ROI_test_1); % 2D fast fourier transform
            ROI_test(:,:,j) = FFTfilter(data_fill,data_fill).*ROI_test_1(data_fill,data_fill);
            [col_shift(j),row_shift(j),ccmax(j)] = fReg(ROI_test(:,:,j),ROI_ref(:,:,j),ROI.size_pass_2,XCF_mesh, data_fill);
        end
        Shift_X_temp(:,i) = col_shift(:);
        Shift_Y_temp(:,i) = row_shift(:);
        CCmax_2(:,i)   = ccmax(:);
        Shift_X_2(:,:,i) = reshape(Shift_X_temp(:,i),size(ROI.position_X_pass_2));
        Shift_Y_2(:,:,i) = reshape(Shift_Y_temp(:,i),size(ROI.position_X_pass_2));
    end
    
elseif method ==2
    
    for i =1:Image_num-1
        %%%% prepare two consecutive images Han filter, 2D-FFT
        parfor j=1:ROI.num_pass_2
            % For the first image
            Image_ref     = Image_re(:,:,i);
            ROI_ref_1     = Image_ref((ROI.coordinator_pass_2{j}(1,2)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{j}(1,2)+ROI.size_pass_2/2-1),(ROI.coordinator_pass_2{j}(1,1)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{j}(1,1)+ROI.size_pass_2/2-1));
            % zero mean and normalise standard deviation
            ROI_ref_1     = (ROI_ref_1- mean(ROI_ref_1(:)))./std(ROI_ref_1(:));
            ROI_ref_1     = ROI_ref_1.*hfilter; % han filtering
            ROI_ref_1     = fft2(ROI_ref_1); % 2D fast fourier transform
            data_fill     =[1:(Filters_setting(3)+Filters_setting(4)),ROI.size_pass_2-(Filters_setting(3)+Filters_setting(4)-1):ROI.size_pass_2];
            ROI_ref       = FFTfilter(data_fill,data_fill).*ROI_ref_1(data_fill,data_fill); % apply high and low frequence filter
            % for the following image
            Image_test    = Image_re(:,:,i+1);
            ROI_test_1 = Image_test((ROI.coordinator_pass_2{j}(1,2)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{j}(1,2)+ROI.size_pass_2/2-1),(ROI.coordinator_pass_2{j}(1,1)-ROI.size_pass_2/2):(ROI.coordinator_pass_2{j}(1,1)+ROI.size_pass_2/2-1));
            % zero mean and normalise standard deviation
            ROI_test_1 = (ROI_test_1- mean(ROI_test_1(:)))./std(ROI_test_1(:));
            ROI_test_1 = ROI_test_1.*hfilter; % han filtering
            ROI_test_1 = fft2(ROI_test_1); % 2D fast fourier transform
            ROI_test   = FFTfilter(data_fill,data_fill).*ROI_test_1(data_fill,data_fill);
            [col_shift(j),row_shift(j),ccmax(j)] = fReg(ROI_test,ROI_ref,ROI.size_pass_2,XCF_mesh, data_fill);
        end
        Shift_X_temp(:,i) = col_shift(:);
        Shift_Y_temp(:,i) = row_shift(:);
        CCmax_2(:,i)   = ccmax(:);
        Shift_X_2(:,:,i) = reshape(Shift_X_temp(:,i),size(ROI.position_X_pass_2));
        Shift_Y_2(:,:,i) = reshape(Shift_Y_temp(:,i),size(ROI.position_X_pass_2));
        if i>1
            Shift_X_2(:,:,i) = Shift_X_2(:,:,i-1)+Shift_X_2(:,:,i);
            Shift_Y_2(:,:,i) = Shift_Y_2(:,:,i-1)+Shift_Y_2(:,:,i);
        end
    end
end