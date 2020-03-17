% main scipt of this DIC code
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1 sorting image number
Im_dir = dir(Image_folder); % Sets up image folder


% FIX: removes all non-image stuff from Im_dir structure, BP 2019.08.19
toRemove = zeros(1,length(Im_dir));
for ii = 1:length(Im_dir)
    if Im_dir(ii).isdir == 1
        toRemove(ii) = 1;
    elseif strncmp(Im_dir(ii).name(1:end),'DIC_Results',11) == 1
        %FIX: prevents this from reading previous results
        toRemove(ii) = 1;
    end
end
Im_dir(logical(toRemove)) = [];


Im_format = Im_dir(end).name(end-2:end); 
if strcmp(Im_format, 'tif')==1 
[Image_list,Image_num] = fDIC_image_sorting(Image_folder,'.tif', Image_Interval);
elseif strcmp(Im_format, 'jpg')==1 
[Image_list,Image_num] = fDIC_image_sorting(Image_folder,'.jpg', Image_Interval);
end
disp(['done 1 - sorting image numbers'])
disp(duration([0, 0, toc]))
disp('-------------------------')
%% Step2: convert the RBG to gray scale and correct the mosaic grids
% use the parallel analysis to accelerate computational speed
if imagecon ==1
    for i=1:Image_num
        Image_original = imread([Image_folder,Image_list{i}]);
        imwrite(Image_original,[Image_folder,Image_list{i}]);
    end
else
    fDIC_Image_convert(Image_folder,Image_num,Image_list)
end
disp('done 2 - convert RGB to gray')
disp(duration([0, 0, toc]))
disp('-------------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step3: set up ROIs
% set up the filter windows, boundary
filters                    = round([log2(ROI.size_pass_1)/2,log2(ROI.size_pass_1)/4,2*log2(ROI.size_pass_1),log2(ROI.size_pass_1)]);
Image_first                = im2double(imread([Image_folder,Image_list{1}]));
[Filters_setting,boundary] = gBounROI(Image_first,filters,ROI.size_pass_1); %this opens filter settings
Filters_setting            = round(Filters_setting); %this is needed within EBSD script
[FFTfilter,hfilter]        = fFilters(ROI.size_pass_1,Filters_setting); %converst to fourier space

% spread out the subregions in the defined ROI and set up the position of ROIs

[ROI.position_X_pass_1, ROI.position_Y_pass_1,ROI.num_x_pass_1,ROI.num_y_pass_1, ROI.coordinator_pass_1, ROI.num_pass_1] = fDIC_ROI_position(ROI.size_pass_1,ROI.overlap_pass_1,boundary);

disp('done 3 - set up ROIs')
disp(duration([0, 0, toc]))
disp('-------------------------')

% Step 4:  perform the XCF and determine shift in x  shift in y and peak height

[Shift_X_1,Shift_Y_1,CCmax_1] = fDIC_xcf(Image_folder,Image_list,Image_num,ROI,Filters_setting,XCF_mesh,hfilter,FFTfilter);

disp('done 4 - 1st round of XC')
disp(duration([0, 0, toc]))
disp('-------------------------')

%% Step 5: correct the image shifts and rotation based on the first pass shifts measurements.
if Image_correction ==1
    [Image_re] = fDIC_Image_correction(Image_folder, Image_list, Shift_X_1, Shift_Y_1, ROI.position_X_pass_1, ROI.position_Y_pass_1);
disp('done 5 - correct for shifts and rotations from 1st round of XC')
disp(duration([0, 0, toc]))
disp('-------------------------')  
%% Step 6: set up ROIs
% set up the filter windows, boundary
[FFTfilter,hfilter]        = fFilters(ROI.size_pass_2,Filters_setting);
[ROI.position_X_pass_2, ROI.position_Y_pass_2,ROI.num_x_pass_2,ROI.num_y_pass_2, ROI.coordinator_pass_2, ROI.num_pass_2] = fDIC_ROI_position(ROI.size_pass_2,ROI.overlap_pass_2,boundary);
disp('done 6 - set up ROIs for 2nd round of XC')
disp(duration([0, 0, toc]))
disp('-------------------------')
%% step 7: re-image cross correlation
    [Shift_X_2,Shift_Y_2,CCmax_2] = fDIC_xcf_2(Image_re,Image_num,Image_Interval, DIC_method,ROI,Filters_setting,XCF_mesh,hfilter,FFTfilter);
disp('done 7 - 2nd round of XC')
disp(duration([0, 0, toc]))
disp('-------------------------')
end
%% step 8: using user defined method to determine the strain values
Image_ref = imread([Image_folder Image_list{1}]);
Image_size = size(Image_ref);
if Image_correction ==1
    [Strain, Rotation,Strain_ef, F] = fDIC_StrainCalc(Shift_X_2, Shift_Y_2, ROI.position_X_pass_2, ROI.position_Y_pass_2, Image_size,strain_method,Image_num, DIC_method, Image_Interval);
else
    [Strain, Rotation, Strain_ef, F] = fDIC_StrainCalc(Shift_X_1, Shift_Y_1, ROI.position_X_pass_1, ROI.position_Y_pass_1, Image_size,strain_method,Image_num, DIC_method, Image_Interval);
end
disp('done 8 - calculate strains')
disp(duration([0, 0, toc]))
disp('-------------------------')
%%  Step 9 sort out data
[Data, Image_set, Setting] = fData_sort(Image_folder,Image_list,Image_num,Image_Interval,Image_ref,Image_re,...
                         Strain,Rotation,Strain_ef,F,Shift_X_2,Shift_Y_2,CCmax_2,...
                         FFTfilter,Filters_setting,ROI,DIC_method,Image_correction,XCF_mesh,imagecon, ...
                         hfilter,boundary,strain_method);
disp('done 9 - sort out the data')
disp(duration([0, 0, toc]))
disp('-------------------------')
%% step 10: display the final Shift in X and Y directions as arrow
%if DIC_method ==1
%    fDIC_figure(ROI.position_X, ROI.position_Y,Shift_X,Shift_Y,strain_11, strain_22, strain_12,CCmax,Image_ref,ROI,Image_folder, Image_num,DIC_method, Image_Interval,strain_range)
%else
%    %     fDIC_figure(ROI.position_X, ROI.position_Y,Shift_X,Shift_Y,strain_11, strain_22, strain_12,CCmax,Image_ref,ROI,Image_folder, Image_num,DIC_method, Image_Interval,strain_range)
%    fDIC_videos(Data, ROI.position_X_pass_2, ROI.position_Y_pass_2, Image_set, strain_range, rotation_range,ROI)
%end

% % sort out data again 
clear Image_folder Image_list Image_num Image_Interval Image_ref Image_re Image_size filter Image_first
clear Strain Rotation Strain_ef Shift_X_1 Shift_Y_1 Shift_X_2 Shift_Y_2 CCmax_1 CCmax_2
clear FFTfilter Filters_setting DIC_method Image_correction XCF_mesh imagecon Image_folder hfilter boundary strain_method
save([Image_set.Folder,'DIC_Results_ROI_', num2str(ROI.size_pass_2),'_',num2str(ROI.overlap_pass_2*100),'.mat'])

delete(gcp('nocreate'))
disp('done 10 - tidy up at the end')
disp(duration([0, 0, toc]))
disp('-------------------------')
disp('FINISHED')