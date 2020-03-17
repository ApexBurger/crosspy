function [Data, Image_set, Setting] = fData_sort(Image_folder,Image_list,Image_num,Image_Interval,Image_ref,Image_re,...
    Strain,Rotation,Strain_ef,F,Shift_X_2,Shift_Y_2,CCMax_2,...
    FFTfilter,Filters_setting,ROI,DIC_method,Image_correction,XCF_mesh,imagecon, ...
    hfilter,boundary,strain_method)
% function to tidy up variables
% created by jun jiang 01/2015

Data.Shift_X     = Shift_X_2;
Data.Shift_Y     = Shift_Y_2;
Data.CCmax       = CCMax_2;
Data.Strain      = Strain;
Data.Rotation    = Rotation;
Data.Strain_ef   = Strain_ef;
Data.F           = F;

Image_set.Folder = Image_folder;
Image_set.List   = Image_list;
Image_set.Num    = Image_num;
Image_set.Interval = Image_Interval;

Image_set.Ref     = Image_ref;
Image_set.Re      = Image_re;

Setting.FFTfilter = FFTfilter;
Setting.Filters_setting = Filters_setting;
Setting.ROI    = ROI;
Setting.method = DIC_method;
Setting.Correct_image = Image_correction;
Setting.mesh_size = XCF_mesh;
Setting.imagecon = imagecon;
Setting.hfilter = hfilter;
Setting.boundary = boundary;
Setting.strain_method = strain_method;
