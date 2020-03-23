%%% Design the matlab DIC script for 2-D strain measurements%%%
%%% created by Jun Jiang ----16/01/2014

%%load the saved tiff images into matlab
isOpen = isempty(gcp('nocreate'));
if isOpen~=1
    parpool 
else 
    delete(gcp('nocreate'))
     parpool 
end
clear all
close all
addpath code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%User input%%%%%%%%%%%%%%%%%%
Image_folder    ='D:\DIC\crosspy\data\Tom\'; % '2' is using the consetutive image as the refernce image
DIC_method      = 2; % '1' is used as the first image as the reference image (suitable for small deformation)
% '2' is used as the previous image as the reference image (prefered for large deformation)
Image_Interval  = 1; %define the interval between the reference and test pattern for Method 2
ROI.size_pass_1    = 200; % pass  to fix rigid body rotation
ROI.overlap_pass_1 = 70; % in percent
ROI.size_pass_2    = 60; % pass which determines final strain resolution
ROI.overlap_pass_2 = 80; % in percent
ROI.overlap_pass_1 = ROI.overlap_pass_1/100;
ROI.overlap_pass_2 = ROI.overlap_pass_2/100;
XCF_mesh        = 500; % XCF mesh size default value  is 250
strain_method   = 4;%1---9nodes FEA; 2---8 nodes FEA; 3---4nodes FAE; 4---least squre poly fitting quaratic in the bulk and linear at the edges and corners
% method 4 is recommanded as it used all ROIs and applied but first and
% second orders of polynomial at the edges and bulk of image
imagecon         = 0; % recorded from different camera 1==without conversion
Image_correction = 1; % 1 for switch on correct image shifts and rotatio
strain_range   = [-0.5 0.5];% display strain range the noise floor is ~0.01
rotation_range = [-0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run  main functions 
main






