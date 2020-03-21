%%DIC process
% plots the figures from HR_DIC
clear all
%% settings
% automatically calculate colourscale min and max values from data
% 28/08/19 doesn't work yet! - Might never work... 
c_lim_auto = 0;

% get the colourmaps 
addpath('Colormaps')
%% get the files
% upload the .mat file 
raw_data_dir = 'D:\DIC\crosspy\data\';
[sample_file,sample_path] = uigetfile([raw_data_dir '\*.mat'],'Where the raw DIC results are');

%% 
full_file = [sample_path sample_file];


% get sample name
sample_name = sample_file(1:end-4);
disp(['Processing sample ' sample_name])

% get the DIC results
Res_data = load(full_file);
Res_data.sample = sample_name;
disp(['Using results file ' sample_file])

% make new directory
this_results_dir = [sample_path sample_name];
if ~exist(this_results_dir, 'dir')
    mkdir(this_results_dir)
end


    
%% plot some maps
% strain
comp_1 = 1; %i.e. x = 1, y = 2, z = 3 (remember z is bollocks as it's out of plane)
comp_2 = 1;

% which components
if c_lim_auto == 1
    c_max = max(max(Res_data.Data.Strain{1}(:,:,comp_1,comp_2)));
    c_min = min(min(Res_data.Data.Strain{1}(:,:,comp_1,comp_2)));
else
    c_max = 0.01;
    c_min = 0;
end

figure
% plot
%imagesc(Res_data.Data.Strain{1}(:,:,comp_1,comp_2))
imagesc(Res_data.Data.Strain_ef)
% settings
caxis([c_min c_max])
fig_name = ['e_' num2str(comp_1) '_,_' num2str(comp_2)];
title(fig_name)
fig_name = erase(fig_name,'_');
fig_name = erase(fig_name,',');
axis off
axis image
colormap jet
colorbar
disp(['Plotting ' fig_name])

% print the figures
% in tif format for now
print([this_results_dir '\' fig_name],'-dtiff')

% in fig format for later editing
savefig([this_results_dir '\' fig_name])


disp('All done!')





%% playing around with denoising these images
% figure
% imagesc(wiener2(Res_data.Data.Strain{1}(:,:,comp_1,comp_2),[5,5]))
% caxis([c_min c_max])
% fig_name = ['e_' num2str(comp_1) '_,_' num2str(comp_2) ' denoised'];
% title(fig_name)
% fig_name = erase(fig_name,'_');
% fig_name = erase(fig_name,',');
% axis off
% axis image
% colormap magma
% colorbar