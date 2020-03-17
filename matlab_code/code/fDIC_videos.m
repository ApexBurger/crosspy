function fDIC_videos(Data, position_X, position_Y, Image_set, strain_range, rotation_range,ROI)
% convert the strain/ displacement maps into a gif video
% created by Jun Jiang ----04/2014
% partition data ointo individual terms
Image_folder   = Image_set.Folder;
Image_list     = Image_set.List;
Image_num      = Image_set.Num;
Image_Interval = Image_set.Interval;
Image_ref      = Image_set.Ref;
Image_re       = Image_set.Re;

for n=1: size(Data.Strain,2)   
    S11(:,:,n) = Data.Strain{n}(:,:,1,1);
    S12(:,:,n) = Data.Strain{n}(:,:,1,2);
    S13(:,:,n) = Data.Strain{n}(:,:,1,3);
    
    S21 (:,:,n)= Data.Strain{n}(:,:,2,1);
    S22(:,:,n) = Data.Strain{n}(:,:,2,2);
    S23(:,:,n) = Data.Strain{n}(:,:,2,3);
    
    S31(:,:,n) = Data.Strain{n}(:,:,3,1);
    S32(:,:,n) = Data.Strain{n}(:,:,3,2);
    S33(:,:,n) = Data.Strain{n}(:,:,3,3);
    
    W11(:,:,n) = real(Data.Rotation{n}(:,:,1,1));
    W12(:,:,n) = real(Data.Rotation{n}(:,:,1,2));
    W13(:,:,n) = real(Data.Rotation{n}(:,:,1,3));
    
    W21(:,:,n) = real(Data.Rotation{n}(:,:,2,1));
    W22(:,:,n) = real(Data.Rotation{n}(:,:,2,2));
    W23(:,:,n) = real(Data.Rotation{n}(:,:,2,3));
    
    W31(:,:,n) = real(Data.Rotation{n}(:,:,3,1));
    W32(:,:,n) = real(Data.Rotation{n}(:,:,3,2));
    W33(:,:,n) = real(Data.Rotation{n}(:,:,3,3));
end

Shift_X = Data.Shift_X;
Shift_Y = Data.Shift_Y;

S_ef        = Data.Strain_ef;
% 1. Plot [3x3] Strain tensor
for k=1:size(Data.Strain,2)
    figure ('color',[1 1 1])
    subplot(3,3,1) % plot S11
    imagesc(S11(:,:,k));
    axis image, axis off
    caxis(strain_range), colormap('jet')
    title('\epsilon^p_1_1', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,2) % plot S12
    imagesc(S12(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_1_2', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,3)  % plot S21
    imagesc(S13(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_1_3', 'fontsize', 18, 'fontweight','bold')
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)
    
    subplot(3,3,4)  % plot S22
    imagesc(S21(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_2_1', 'fontsize', 18, 'fontweight','bold')
   
    subplot(3,3,5)  % plot S22
    imagesc(S22(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_2_2', 'fontsize', 18, 'fontweight','bold')
     
    subplot(3,3,6)  % plot S22
    imagesc(S23(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_2_3', 'fontsize', 18, 'fontweight','bold')
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)
   
    subplot(3,3,7)  % plot S22
    imagesc(S31(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_3_1', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,8)  % plot S22
    imagesc(S32(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_3_2', 'fontsize', 18, 'fontweight','bold')       
       
    subplot(3,3,9)  % plot S22
    imagesc(S33(:,:,k));
    axis image, axis off
    caxis(strain_range),colormap('jet')
    title('\epsilon^p_3_3', 'fontsize', 18, 'fontweight','bold')
     Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)
    
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)   
    set(gcf,'position',[300 150 1200 800])
%     suptitle('Plastic strain')
    
    h = gcf;
    f_S = getframe(h);
     [~,map_temp]    = rgb2ind(f_S.cdata,256,'nodither');
     map_S(:,:,k)    = map_temp(1:98,:);
     im_S(:,:,1,k) = rgb2ind(f_S.cdata,map_S(:,:,k),'nodither');
%     im_S_test(:,:,k) = rgb2ind(f_S.cdata,map_S,'nodither');
end
% Write out the image file.
imwrite(im_S,map_S,([Image_folder,'Strain_ROI_',num2str(ROI.size_pass_2),'_',num2str(ROI.overlap_pass_2*100),'%.gif']),'DelayTime',0.1,'LoopCount',inf);

% plot [2x2] Rotation tensor
for k=1:size(Data.Strain,2)
    figure('color',[1 1 1])
    subplot(3,3,1) % plot W11
    imagesc(W11(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_1_1', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,2) % plot W12
    imagesc(W12(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_1_2', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,3)  % plot W21
    imagesc(W13(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_1_3', 'fontsize', 18, 'fontweight','bold')
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)
    
    subplot(3,3,4)  % plot W22
    imagesc(W21(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_2_1', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,5)  % plot W22
    imagesc(W22(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_2_2', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,6)  % plot W22
    imagesc(W23(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_2_3', 'fontsize', 18, 'fontweight','bold')
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)
    
    subplot(3,3,7)  % plot W22
    imagesc(W31(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_3_1', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,8)  % plot W22
    imagesc(W32(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_3_2', 'fontsize', 18, 'fontweight','bold')
    
    subplot(3,3,9)  % plot W22
    imagesc(W33(:,:,k));
    axis image, axis off
    caxis(rotation_range),colormap('jet')
    title('\omega_3_3', 'fontsize', 18, 'fontweight','bold')
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set (gca,'Position',Orig_Posi)
    
    set(gcf,'position',[300 150 1200 800])

%     suptitle('Plastic Rotation')
    
    h= gcf;
    f_W = getframe(h);
    [~,map_temp] = rgb2ind(f_W.cdata,256,'nodither');
    map_W(:,:,k)  = map_temp(1:98,:);
    im_W(:,:,1,k) = rgb2ind(f_W.cdata,map_W(:,:,k),'nodither');
end
% Write out the image file.
imwrite(im_W,map_W,([Image_folder,'Rotation_ROI_', num2str(ROI.size_pass_2),'_',num2str(ROI.overlap_pass_2*100),'%.gif']),'DelayTime',0.1,'LoopCount',inf);

% plot effective strain
for k=1: size(Data.Strain,2)
    figure('color',[1 1 1]) % plot W22
    imagesc(S_ef(:,:,k));
    axis image, axis off
    caxis([0 strain_range(2)]),colormap('jet')
    title('Effective Strain (\epsilon_e_f_f_e_c_t_i_v_e)', 'fontsize', 18, 'fontweight','bold')
    
    Orig_Posi = get(gca,'Position'); % sort out the image resize problem
    colorbar('fontsize',14, 'fontweight','bold','location','EastOutside')
    set(gca,'Position',Orig_Posi)    
    set(gcf,'position',[300 150 1200 800])
    
    h= gcf;
    f_S_ef = getframe(h);
    [~, map_temp] = rgb2ind(f_S_ef.cdata,256,'nodither');
    map_S_ef(:,:,k) = map_temp(1:413,:); 
    im_S_ef(:,:,1,k) = rgb2ind(f_S_ef.cdata,map_S_ef(:,:,k),'nodither');
end
% Write out the image file.
imwrite(im_S_ef,map_S_ef,([Image_folder,'Strain_effective_ROI_', num2str(ROI.size_pass_2),'_',num2str(ROI.overlap_pass_2*100),'%.gif']),'DelayTime',0.1,'LoopCount',inf);

%% displacement video
for k=1: size(Data.Strain,2)
    figure('color',[1 1 1])
    imagesc(Image_re(:,:,k+1)), colormap('gray'), axis off
    caxis([0 1])
    hold on
    Shift_X_temp = Shift_X(:,:,k);
    Shift_Y_temp = Shift_Y(:,:,k);
    quiver(position_X(:), position_Y(:),Shift_X_temp(:), Shift_Y_temp(:),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',10,'linewidth',3), axis image
    hold off
    title('Displacement', 'fontsize', 18, 'fontweight','bold')
    h= gcf;
    f_D = getframe(h);
    [~,map_D(:,:,k)] = rgb2ind(f_D.cdata,256,'nodither');
    im_D(:,:,1,k) = rgb2ind(f_D.cdata,map_D(:,:,k),'nodither');
end
% Write out the image file.
imwrite(im_D,map_D,([Image_folder,'Displacement_ROI_',num2str(ROI.size_pass_2),'_',num2str(ROI.overlap_pass_2*100),'%.gif']),'DelayTime',1,'LoopCount',inf);

%% Test image video
% for k=1: size(Data.Strain,2)
%     figure('color',[1 1 1])
%     imagesc(Image_re(:,:,k)), colormap('gray'), axis image, axis equal, axis off
%     title('Test Images', 'fontsize', 16, 'fontweight','bold')
%     caxis([0 1])
%     h= gcf;
%     f_test = getframe(h);
%     [~,map_test(:,:,k)] = rgb2ind(f_test.cdata,256,'nodither');
%     im_test(:,:,1,k) = rgb2ind(f_test.cdata,map_test(:,:,k),'nodither');
% end
% Write out the image file.
% imwrite(im_test,map_test,([Image_folder,'Image_ROI_', num2str(ROI.size_pass_2),'_',num2str(ROI.overlap_pass_2*100),'%.gif']),'DelayTime',1,'LoopCount',inf);

close all

