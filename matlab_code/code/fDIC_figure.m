function fDIC_figure(position_X, position_Y,Shift_X,Shift_Y, strain_11, strain_22, strain_12,CCmax, Image_ref, ROI, Image_folder,Image_num, DIC_method, Image_interval, strain_range)
% this function is used for plotting figure such as shift quiver plot and
% strain mapps
if DIC_method ==1
    Image_num_DIC     = Image_num;
else
    Image_num_DIC     = ceil(Image_num/Image_interval)-1;
end
    
for i=1:Image_num_DIC
    % plot displacement map
    figure('color',[1 1 1])
    imagesc(Image_ref), colormap('gray'), axis off
    hold on
    Shift_X_temp = Shift_X(:,:,i);
    Shift_Y_temp = Shift_Y(:,:,i);
    CCmax_temp   = CCmax(:,i);
    quiver(position_X(:), position_Y(:),Shift_X_temp(:), Shift_Y_temp(:),'marker','o','markerfacecolor','k','markeredgecolor','k','markersize',18,'linewidth',3), axis image
    hold off
    F1= gcf;
    imwrite(F1,[Image_folder,'displacement','_', num2str(i),'.tif'])
    % plot correlation peak height map
    figure; scatter(position_X(:),position_Y(:),20,CCmax_temp(:),'filled'); axis ij;
    F2=gcf;
    imwrite(F2,[Image_folder,'PH','_', num2str(i),'.tif'])
    axis off
    % plot strain maps
    figure('color',[1 1 1])
    subplot(2,2,1)
    imagesc(strain_11(:,:,i)), axis image,  axis off
    caxis(strain_range);
    title('\epsilon _1_1','fontsize',12,'fontweight','bold')
    subplot(2,2,4)
    imagesc(strain_22(:,:,i)), axis image, axis equal, axis off
    caxis(strain_range);
    title('\epsilon _2_2','fontsize',12,'fontweight','bold')
    subplot(2,2,3)
    imagesc(strain_12(:,:,i)), axis image, axis equal, axis off
    caxis(strain_range);
    title('\epsilon _1_2','fontsize',12,'fontweight','bold')
    subplot(2,2,2)
    text(0.1,0.35,['ROI size is', ROI.size],'fontsize',10)
    text(0.1,0.55,['ROI overlapping is', ROI.overlap*100,'%'], 'fontsize',10)
    
    axis off
    colorbar('fontsize',12,'fontweight','bold')
    caxis(strain_range);
    F3=gcf;
    imwrite(F3,[Image_folder,'2D_srain','_', num2str(i),'.tif'])
end

