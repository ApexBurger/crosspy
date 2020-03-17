function fDIC_Image_convert(Image_folder, Image_num,Image_list)
% remove the mosaic format and save image as tiff in the original folder,
% no compression is used
for i=1:Image_num 
    Image_original = imread([Image_folder Image_list{i}]);
%    Image_demosaic = demosaic(Image_original,'gbr');
if size(size(Image_original),2)>2
    Image_gray     = rgb2gray(Image_original); 
    imwrite(Image_gray,[Image_folder,Image_list{i}]);
else 
    imwrite(Image_original,[Image_folder,Image_list{i}]);
end 
end