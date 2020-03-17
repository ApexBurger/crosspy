function [Image_list, Image_num] = fDIC_image_sorting(Image_folder,image_extension, Image_interval)
% this function is to put image into a correct list 
Image_list_temp      = dir(Image_folder);

num_images=length(Image_list_temp);
image_ok=false(num_images,1);
image_name_full=cell(num_images,1);


for n=1:num_images
    
%check dir
    if Image_list_temp(n).isdir == 0
        image_name=Image_list_temp(n).name;
        %check that file name has more chars than the extension
        if length(image_name) > length(image_extension)
            %check the file extension
            
            if strcmpi(image_name(end-length(image_extension)+1:end),image_extension)
                image_ok(n)=true;
                image_name_full{n}=image_name;
            end
        end
    end
end

Image_list_name=image_name_full(image_ok);

% correct the number sorting problem
Image_list_t = sort_nat(Image_list_name);
Image_num  = ceil(size(Image_list_name,1)/Image_interval);
for n=1: Image_num 
    Image_list{n,1} = Image_list_t{ceil((n-1)*Image_interval+1)};
end 