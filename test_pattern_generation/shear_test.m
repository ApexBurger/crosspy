clear all


region_x=10; % dimensions of region to be patterned in um
region_y=10;

speckles=1; %speckles on or off
speckle_width=0.05; %in um
subset_dimension=0.4; %square subset dimension in um: There should be 3x3 speckles per subset. ie 9 per this no squared.
req_speckles=9; %no of required speckles per subset: 9? 5?

% speckle pattern generation parameters:
no=3; %no of speckles that affect shifts
%too high and av_distance gets too big; too low and not enough shift.
thresh_high=0.4; %of +/- of equilibrium separation that will shift
thresh_low=0.4;
mag=0.2; %of shift of micro-adjustments
its=8000; %no of iterations for micro-adjustments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%speckle shape for output image
shape=...
[0,0,1,0,0;
0,1,1,1,0;
1,1,1,1,1;
0,1,1,1,0;
0,0,1,0,0;];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stitch image parameters:
img_size=[1768, 2048]; %rows and columns so x, y reversed [884, 1024]
d_um2=req_speckles./(subset_dimension^2); %density of speckles: no. per um^2
pixel_width=speckle_width./3; %3 pixels per speckle determines recommended image size (in um per pixel)
hfw=img_size(2).*pixel_width; % j dimension is x - determines recommended hfw of stitch images.
%
% Initialise things in terms of pixels rather than space
region_i=round(region_y./pixel_width);
region_j=round(region_x./pixel_width);
d_pixels2=d_um2.*(pixel_width.^2); %density of speckles in no. per sq. pixel.= no. per um2 * pixelsperum^2
speckle_separation=round(sqrt(1/d_pixels2)); %speckle separation in pixels

%%%%%%%%%%%%%%%%%%%
% Set up region
%%% set up lattice of speckle centres with density set by d (and therefore speckle separation)
%%% each speckle then can move by a random number of pixels in x and y
%%% amount by which each pixel can move is 1/3 of ordered distance

% region is pattern without shape convolution
region=zeros(region_i,region_j);

for i=1:region_i;
    for j=1:region_j;
        
        % random offset between +/- 1/3 of speckle separation
        if mod(i,speckle_separation)==0 & mod(j,speckle_separation)==0;
            
           offset1=randi([-round(speckle_separation/3), round(speckle_separation/3)]);
           offset2=randi([-round(speckle_separation/3), round(speckle_separation/3)]);
           region(i+offset1-round(speckle_separation/3),j+offset2-round(speckle_separation/3))=1; 
        end
        j=j+1;
    end
    i=i+1;
end

clear offset1 offset2 i j a speckles_init 
original_region=region;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimise speckle pattern
%%%% do micro-adjustments on random speckles
f=waitbar(1./its, 'optimising...');

for index=1:its;
    waitbar(index./its);

    %[D,idx]=bwdist(image);
    [speckle_i,speckle_j]=find(region); %find non-zero elements and their indices
    speckle_locs=[speckle_i, speckle_j];
    clear speckle_i speckle_j
    clear distances distances_sorted distances_sorted_cropped speckles_closest sq_distances vectors net_vector loc speckle_pos_new;
    
% trial a random pixel and get closest speckles
    trial_i=randi([30,region_i-30]);
    trial_j=randi([30,region_j-30]);
    loc=[trial_i,trial_j];
    clear trial_i trial_j
    
    distances=speckle_locs-loc; %distances between trial location and closest speckles
    sq_distances=sqrt((distances(:,1)).^2+(distances(:,2)).^2); %square distances
    [distances_sorted(:,1),distances_sorted(:,2)]=sort(sq_distances); %outputs min. square distances and the index in the array
    distances_sorted_cropped=distances_sorted(1:no,:); %crops array to the speckle itself and no. closest
    
    for b=1:no; %generate closest speckles i,j indices
        speckles_closest(b,:)=speckle_locs(distances_sorted(b,2),:);
    end
    
    av_distance=mean(distances_sorted_cropped(:,1));

    %%%%% make adjustments to speckles
    if av_distance>(1+thresh_high)*speckle_separation;
        %move closest speckles to be slightly closer to trial location
        for b=1:no; %consider the requisite number of speckles specified at top
            vector=loc-speckles_closest(b,:); %vector from closest speckle to trial location
            region(speckles_closest(b,1),speckles_closest(b,2))=0; %remove considered speckle
            region(speckles_closest(b,1)+round(mag.*vector(1)),speckles_closest(b,2)+round(mag.*vector(2)))=1; %move considered speckle a bit closer
        end
    end
    
    if av_distance<(1-thresh_low)*speckle_separation;
        %move closest speckles to be slightly further from trial location
        for b=1:no; %consider the requisite number of speckles specified at top
            vector=loc-speckles_closest(b,:); %vector from closest speckle to trial location
            if speckles_closest(b,1)-round(mag.*vector(1))>0&speckles_closest(b,2)-round(mag.*vector(2))>0;
                region(speckles_closest(b,1),speckles_closest(b,2))=0; %remove considered speckle
                region(speckles_closest(b,1)-round(mag.*vector(1)),speckles_closest(b,2)-round(mag.*vector(2)))=1; %move considered speckle a bit further
            end
        end
    end
    
end
close(f)

clear a b speckles_closest index speckle_locs speckle_separation sq_distances loc speckle_pos_new c distances distances_sorted distances_sorted_cropped net_vector vectors f vector av_distance

%%%%%%%%%%%%%%%%

%%
% convolution of shape speckle at each centre (3x3 pixels optimal)
pattern_1=conv2(shape,original_region); % in case want to compare optimised / unoptimised
pattern_1=(pattern_1>0);

pattern=pattern_1;

shear_band_angle=[45,120,63,129,55];
shear_band_loc=[40,30,150,500,400];

sb_divide=zeros(size(pattern_1));

sb_x=5;
sb_y=3;

shear=zeros(size(sb_divide,2),size(sb_divide,1),length(shear_band_angle));
pattern_2=zeros(size(pattern));

for n=1:length(shear_band_angle)

for i=1:size(sb_divide,1)
    for j=1:size(sb_divide,2)
        
        %if i>shear_band_loc
            j_lim=(i-shear_band_loc(n))*tand(shear_band_angle(n));
            if j>j_lim
                try
                    pattern_2(i,j)=pattern(i+2,j+2);
                    shear(i,j,n)=1;
                end
            
            else
                pattern_2(i,j)=pattern(i,j);
            end
    end
end

pattern=pattern_2;

end

%% save patterns
pattern_1=pattern_1(1:size(pattern,1),1:size(pattern,2));

original=imgaussfilt(double(pattern_1));
deformed=imgaussfilt(double(pattern));

%%
imagesc(-original);
pbaspect([1,1,1])
axis off
colormap gray
print('original','-dpng','-r300');

imagesc(-deformed);
pbaspect([1,1,1])
axis off
colormap gray
print('deformed','-dpng','-r300');


