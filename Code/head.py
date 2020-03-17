# Head script for DIC

# Preamble
import dpyc as dic

#%% Step 1 - Import images and convert to grayscale

im1 = dic.im(file1)
im2 = dic.im(file2)

im1 = im1.convg()
im2 = im2.convg()
# The class should contain the image number

# sort by DIm.imnum

print('Finished sorting images and converting to grayscale')
#%% Step 2 - Set up ROIs
# filter windows and boundaries

filters = (roi[size_pass_1])
filters_setting, boundary = dic.gboundroi(filters)

fftfilter, hfilter = dic.filters(roi1, filters_setting)

roi1 = dic.roi1()
#%% Step 3 - Perform xcf and determine shift in x, shift in y
# and peak height

xcf1 = dic.xcf(im1,im2)

print('Finished 1st round of XC')

#%% Step 4 - Correct image shifts and rotation
if image_correction == 1:
    im_re = dic.im_correction()

print('Finished image corrections')
#%% Step 5 - Set up filter windows, boundary

roi2 = dic.roi2()

print('Finished setting up ROIs for 2nd round of XC')

#%% Step 6 - 2nd round of xc

xcf2 = dic.xcf2(im_re)

print('Finished second round of XC')

#%% Step 7 - Determine strain 

strains = dic.straincalc(xcf2)
print('Finished calculating strains')

#%% Step 8 - Sort out data?

