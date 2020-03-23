#%%
%load_ext autoreload
%autoreload 2
#%%
import os as o
#o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
o.chdir('D:/DIC/crosspy/Code')
from Classes import *
from imprep_functions import *
from XCF_functions import *
from pathlib import Path
import matplotlib.pyplot as plt 
from ImageCorrection_functions import *

folder_path = Path(r"D:/DIC/crosspy/data/Tom")
Images = Imset(folder_path,'tif')

fig = plt.figure()
plt.imshow(Images.imload([1]))

# %% Instantiate and run the DIC

#fft filter settings: high pass, high pass width, low pass, low pass width
filter_settings=[4,2,15,8]
roi_1stpass = dict(size_pass = 200, overlap_pass = 0.7, xcf_mesh=500)

#build the dic class (but don't run it yet):
dic_1stpass=DIC(Images,roi_1stpass,filter_settings)

#run the dic on specified images within the stack, and get displacements:
dic_1stpass.run_sequential() #figures out imnos as consecutive images

#can also do:
#dx_maps, dy_maps, ph_maps = dic_1stpass.run_cumulative() #figures out imnos as 0 and sequential

#%% Show the displacements
dic_1stpass.plot_results()

#%% Image correction
Images_cor =im_correct(Images, dic_1stpass.dx_maps, dic_1stpass.dy_maps, dic_1stpass.x_pos, dic_1stpass.y_pos)

#%% Second pass
roi_2ndpass = dict(size_pass = 200, overlap_pass = 0.7, xcf_mesh=500)
dic_2ndpass=dic(Images_cor,roi_2ndpass,filter_settings)

dx_maps, dy_maps, ph_maps = dic_2ndpass.run_sequential() #figures out imnos as consecutive images

#%% Image correction - ALEX

image_cor = im_correct(Images, shift_x, shift_y, ss_locations)

#%% calculate strain - ALEX
