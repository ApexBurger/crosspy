#%%
%load_ext autoreload
%autoreload 2

roi = dict(size_pass = 200, overlap_pass = 0.7, xcf_mesh=500)
filter_settings=[4,2,15,8] #fft filter settings: high pass, high pass width, low pass, low pass width

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

# %% Instantiate and run the DIC

#build the dic class (but don't run it yet):
dic_1stpass=DIC(Images,roi,filter_settings)

#run the dic on specified images within the stack, and get displacements:
dx_map, dy_map, ph_map = dic_1stpass.run(imnos=[0,1]) 

#can also do:
#dx_maps, dy_maps, ph_maps = dic_1stpass.run_sequential() #figures out imnos as consecutive images
#dx_maps, dy_maps, ph_maps = dic_1stpass.run_cumulative() #figures out imnos as 0 and sequential

#%% Generate filters - ALEX

fftfil, hfil = Imset.gen_filters(roi['size_pass'])

#%% Run cross correlation - TOM

shift_x, shift_y, peakcc = xcf()

#%% Image correction - ALEX

image_cor = im_correct(Images, shift_x, shift_y, ss_locations)

#%% Re-run cross correlation with smaller subset size - TOM

#%% calculate strain - ALEX
