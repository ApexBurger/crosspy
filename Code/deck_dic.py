#%%
%load_ext autoreload
%autoreload 2

roi = dict(size_pass = 200, overlap_pass = 0.7, xcf_mesh=500)
filter_settings=[4,2,15,8] #fft filter settings: high pass, high pass width, low pass, low pass width

import os as o
o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
from Classes import *
from imprep_functions import *
from XCF_functions import *
from pathlib import Path
import matplotlib.pyplot as plt 

folder_path = Path(r"/Users/tom/Documents/GitHub/crosspy/data")
Images = Imset(folder_path,'tif')

# %% Instantiate and run the DIC

#build the dic class (but don't run it yet):
dic_1stpass=dic(Images,roi,filter_settings)

#run the dic on specified images within the stack, and get displacements:
dx_map, dy_map, ph_map = dic_1stpass.run(imnos=[0,1]) 

#can also do:
#dx_maps, dy_maps, ph_maps = dic_1stpass.run_sequential() #figures out imnos as consecutive images
#dx_maps, dy_maps, ph_maps = dic_1stpass.run_cumulative() #figures out imnos as 0 and sequential

# %%
