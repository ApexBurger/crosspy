#Notes / plans:

# displacements1,rotation1=get_displacements(Imset,start_im=0,end_im=1,FFT1=1234,FFT2=1234,ss=1234) #array

# strain1=get_strain(displacement1,rotation1,strain_method=0) #array
# plot_strains(strain_1)

#%%
%load_ext autoreload
%autoreload 2

roi = dict(size_pass = 200, overlap_pass = 0.7, xcf_mesh=250)
filters_settings=[4,2,15,8] #fft filter settings: high pass, high pass width, low pass, low pass width

import os as o
o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
from Classes import *
from DataPrep_Functions import *
from fDIC_functions import *
from pathlib import Path
import matplotlib.pyplot as plt 

folder_path = Path(r"/Users/tom/Documents/GitHub/crosspy/data")
Images = Imset(folder_path,'tif')

#examples: these dont need to be in deck
ims = Images.imload([0,1])
ss_locations=gen_ROIs(ims.shape[0:2],roi)

#grab a couple of subsets
subset1=get_subset(ims,roi,ss_locations,0,0)
subset2=get_subset(ims,roi,ss_locations,0,1)

#calculate an xcf peak height
dx,dy,ph=fxcorr(subset1,subset2,roi,filters_settings)

# %%
