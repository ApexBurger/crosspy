# folder = 'asdf'
# Imset1 = Imset(folder)

# Imset.prepare(gray=True)



# displacements1,rotation1=get_displacements(Imset,start_im=0,end_im=1,FFT1=1234,FFT2=1234,ss=1234) #array

# strain1=get_strain(displacement1,rotation1,strain_method=0) #array
# plot_strains(strain_1)

#Settings:

roi = dict(
    size_pass = 200,
    overlap_pass = 70/100

)

#%%
%load_ext autoreload
%autoreload 2

import os as o
o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
from Classes import *
from pathlib import Path

folder_path = Path(r"/Users/tom/Documents/GitHub/crosspy/data")
Images = Imset(folder_path,'tif')

#example: this doesn't need to be in deck
ims = Images.imload([0,1])

#%% Generate filters

fftfil, hfil = Imset.gen_filters(roi[size_pass])
