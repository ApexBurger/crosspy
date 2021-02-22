#%%

import matplotlib.pyplot as plt
import numpy as np 
from crosspy import plan_ffts, gen_filters, freg, fxcorr
from crosspy import DIC,Imset

# %% create some data

f = np.random.rand(32,32)
g = np.random.rand(32,32)

# Fine pass settings
fine = 32 # in pixels
overlap = 90 # in percentage

Images = np.dstack((f,g))
# fft filter settings: high pass, high pass width, low pass, low pass width
filter_settings = [4,2,15,8]
roi_1stpass = dict(size_pass = fine, overlap_percentage = overlap, xcf_mesh=fine*2)
d = DIC(Images,roi_1stpass,filter_settings)

# %%
prepffts = plan_ffts(d)    
fxcorr(f,g,d,prepared_ffts=prepffts)
# %%

def denom_precalc(f,g,d):
    """
    """
    pass