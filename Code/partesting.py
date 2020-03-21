#%%
%load_ext autoreload
%autoreload 2
#%%
import os as o
o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
from Classes import *
from imprep_functions import *
from XCF_functions import *
from pathlib import Path

from numba import jit

folder_path = Path(r"/Users/tom/Documents/GitHub/crosspy/data/Siyang")
Images = Imset(folder_path,'tif')

# %% Instantiate and run the DIC

#fft filter settings: high pass, high pass width, low pass, low pass width
filter_settings=[4,2,15,8]
roi_1stpass = [200,50,400] #subset, overlap percentage, xcf mesh

#%%build the dic class (but don't run it yet):
d=dic(Images,roi_1stpass,filter_settings)

#%% test some functions
gen_ROIs((1700,1700),roi_1stpass)
ref=get_subset(d.ims,200,d.ss_locations,0,1)
test=get_subset(d.ims,200,d.ss_locations,0,1)

dxs[subset_n],dys[subset_n],phs[subset_n]=fxcorr(ref,test,roi,fftfil,hfil,filters_settings)

#%%

def run(d,imnos=[0,1]):
        
    phs=np.zeros(d.n_subsets)
    dxs=np.zeros(d.n_subsets)
    dys=np.zeros(d.n_subsets)

    for subset_n in range(0,d.n_subsets):
        #grab the reference and test subsets, and get subpixel registration
        ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
        test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
        dxs[subset_n],dys[subset_n],phs[subset_n]=fxcorr(ref,test,d)

        #translate best_dxs etc back onto image grid
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

    return dx_map,dy_map,ph_map

#%%
import time
start=time.time()
run(d,imnos=[0,1])
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))

import time
start=time.time()
run(d,imnos=[0,1])
end=time.time()
print("Elapsed (with compilation) = %s" % (end - start))

# %%
