
#%%
%load_ext autoreload
%autoreload 2

#%%
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy import DIC, Imset

import numpy as np 
import time

import os as o 
o.chdir(r'/Users/tom/Documents/GitHub/crosspy/testing')

import partesting_functions as p 

import pyfftw
import multiprocessing

if __name__=='__main__':

    #folder_path = Path(r"C:\Users\tpm416\Documents\GitHub\crosspy\data\Siyang")
    folder_path=Path(r'/Users/tom/Documents/GitHub/crosspy/data/Siyang/')
    Images = Imset(folder_path,'tif')

    # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = 200, overlap_percentage = 70, xcf_mesh=250)

    # build the dic class (but don't run it yet):
    dic_1stpass=DIC(Images,roi_1stpass,filter_settings)
    
    def run_DIC(d,imnos=[0,1],par=False,cores=None):
        
        #preallocate for this DIC pair
        phs=np.zeros(d.n_subsets)
        dxs=np.zeros(d.n_subsets)
        dys=np.zeros(d.n_subsets)

        prepared_ffts=p.plan_ffts(d)
        original_ffts=[np.fft.fft2,np.fft.ifft2,np.fft.fft2,np.fft.ifft2]

        pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
        pyfftw.interfaces.cache.enable()

        np_ffts=[pyfftw.interfaces.numpy_fft.fft2,pyfftw.interfaces.numpy_fft.ifft2,pyfftw.interfaces.numpy_fft.fft2,pyfftw.interfaces.numpy_fft.ifft2]
        scipy_ffts=[pyfftw.interfaces.scipy_fftpack.fft2,pyfftw.interfaces.scipy_fftpack.ifft2,pyfftw.interfaces.scipy_fftpack.fft2,pyfftw.interfaces.scipy_fftpack.ifft2]

        for subset_n in range(0,d.n_subsets):
            dxs[subset_n],dys[subset_n],phs[subset_n]=p.subset_compare_1(d,imnos,subset_n,np_ffts)

        #translate best_dxs etc back onto image grid
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

        return dx_map,dy_map,ph_map

    t0=time.time()
    print('starting')
    output=run_DIC(dic_1stpass)
    print('finished in ' +str(time.time()-t0))

    plt.imshow(output[1])

# %%
