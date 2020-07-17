import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import functools
import pyfftw
from crosspy.ImagePreparation import *
from crosspy.hs import *
from crosspy.XCF import *
from crosspy.subset_compare import *
import numexpr as ne
from joblib import Parallel, delayed

def run_DIC(d,imnos=[0,1],hs=False, cores=None,ffttype='fftw_numpy'):
    #fft type can be : fftw_numpy (default), fftw_scipy, else defaults to numpy

    #set up numexpr to run with the chosen number of threads
    
    #discontinuity enables or disables slip trace tracking via heaviside filtering


    #preallocate for this DIC pair
    phs=np.zeros(d.n_subsets)
    dxs=np.zeros(d.n_subsets)
    dys=np.zeros(d.n_subsets)

    prepared_ffts=plan_ffts(d,ffttype)

    #enable the pyfftw cache for speedup
    pyfftw.interfaces.cache.enable()   
    
    # start mp pool

    #check for discontinuity tracker
    if hs == True:
        r = np.zeros(d.n_subsets)
        theta = np.zeros(d.n_subsets)
        hson = np.zeros(d.n_subsets)
        results = np.zeros((d.n_subsets,6))

        results = Parallel(n_jobs=cores, verbose=5)(delayed(subset_compare)(d=d,imnos=[0,1],subset_n=i,prepared_ffts=prepared_ffts,hs=True) for i in range(0,d.n_subsets))

        dxs,dys,phs,r,theta,hson = zip(*results)

        # convert to maps
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')
        r_map = np.reshape(r,(d.n_rows,d.n_cols),'F')
        theta_map = np.reshape(theta,(d.n_rows,d.n_cols),'F')
        hson_map = np.reshape(hson,(d.n_rows,d.n_cols),'F')
        return dx_map, dy_map, ph_map, r_map, theta_map, hson_map
    else:
        results = np.zeros((d.n_subsets,3))

        results = Parallel(n_jobs=cores, verbose=5)(delayed(subset_compare)(d=d, imnos=[0,1], subset_n=i, prepared_ffts=prepared_ffts) 
            for i in range(0,d.n_subsets))

        dxs,dys,phs = zip(*results)

        # convert to maps
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

        return dx_map,dy_map,ph_map