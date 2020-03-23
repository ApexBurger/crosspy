import numpy as np
from imprep_functions import *
from XCF_functions import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools

def subset_compare(d,imnos,subset_n):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements
    dxs,dys,phs=fxcorr(ref,test,d)
    return dxs,dys,phs

def subset_compare_par(d,imnos,chunks,cores):
    subset_compare_partial=functools.partial(subset_compare,d,imnos)

    with Pool(cores) as p:
        output = p.map(subset_compare_partial,range(0,len(d.ss_locations)),chunks)
    return output

def run_DIC(d,imnos=[0,1],par=False,chunks=10,cores=None):
    
    #preallocate for this DIC pair
    phs=np.zeros(d.n_subsets)
    dxs=np.zeros(d.n_subsets)
    dys=np.zeros(d.n_subsets)

    #serial version of this function
    if par==False:

        for subset_n in range(0,d.n_subsets):
            dxs[subset_n],dys[subset_n],phs[subset_n]=subset_compare(d,imnos,subset_n)

    #parallel version
    if par==True:
        output = subset_compare_par(d,imnos,chunks,cores)

        for subset_n in range(0,len(output)):
            dxs[subset_n]=output[subset_n][0]
            dys[subset_n]=output[subset_n][1]
            phs[subset_n]=output[subset_n][2]

    #translate best_dxs etc back onto image grid
    dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
    dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
    ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

    return dx_map,dy_map,ph_map