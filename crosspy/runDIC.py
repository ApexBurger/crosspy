import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import functools
import pyfftw
from crosspy.XCF import *
from crosspy.ImagePreparation import *
from crosspy.hs import *
import numexpr as ne

#import crosspy

def subset_compare(d,imnos,subset_n,prepared_ffts,discontinuity=False):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements 
    if discontinuity == True:
        subsets = np.stack((ref,test),axis=2)
        dxs, dys, phs,r,theta, hson = minimise_rt_lstsq(subsets,d, prepared_ffts)
        return dxs,dys,phs,r,theta,hson
    else:
        dxs,dys,phs=fxcorr(ref,test,d,prepared_ffts)
        return dxs,dys,phs



# def subset_compare_batch(d,imnos,batch):

#     #this could be moved outside the function for speedup
#     dx_batch=np.empty((len(batch)))
#     dy_batch=np.empty((len(batch)))
#     ph_batch=np.empty((len(batch)))

#     for i in range(0,len(batch)):
#         subset_n=batch[i]

#         #grab the reference and test subsets, and get subpixel registration
#         ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
#         test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
#         #get the displacements
#         dx,dy,ph=fxcorr(ref,test,d,prepared_ffts)

#         #append these
#         dx_batch[i]=dx
#         dy_batch[i]=dy
#         ph_batch[i]=ph
#         #print(subset_n)

#     return (dx_batch,dy_batch,ph_batch)


# def subset_compare_par(d,imnos,cores,chunk_length):

#     #generate a to-do list
#     inlist=[i for i in range(0,len(d.ss_locations))]

#     # split the to-do list into batches
#     def splitlist(inlist, chunk_length):
#         remainder=len(inlist)%chunk_length
#         n_jobs=int((len(inlist)-remainder)/chunk_length+1)
#         result=[None for i in range(n_jobs)]

#         for i in range(0,n_jobs):
#             result[i]=inlist[i*chunk_length:(i+1)*chunk_length]
#         result[n_jobs-1]=inlist[-remainder:]
#         return result

#     batches=splitlist(inlist,chunk_length)

#     #make a partial function, as can only feed a single argument into the parallelised version
#     subset_compare_batch_partial=functools.partial(subset_compare_batch,d,imnos)

#     with Pool(cores) as p:
#         output_batched = p.map(subset_compare_batch_partial,batches)

#     #unpack the results
#     dxs_batched=[output_batched[i][0] for i in range(0,len(output_batched))]
#     dys_batched=[output_batched[i][1] for i in range(0,len(output_batched))]
#     phs_batched=[output_batched[i][2] for i in range(0,len(output_batched))]

#     dxs=list(dxs_batched[0][:])
#     dys=list(dys_batched[0][:])
#     phs=list(phs_batched[0][:])

#     for i in range(1,len(dxs_batched)):
#             dxs=dxs+list(dxs_batched[i][:])
#             dys=dys+list(dys_batched[i][:])
#             phs=phs+list(phs_batched[i][:])

#     return dxs,dys,phs
    

def run_DIC(d,imnos=[0,1], cores=None,ffttype='fftw_numpy', discontinuity=False):
    #fft type can be : fftw_numpy (default), fftw_scipy, else defaults to numpy

    #set up numexpr to run with the chosen number of threads
    
    #discontinuity enables or disables slip trace tracking via heaviside filtering
    if cores==None:
        cores=multiprocessing.cpu_count()

    ne.set_num_threads(cores)


    #preallocate for this DIC pair
    phs=np.zeros(d.n_subsets)
    dxs=np.zeros(d.n_subsets)
    dys=np.zeros(d.n_subsets)

    prepared_ffts=plan_ffts(d,ffttype)

    #decide the number of cores to allocate
    if cores:
        pyfftw.config.NUM_THREADS = cores
    else:
        pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()

    #enable the pyfftw cache for speedup
    pyfftw.interfaces.cache.enable()   
    
    #check for discontinuity tracker
    if discontinuity == True:
        r = np.zeros(d.n_subsets)
        theta = np.zeros(d.n_subsets)
        hson = np.zeros(d.n_subsets)
        for subset_n in range(0,d.n_subsets):
            dxs[subset_n],dys[subset_n],phs[subset_n],r[subset_n],theta[subset_n],hson[subset_n]=subset_compare(d,imnos,subset_n,prepared_ffts,discontinuity=True)
        
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')
        r_map = np.reshape(r,(d.n_rows,d.n_cols),'F')
        theta_map = np.reshape(theta,(d.n_rows,d.n_cols),'F')
        hson_map = np.reshape(hson,(d.n_rows,d.n_cols),'F')
        
        return dx_map, dy_map, ph_map, r_map, theta_map, hson_map
    else:
        for subset_n in range(0,d.n_subsets):
            dxs[subset_n],dys[subset_n],phs[subset_n]=subset_compare(d,imnos,subset_n,prepared_ffts)

        #translate best_dxs etc back onto image grid
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

        return dx_map,dy_map,ph_map