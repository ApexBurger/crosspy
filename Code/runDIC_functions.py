import numpy as np
from imprep_functions import *
from XCF_functions import *
import matplotlib.pyplot as plt
from multiprocessing import Pool, Process, Queue
import functools

def subset_compare(d,imnos,subset_n):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements
    dxs,dys,phs=fxcorr(ref,test,d)
    print(str(subset_n))
    return dxs,dys,phs

def subset_compare_batch(d,imnos,batch):
    dx_batch=[]
    dy_batch=[]
    ph_batch=[]
    

    for subset_n in range(0,len(batch)):
        #grab the reference and test subsets, and get subpixel registration
        ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
        test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
        #get the displacements
        dx,dy,ph=fxcorr(ref,test,d)

        dx_batch.append(dx)
        dy_batch.append(dy)
        ph_batch.append(ph)
        print(subset_n)

    return dx_batch, dy_batch, ph_batch


def subset_compare_par(d,imnos,cores,chunk_length):

    #generate a to-do list
    inlist=[i for i in range(0,len(d.ss_locations))]

    # split the to-do list into batches
    def splitlist(inlist, chunk_length):
        remainder=len(inlist)%chunk_length
        n_jobs=int((len(inlist)-remainder)/chunk_length+1)
        result=[None for i in range(n_jobs)]

        for i in range(0,n_jobs):
            result[i]=inlist[i*chunk_length:(i+1)*chunk_length]
        result[n_jobs-1]=inlist[-remainder:]
        return result

    batches=splitlist(inlist,chunk_length)

    #make a partial function, as can only feed a single argument into the parallelised version
    subset_compare_batch_partial=functools.partial(subset_compare_batch,d,imnos)

    with Pool(cores) as p:
        output_batched = p.map(subset_compare_batch_partial,batches)

    #unpack the output
    output=[]
    for i in range(len(output_batched)):
        output=output+output_batched[i]

    return output
    

def run_DIC(d,imnos=[0,1],par=False,cores=5,chunk_length=50):
    
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
        output = subset_compare_par(d,imnos,cores,chunk_length)

        for subset_n in range(0,len(output)):
            dxs[subset_n]=output[subset_n][0]
            dys[subset_n]=output[subset_n][1]
            phs[subset_n]=output[subset_n][2]

    #translate best_dxs etc back onto image grid
    dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
    dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
    ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

    return dx_map,dy_map,ph_map