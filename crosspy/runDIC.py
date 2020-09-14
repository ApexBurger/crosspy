import numpy as np
import pyfftw
from crosspy.XCF import plan_ffts
from crosspy.subset_compare import subset_compare
from dask import delayed, compute
from dask.distributed import Client, LocalCluster

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

    cluster = LocalCluster(n_workers=cores, 
                       threads_per_worker=1,
                       memory_limit='16GB')
    client = Client(cluster)

    print(client)
    #check for discontinuity tracker
    if hs == True:
        r = np.zeros(d.n_subsets)
        theta = np.zeros(d.n_subsets)
        hson = np.zeros(d.n_subsets)
        u = np.zeros(d.n_subsets)
        v = np.zeros(d.n_subsets)
        j = np.zeros(d.n_subsets)

        #results = np.zeros((d.n_subsets,9))
        #results = np.zeros((d.n_subsets,9))
        results = []
        results_list = []
        # main loop
        for i in range(d.n_subsets):
            result = delayed(subset_compare)(d=d,imnos=[0,1],subset_n=i,prepared_ffts=prepared_ffts,hs=True)
            results.append(result)
        
        results_list = compute(*results)
        # for i in range(0,d.n_subsets):
        #     results[i,:] = subset_compare(d=d,imnos=[0,1],subset_n=i,prepared_ffts=prepared_ffts,hs=True)
        #     print(i)
        dxs,dys,phs,r,theta,hson, u, v, j = zip(*results_list)

        client.close()
        cluster.close()
        # convert to maps
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')
        r_map = np.reshape(r,(d.n_rows,d.n_cols),'F')
        theta_map = np.reshape(theta,(d.n_rows,d.n_cols),'F')
        hson_map = np.reshape(hson,(d.n_rows,d.n_cols),'F')
        u_map = np.reshape(u,(d.n_rows,d.n_cols),'F')
        v_map = np.reshape(v,(d.n_rows,d.n_cols),'F')
        j_map = np.reshape(j,(d.n_rows,d.n_cols),'F')

        return dx_map, dy_map, ph_map, r_map, theta_map, hson_map, j_map
    else:
        results = []
        results_list = []
        # main loop
        for i in range(d.n_subsets):
            result = delayed(subset_compare)(d=d,imnos=[0,1],subset_n=i,prepared_ffts=prepared_ffts,hs=False)
            results.append(result)
        
        results_list = compute(*results)
        
        dxs,dys,phs = zip(*results_list)

        client.close()
        cluster.close()
        # convert to maps
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

        return dx_map,dy_map,ph_map