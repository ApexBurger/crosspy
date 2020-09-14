import numpy as np
import pyfftw
from crosspy.XCF import plan_ffts
from crosspy.subset_compare import subset_compare
from joblib import delayed, Parallel

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

    #check for discontinuity tracker
    if hs == True:
        r = np.zeros(d.n_subsets)
        theta = np.zeros(d.n_subsets)
        hson = np.zeros(d.n_subsets)
        u = np.zeros(d.n_subsets)
        v = np.zeros(d.n_subsets)
        j = np.zeros(d.n_subsets)

        results = np.zeros((d.n_subsets,9))

        # main loop
        # for i in range(d.n_subsets):
        #     result = delayed(subset_compare)(d=d,imnos=[0,1],subset_n=i,prepared_ffts=prepared_ffts,hs=True)
        #     results.append(result)
        
        results = Parallel(n_jobs=cores, verbose=5)(delayed(subset_compare) \
            (subset_n=i, d=d, imnos=[0,1], prepared_ffts=prepared_ffts, hs=True) for i in range(d.n_subsets))
        
        dxs,dys,phs,r,theta,hson, u, v, j = zip(*results)

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

        return dx_map, dy_map, ph_map, r_map, theta_map, hson_map, u_map, v_map, j_map
    else:
        results = np.empty(shape=(d.n_subsets,9))
        # main loop
        results = Parallel(n_jobs=cores, verbose=5)(delayed(subset_compare) \
            (subset_n=i, d=d, imnos=[0,1], prepared_ffts=prepared_ffts, hs=False) for i in range(d.n_subsets))

        dxs,dys,phs = zip(*results)

        # convert to maps
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

        return dx_map,dy_map,ph_map