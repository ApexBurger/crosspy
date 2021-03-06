# Subset compare is placed in separate file to allow cloudpickling in parallel

from crosspy.ImagePreparation import get_subset
from crosspy.hs import *
import numexpr as ne
import numpy as np

def subset_compare(subset_n,d,imnos,prepared_ffts,hs=False,cc_t=0.,cormeth="efficient"):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements 
    if hs == True:
        #subsets = np.stack((ref,test),axis=2)
        #results = correlate_subsets(ref,test,d,prepared_ffts)
        results = minimise_rt_lstsq(ref,test,cc_t, d, prepared_ffts,cormeth)
        return results
    else:
        dxs,dys,phs= fxcorr(ref,test,d,prepared_ffts,cormeth)
        return dxs,dys,phs