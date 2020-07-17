# Subset compare is placed in separate file to allow cloudpickling in parallel

import functools
import pyfftw
from crosspy.XCF import *
from crosspy.ImagePreparation import *
from crosspy.hs import *
import numexpr as ne


def subset_compare(d,imnos,subset_n,prepared_ffts,hs=False):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements 
    if hs == True:
        subsets = np.stack((ref,test),axis=2)
        dxs, dys, phs,r,theta, hson = minimise_rt_lstsq(subsets,d, prepared_ffts)
        return dxs,dys,phs,r,theta,hson
    else:
        dxs,dys,phs=fxcorr(ref,test,d,prepared_ffts)
        return dxs,dys,phs