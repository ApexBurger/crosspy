import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, Process, Queue
import functools
import crosspy
import multiprocessing

import partesting_XCF_2 as pxcf

import pyfftw

def plan_ffts(d):
    
    subset_shape=d.roi[0]
    #xcf_mesh=

    a = pyfftw.empty_aligned((subset_shape, subset_shape), dtype='complex128')
    b = pyfftw.empty_aligned((subset_shape, subset_shape), dtype='complex128')
    c = pyfftw.empty_aligned((2*subset_shape, 2*subset_shape), dtype='complex128')
    d = pyfftw.empty_aligned((2*subset_shape, 2*subset_shape), dtype='complex128')

    # Over the both axes
    forward_fft = pyfftw.FFTW(a, b, axes=(0,1))
    inverse_fft = pyfftw.FFTW(a,b ,axes=(0,1),direction='FFTW_BACKWARD')
    forward_xcf_fft1 = pyfftw.FFTW(c,d, axes=(0,1))
    inverse_xcf_fft1 = pyfftw.FFTW(c,d, axes=(0,1),direction='FFTW_BACKWARD')

    prepared_ffts=[forward_fft,inverse_fft,forward_xcf_fft1,inverse_xcf_fft1]

    return prepared_ffts

def subset_compare_1(d,imnos,subset_n,prepared_ffts):
    #grab the reference and test subsets, and get subpixel registration
    ref=crosspy.get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=crosspy.get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])

    ref=pyfftw.byte_align(ref)
    test=pyfftw.byte_align(test)

    #get the displacements
    dxs,dys,phs=pxcf.fxcorr(ref,test,d,prepared_ffts)
    #print(str(subset_n))
    return dxs,dys,phs