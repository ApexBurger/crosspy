# Following Britton et al ; Bergsmo & McAuliffe 2020

import numpy as np 
import numpy.fft 
import multiprocessing
import pyfftw 
import numexpr as ne
from crosspy.ImagePreparation import get_subset
import cv2 as cv

def plan_ffts(d,ffttype='fftw_numpy'):
    
    subset_shape=d.roi[0]

    a = pyfftw.empty_aligned((subset_shape, subset_shape), dtype='complex128')
    b = pyfftw.empty_aligned((subset_shape, subset_shape), dtype='complex128')
    c = pyfftw.empty_aligned((2*subset_shape, 2*subset_shape), dtype='complex128')
    d = pyfftw.empty_aligned((2*subset_shape, 2*subset_shape), dtype='complex128')

    # Over both axes
    forward_fft = pyfftw.FFTW(a, b, axes=(0,1))
    inverse_fft = pyfftw.FFTW(a,b ,axes=(0,1),direction='FFTW_BACKWARD')
    forward_xcf_fft1 = pyfftw.FFTW(c,d, axes=(0,1))
    inverse_xcf_fft1 = pyfftw.FFTW(c,d, axes=(0,1),direction='FFTW_BACKWARD')

    prepared_ffts=[forward_fft,inverse_fft,forward_xcf_fft1,inverse_xcf_fft1]

    if ffttype=='fftw_numpy':
        prepared_ffts=[pyfftw.interfaces.numpy_fft.fft2,pyfftw.interfaces.numpy_fft.ifft2,pyfftw.interfaces.numpy_fft.fft2,pyfftw.interfaces.numpy_fft.ifft2]
    elif ffttype=='fftw_scipy':
        prepared_ffts=[pyfftw.interfaces.scipy_fftpack.fft2,pyfftw.interfaces.scipy_fftpack.ifft2,pyfftw.interfaces.scipy_fftpack.fft2,pyfftw.interfaces.scipy_fftpack.ifft2]
    else:
        prepared_ffts=[np.fft.fft2,np.fft.ifft2,np.fft.fft2,np.fft.ifft2]

    return prepared_ffts


def gen_filters(roi, filter_settings=[4,2,16,32]):
    #genearte FFT filters

    # Inputs:
    #   size_pass = subset size - (128, 256 etc.)
    #   fpasset = [high pass cut off, high pass width, low pass cut off, low pass width]
    # outputs:
    #   fftfilter = non idea (gaussian) band pass filter
    #   hfilter = hanning filter

    size_pass=roi[0]

    pi = np.pi
    cos = np.cos
    dot = np.dot
    sqrt = np.sqrt
    exp = np.exp

    lcutoff = filter_settings[2]
    lwidth = filter_settings[3]/2
    hcutoff = filter_settings[0]
    hwidth = filter_settings[1]/2

    if lcutoff < hcutoff:
        print('low pass filter smaller than high pass filter')
        return

    # generate square grid
    u = range(0, size_pass)

    # meshgrid function
    meshv, meshu = np.meshgrid(u,u)

    meshvf = meshv-size_pass/2+0.5
    meshuf = meshu-size_pass/2+0.5

    # create Hann window
    hfilter = (cos(((pi*(meshuf)/size_pass)))*(cos((pi*meshvf/size_pass))))

    # create fft filter
    distf = sqrt((meshvf*meshvf)+(meshuf*meshuf))

    # lowpass
    lfftfilter = exp(-((distf-lcutoff)/(sqrt(2)*lwidth/2))**2)
    lfftfilter[distf > (lcutoff+2*lwidth)] = 0
    lfftfilter[distf < lcutoff] = 1


    # highpass
    hfftfilter = exp(-((hcutoff-distf)/(sqrt(2)*hwidth/2))**2)
    hfftfilter[distf < (hcutoff-2*hwidth)] = 0
    hfftfilter[distf > hcutoff] = 1

    # combine
    fftfilter = hfftfilter*lfftfilter
    fftfilter = np.fft.fftshift(fftfilter)
    
    return fftfilter, hfilter


def freg(ROI_test,ROI_ref,XCF_roisize,XCF_mesh,data_fill,prepared_ffts):
    
    #FREG Register two FFTs to subpixel accuracy
    #a reduced form of the code submitted to the matlab file exchange
    #http://www.mathworks.co.uk/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation
    
    #reported in the literature:
    #Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup
    #"Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
    
    #modified to handle the filtered FFT sizes
    #TBB 2012
    
    #ported to python
    #TPM 2020

    forward_fft=prepared_ffts[0]
    inverse_fft=prepared_ffts[1]
    forward_xcf_fft1=prepared_ffts[2]
    inverse_xcf_fft1=prepared_ffts[3]
    
    CC=np.zeros((XCF_roisize*2,XCF_roisize*2),dtype='complex')
    red_roisize=len(data_fill)

    ind1=int(XCF_roisize-np.floor(red_roisize/2))
    ind2=int(XCF_roisize+np.floor((red_roisize-1)/2))+1
    
    #fft shifting required to make it cross-correlation rather than convolution
    f=np.fft.fftshift(ROI_test*np.conj(ROI_ref))
    CC[ind1:ind2,ind1:ind2]= f

    #inverse FFT
    CC=inverse_xcf_fft1(CC)
    
    #get the maximum of the XCF and its location
    chi=np.amax(CC)
    loc=np.argmax(CC)
    rloc,cloc=np.unravel_index(loc,CC.shape)

    #get the shift in the original pixel grid from position of XCF peak
    XCF_roisize2=2*XCF_roisize

    if rloc > XCF_roisize:
        row_shift = rloc - XCF_roisize2
    else:
        row_shift = rloc

    if cloc > XCF_roisize:
        col_shift = cloc - XCF_roisize2
    else:
        col_shift = cloc
    
    row_shift = row_shift/2
    col_shift = col_shift/2

    #DFT computation
    row_shift = round(row_shift*XCF_mesh)/XCF_mesh
    col_shift = round(col_shift*XCF_mesh)/XCF_mesh
    dftshift = np.floor(np.ceil(XCF_mesh*1.5)/2)

    #matrix multiply DFT around current shift estimate
    roff=dftshift-row_shift*XCF_mesh
    coff=dftshift-col_shift*XCF_mesh

    #computer kernels and obtain DFT by matrix products
    prefac=-1j*2*np.pi/(XCF_roisize*XCF_mesh)

    #speed up kernel generation for reduced filtered FFT
    c1=[i for i in range(0,(XCF_roisize))]
    c_i=np.fft.ifftshift(c1).T - np.floor(XCF_roisize/2)
    c_i=c_i[data_fill]
    c_i=np.expand_dims(c_i,1)

    r_i=np.fft.ifftshift(c1) - np.floor(XCF_roisize/2)
    r_i=r_i[data_fill]
    r_i=np.expand_dims(r_i,0)

    m=np.array([i for i in range(0,int(np.ceil(XCF_mesh*1.5)))])
    m1=np.expand_dims(m,0)
    m2=np.expand_dims(m,1)

    arg1=(c_i)@(m1-coff)
    arg2=(m2-roff)@(r_i)
    
    kernc=ne.evaluate('exp(prefac*arg1)')
    kernr=ne.evaluate('exp(prefac*arg2)')
    kern = ROI_ref*ne.evaluate('conj(ROI_test)')
    arg3=kernr@kern@kernc
    CC2 = ne.evaluate('conj(arg3)')

    #locate maximum and map back to original pixel grid
    CCmax=np.abs(np.amax(CC2))
    loc=np.argmax(CC2)
    loc1,loc2=np.unravel_index(loc,CC2.shape)

    rloc=loc1-dftshift-1
    cloc=loc2-dftshift-1

    row_shift=row_shift+rloc/XCF_mesh
    col_shift=col_shift+cloc/XCF_mesh

    #bf1=np.sum(ROI_test.flatten()*np.conjugate(ROI_test.flatten()))
    #bf2=np.sum(ROI_ref.flatten()*np.conjugate(ROI_ref.flatten()))

    return col_shift, row_shift, CCmax#/np.sqrt(float(bf1*bf2))

def fxcorr(subset1,subset2,d,prepared_ffts, cormeth="efficient"):
    if cormeth == "efficient":
        forward_fft=prepared_ffts[0]
        inverse_fft=prepared_ffts[1]

        roi=d.roi    
        fftfil=d.fftfilter
        hfil=d.hfilter
        filters_settings=d.filter_settings

        #h-filter the subsets and generate the fft filter
        hfil = 1. # THIS IS AN ATTEMPT TO REMOVE H FILTER FOR HEAVISIDE
        subset1_filt=hfil*subset1
        subset2_filt=hfil*subset2
        

        #FFT the subsets
        f_s1=forward_fft(subset1_filt)
        f_s2=forward_fft(subset2_filt)

        fill1=(filters_settings[2]+filters_settings[3])
        fill2=(roi[0]-(filters_settings[2]+filters_settings[3]-1))
        fill3=roi[0]
        data_fill=np.array([i for i in range(fill1)]+[i for i in range(fill2-1,fill3)])

        #grab the reduced filter
        fftfil_red=fftfil[data_fill,:]
        fftfil_red=fftfil_red[:,data_fill]

        #grab the reduced subsets
        f_s1_red=f_s1[data_fill,:]
        f_s1_red=f_s1_red[:,data_fill]

        f_s2_red=f_s2[data_fill,:]
        f_s2_red=f_s2_red[:,data_fill]

        #perform the filtering (point-wise multiplication)
        ROI_ref=f_s1_red*fftfil_red
        ROI_test=f_s2_red*fftfil_red

        col_shift, row_shift, ccmax = freg(ROI_test,ROI_ref,roi[0],roi[2],data_fill,prepared_ffts)

    elif cormeth=="cv":
        """ This function uses openCV to upsample and correlate subsets
        It is fast but the subpixel accuracy is inefficient as upsampling
        is done in real space over the entire subset
        
        To improve, look at what ncorr does...
        or, write a custom function which finds the interger shift location, then performs
        a "mini" cross correlation on a smaller subset, upsampled
        """
        # upsampling factor
        f = d.roi[2]/100

        # upsample subsets by bicubic interpolation
        ref_up = interp_subset_cv(subset1,f).astype(np.float32)
        test_up = interp_subset_cv(subset2,f).astype(np.float32)
        # register shifts
        col_shift,row_shift,ccmax,_ = reg_cv(ref_up,test_up,f, method='cv.TM_CCOEFF_NORMED')
    elif cormeth == "pixelic":
        # register shifts
        col_shift,row_shift,ccmax,_ = reg_cv(subset1,subset2,f=1.0, method='cv.TM_CCOEFF_NORMED')

    return col_shift, row_shift, ccmax


def reg_cv(a,b,f, method):
    """registers shifts from correlogram
    """
    res = correlogram_cv(a, b, method)
    
    if method == 'cv.TM_SQDIFF_NORMED' or method == 'cv.TM_SQDIFF':
        cc = np.abs(np.amin(res))
        loc=np.argmin(res)
        loc1,loc2=np.unravel_index(loc, res.shape)

        row_shift=loc1-a.shape[0]
        col_shift=loc2-a.shape[0]
    else:
        cc = np.abs(np.max(res))
        loc=np.argmax(res)
        loc1,loc2=np.unravel_index(loc, res.shape)

        row_shift=loc1-a.shape[0]
        col_shift=loc2-a.shape[0]

    if abs(col_shift/f) > a.shape[0] or abs(row_shift/f) > a.shape[1]:
        col_shift = 0
        row_shift = 0
        print("Subset shift greater than subset size, setting to 0")
    
    return -col_shift/f, -row_shift/f, cc, res

def interp_subset_cv(array,factor):
    """
    Upsamples input array using bicubic interpolation
    """
    dimx = int(array.shape[1]*factor)
    dimy = int(array.shape[0]*factor)
    image = cv.resize(array.astype(np.float32), ( dimx, dimy ), interpolation = cv.INTER_CUBIC )
    return image

def correlogram_cv(a, b, method):
    """
        Returns correlogram of a and b arrays using openCV library
        This library is extremely efficient.
        
    """
    # zero pad reference image 
    img = np.pad(a, a.shape[0])
    img2 = img.copy()
    template = b
    w, h = template.shape[::-1]
    img = img2.copy()
    meth = eval(method)

    # Obtain correlogram through square-difference normed
    return cv.matchTemplate(img,template,meth)