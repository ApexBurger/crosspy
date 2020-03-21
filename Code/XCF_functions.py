# Following Britton et al ; Bergsmo & McAuliffe 2020

import numpy as np 
import numpy.fft 
from numba import jit

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

def freg(ROI_test,ROI_ref,XCF_roisize,XCF_mesh,data_fill):
    CC=np.zeros((XCF_roisize*2,XCF_roisize*2),dtype='complex')
    red_roisize=len(data_fill)

    ind1=int(XCF_roisize-np.floor(red_roisize/2))
    ind2=int(XCF_roisize+np.floor((red_roisize-1)/2))+1
    #fft shifting required to make it cross-correlation rather than convolution
    f=np.fft.fftshift(ROI_test*np.conj(ROI_ref))
    CC[ind1:ind2,ind1:ind2]= f

    #inverse FFT
    CC=np.fft.ifft2(CC)
    
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

    kernc=np.exp(prefac*(c_i)@(m1-coff))
    kernr=np.exp(prefac*(m2-roff)@(r_i))
    kern = ROI_ref*np.conj(ROI_test)
    CC2 = np.conj(kernr@kern@kernc)

    #locate maximum and map back to original pixel grid
    CCmax=np.abs(np.amax(CC2))
    loc=np.argmax(CC2)
    loc1,loc2=np.unravel_index(loc,CC2.shape)

    rloc=loc1-dftshift-1
    cloc=loc2-dftshift-1

    row_shift=row_shift+rloc/XCF_mesh
    col_shift=col_shift+cloc/XCF_mesh

    return col_shift, row_shift, CCmax

def fxcorr(subset1,subset2,d):

    roi=d.roi    
    fftfil=d.fftfilter
    hfil=d.hfilter
    filters_settings=d.filter_settings

    #h-filter the subsets and generate the fft filter
    subset1_filt=hfil*subset1
    subset2_filt=hfil*subset2

    #FFT the subsets
    f_s1=np.fft.fft2(subset1_filt)
    f_s2=np.fft.fft2(subset2_filt)

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

    col_shift, row_shift, CCmax = freg(ROI_test,ROI_ref,roi[0],roi[2],data_fill)
    return col_shift, row_shift, CCmax
