import numpy as np 
import numpy.fft 

def freg(ROI_test,ROI_ref,XCF_roisize,XCF_mesh,data_fill):
    
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