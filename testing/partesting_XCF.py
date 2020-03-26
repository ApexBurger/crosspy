import numpy as np
from freg import freg

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