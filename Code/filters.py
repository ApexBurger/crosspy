# Filters
# This file contains imaging filters for the main DIC script
# Written by Alexander Bergsmo 2020
# Imperial College London 
# ab12215@ic.ac.uk
import numpy as np

def xfilters(roisize, fpassset):
    # Inputs:
    #   roisize = subset size - (128, 256 etc.)
    #   fpasset = [high pass cut off, high pass width, low pass cut off, low pass width]
    # outputs:
    #   fftfilter = non idea (gaussian) band pass filter
    #   hfilter = hanning filter
    pi = np.pi
    cos = np.cos
    dot = np.dot
    sqrt = np.sqrt
    exp = np.exp

    lcutoff = fpassset[2]
    lwidth = fpassset[3]/2
    hcutoff = fpassset[0]
    hwidth = fpassset[1]/2

    if lcutoff < hcutoff:
        print('low pass filter smaller than high pass filter')
        return

    # generate square grid
    u = range(1, roisize)
    # meshgrid function
    meshv, meshu = np.meshgrid(u,u)
    meshvf = meshv-roisize/2-0.5
    meshuf = meshu-roisize/2-0.5


    # create Hanning Filter
    hfilter = (cos(((pi*(meshuf)/roisize)))*(cos((pi*meshvf/roisize))))

    # create fft filter

    distf = sqrt(dot(meshvf*meshvf)+(meshuf*meshuf))

    # lowpass
    lfftfilter = exp(-((distf-lcutoff)/sqrt(2)*lwidth/2))**2)
    # need to iterate through the array to satisfy expression:
    # distf > (lcufoff+2*lwidth) = 0 and
    # distf < lcutoff = 1

    # highpass

    hfftfilter = exp(-((hcutoff-distf)/sqrt(2)*hwidth/2))**2)
    # need to iterate through array to satisfy
    # distf < (hcutoff-2*hwidth) = 2 and
    # distf > hcutoff = 1

    # combine

    fftfilter = hfftfilter*lfftfilter
    fftfilter = fftshift(fftfilter) #THIS IS A MATLAB FUNCTION - find python equivalent
    
    return fftfilter, hfilter
    return fftfilter, hfilter