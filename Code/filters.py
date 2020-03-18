# Filters
# This file contains imaging filters for the main DIC script
# Current filters:
#   - Han window
#   - FFT filter
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
    u = range(0, roisize)
    # meshgrid function
    meshv, meshu = np.meshgrid(u,u)
    meshvf = meshv-roisize/2-0.5
    meshuf = meshu-roisize/2-0.5


    # create Hann window
    hfilter = (cos(((pi*(meshuf)/roisize)))*(cos((pi*meshvf/roisize))))

    # create fft filter

    distf = sqrt((meshvf*meshvf)+(meshuf*meshuf))

    # lowpass
    lfftfilter = exp(-((distf-lcutoff)/sqrt(2)*lwidth/2)**2)
    # need to iterate through the array to satisfy expression:
    lfftfilter[distf > (lcutoff+2*lwidth)] = 0
    lfftfilter[distf < (lcutoff+2*lwidth)] = 1
    # distf < lcutoff = 1

    # highpass

    hfftfilter = exp(-((hcutoff-distf)/sqrt(2)*hwidth/2))**2
    hfftfilter[distf < (hcutoff-2*hwidth)] = 0
    hfftfilter[distf > hcutoff] = 1

    # combine

    fftfilter = hfftfilter*lfftfilter
    fftfilter = np.fft.fftshift(fftfilter) #THIS IS A MATLAB FUNCTION - find python equivalent
    
    return fftfilter, hfilter
