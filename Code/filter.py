# Filter functions

def dicfilter(roisize, fpassset):
    # Creates filters for cross correlation

    pi = np.pi
    cos = np.cos
    dot = np.dot
    sqrt = np.sqrt
    exp = np.exp

    lcutoff = fpasset[2]
    lwidth = fpassset[3]/2
    hcutoff = fpassset[0]
    hwidth = fpassset[1]/2

    if lcutoff < hcutoff:
        print('low pass filter smaller than high pass filter')
        return

    # generate square grid
    u = range(1, roisize)
    # meshgrid function
    meshv, meshu = meshgrid(u,u)
    meshvf = meshv-roisize/2-0.5
    meshuf = meshu-roisize/2-0.5


    # create Hfilter
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