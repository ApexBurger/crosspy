import numpy as np

def gen_filters(size_pass, filter_settings=[4,2,16,32]):
    # Inputs:
    #   size_pass = subset size - (128, 256 etc.)
    #   fpasset = [high pass cut off, high pass width, low pass cut off, low pass width]
    # outputs:
    #   fftfilter = non idea (gaussian) band pass filter
    #   hfilter = hanning filter

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
    meshvf = meshv-size_pass/2-0.5
    meshuf = meshu-size_pass/2-0.5

    # create Hann window
    hfilter = (cos(((pi*(meshuf)/size_pass)))*(cos((pi*meshvf/size_pass))))

    # create fft filter
    distf = sqrt((meshvf*meshvf)+(meshuf*meshuf))

    # lowpass
    lfftfilter = exp(-((distf-lcutoff)/sqrt(2)*lwidth/2)**2)
    lfftfilter[distf > (lcutoff+2*lwidth)] = 0
    lfftfilter[distf < (lcutoff+2*lwidth)] = 1


    # highpass
    hfftfilter = exp(-((hcutoff-distf)/sqrt(2)*hwidth/2))**2
    hfftfilter[distf < (hcutoff-2*hwidth)] = 0
    hfftfilter[distf > hcutoff] = 1

    # combine
    fftfilter = hfftfilter*lfftfilter
    fftfilter = np.fft.fftshift(fftfilter)
    
    return fftfilter, hfilter

def gen_ROIs(imshape,size_pass,overlap):
    
    rows=imshape[0]
    cols=imshape[1]
    spacing = round(size_pass*(1-overlap))
    
    #size_pass = subset size = 128,256 etc, overlap = proportion, eg 0.5, image1 and image2 must be same shape
    col_remainder=cols%spacing
    row_remainder=rows%spacing

    n_col_sets=int((cols-col_remainder)/spacing) #number of subsets in horizontal direction
    n_row_sets=int((rows-row_remainder)/spacing) #number of subsets in vertical direction

    #generate an array of subset locations [top L corner row, top L corner col] - can get bot R row and col from known subset size
    ss_locations=np.zeros((n_row_sets*n_col_sets,2))
    for c in range(0,n_col_sets):
        for r in range(0,n_row_sets):
            i = c*n_row_sets+r #index in this list
            #print(c,r)
            ss_locations[i,:]=np.array([r*spacing,c*spacing])
    
    return ss_locations


