
@jit()
def freg_red(roi_test, roi_ref, xcf_roisize, xcf_mesh, datafill, prepared_ffts):
    # reduced FREG with no subpixel registration

    # ROI_test = test sub image
    # ROI_ref = reference sub image
    # XCF_roisize = size of sub image
    # XCF_mesh = this determines the subpixel registration - we will do away with this
    # datafill = filtering array
    # prepared_ffts = prepared fft arrays for FFTW

    # assign prepared ffts arrays
    forward_fft=prepared_ffts[0]
    inverse_fft=prepared_ffts[1]
    forward_xcf_fft1=prepared_ffts[2]
    inverse_xcf_fft1=prepared_ffts[3]

    # initialise arrays
    cc = np.zeros((xcf_roisize*2,xcf_roisize*2), dtype="complex")
    red_roisize = len(datafill) # reduced

    # indices
    ind1 = int(xcf_roisize-np.floor(red_roisize/2))
    ind2 = int(xcf_roisize+np.floor(red_roisize-1/2))+1

    # FFT shifting to enable cross-correlation
    f = np.fft.fftshift(roi_test*np.conj(roi_ref))

    # assign centre according to the reduced roisize
    cc[ind1:ind2,ind1:ind2] = f 

    # inverse fft
    cc = inverse_xcf_fft1(cc)

    # get max of xcf and location
    chi = np.amax(cc)
    loc = np.argmax(cc)
    rloc,cloc = np.unravel_index(loc,cc.shape)

    # get shift in original pixel grid from position of xcf peak

    xcf_roisize2 = 2*xcf_roisize

    if rloc > xcf_roisize:
        row_shift = rloc - xcf_roisize2
    else:
        row_shift = rloc
    
    if cloc > xcf_roisize:
        col_shift = cloc - xcf_roisize2
    else:
        col_shift = cloc

    row_shift = row_shift/2
    col_shift = col_shift/2

    # Discrete fourier transform

    row_shift = round(row_shift*xcf_mesh)/xcf_mesh # row shift as interpreted by subpixel mesh
    col_shift = round(col_shift*xcf_mesh)/xcf_mesh
    dft_shift = np.floor(np.ceil(XCF_mesh*1.5)/2)

    # matrix multiply dft around current shift estimate
    roff = dft_shift-row_shift*xcf_mesh
    coff = dft_shift-col_shift*xcf_mesh

    # s
    c1 = [i for i in range(0,(xcf_roisize))]
    c_i = np.fft.ifftshift(c1).T - np.floor(xcf_roisize/2)
    # obtain filtered result
    c_i = c_i[datafill]
    # expand to obtain conforming results
    c_i = np.expand_dims(c_i,1)

    r_i = np.fft.ifftshift(c1) - np.floor(xcf_roisize/2)
    r_i = r_i[datafill]
    r_i = np.expand_dims(c_i,1)

    m = np.array([i for i in range(0,int(np.ceil(xcf_mesh*1.5)))])
    m1 = np.expand_dims(m,0)
    m2 = np.expand_dims(m,1)

    arg1 = (c_i)@(m1-coff)
    arg2 = (m2-roff)@(r_i)

    kernc = ne.evaluate('exp(prefac*arg1)')
    kernr = ne.evaluate('exp(prefac*arg2)')
    kern = roi_ref * ne.evaluate('conj(roi_test)')
    arg3 = kernr@kern@kernc