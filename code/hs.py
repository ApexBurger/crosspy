# Heaviside treatment functions

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy.XCF import fxcorr
import os
import cv2 as cv
from numba import double, jit, njit, vectorize
from numba import int32, float32, uint8, float64, int64, boolean
from scipy.signal import fftconvolve
from numpy import cos, sin, radians

@njit#@vectorize(["boolean (float32, float32, float32, float32)"])
def line_test(x,y, r, theta):
    # function to test if a pixel lies above or below a line segment
    # Find start and end points of line segments for different cases of theta and r
    # Line segment to interogation point
    q1x = x
    q1y = y
    # Vector magnitude cases
    theta = theta % 360
    if r == 0:
        r = 1e-8
    # Rotation cases
    if theta == 0. or theta == 360.: # vertical to right
        x1 = r
        x2 = q1x
        if x2 > x1:
            return False
        else:
            return True
    elif theta == 90.: # horizontal line above
        y1 = r
        y2 = q1y
        if y2>y1:
            return False
        else:
            return True
    elif theta == 180.: # vertical to left
        x1 = -r
        x2 = q1x
        if x2 > x1:
            return True
        else:
            return False
    elif theta == 270.: # horizontal below
        y1 = -r
        y2 = q1y
        if y2 < y1:
            return False
        else:
            return True
    elif theta>0 and theta<180:
        theta = radians(theta)
        # Tangent line segment
        t1x = r*cos(theta)
        t1y = r*sin(theta)
        m = -1*(cos(theta)/sin(theta))
        c = t1y - m*t1x
        y1 = q1y
        y2 = m*q1x + c
        if y1>y2:
            return False
        else:
            return True
    elif theta>180 and theta<360:
        theta = radians(theta)
        # Tangent line segment
        t1x = r*cos(theta)
        t1y = r*sin(theta)
        m = -1*cos(theta)/sin(theta)
        c = t1y - m*t1x
        y1 = q1y
        y2 = m*q1x + c
        if y1<y2:
            return False
        else:
            return True

def gen_hsfilter(a,b, r, theta):
    # preallocate arrays
    subsets = np.stack((a,b),axis=2)
    filter_arr = np.zeros((subsets.shape[0],subsets.shape[0]), dtype=np.bool)
    xc = filter_arr.shape[0]/2
    yc = filter_arr.shape[1]/2
    xc = xc
    yc = yc
    r = np.array([r],dtype=np.float32)
    theta = np.array([theta],dtype=np.float32)
    # Create x and y coordinates which are centred
    x_length = np.linspace(-xc, xc,subsets.shape[0], dtype=np.float32)
    y_length = np.linspace(-yc,yc,subsets.shape[0], dtype=np.float32)
    xs,ys = np.meshgrid(x_length,y_length)
    hsfilter = filt_loop(subsets, r[0], theta[0], xs, ys, filter_arr)
    return hsfilter


@njit
def filt_loop(subsets, r, t, xs, ys, filter_arr):
    """ Loops quickly through every pixel indice and looks at position relative
        to discontinuity vector. If beyond line segment, then it returns false
        for that array value
    """
    xs = xs.astype(np.float32)
    ys = ys.astype(np.float32)
    flag = np.bool
    # iterate pixel by pixel
    for col in range(subsets.shape[0]):
        for row in range(subsets.shape[1]):
            #rasters through columns and rows for a given coordinate in xy
            # Note that y axis is mirrored
            x = xs[row, col]
            y = np.multiply(ys[row, col],-1)
            # Test if pixel is beyond the discontinuity line
            filter_arr[row,col]  = line_test(float32(x), float32(y), float32(r), float32(t))
                
    return filter_arr

def correlogram(a, b, method):
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

def reg(a,b, method = 'cv.TM_SQDIFF_NORMED'):
    """registers shifts from correlogram
    """
    res = correlogram(a, b, method)
    
    if method == 'cv.TM_SQDIFF_NORMED':
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

    if abs(col_shift) > a.shape[0]/2 or abs(row_shift) > a.shape[1]/2:
        col_shift = 0
        row_shift = 0
    
    return col_shift, row_shift, cc

@jit
def local_minima(array2d):
    """ Returns the local minima defined by local variation change.
        Square diff normalised typically filters noise extremely well.
        This is therefore only a suitable solution if your correlogram 
        does not have too many minor peaks.
    """
    local_minima = [ index 
                     for index in np.ndindex(array2d.shape)
                     if index[0] > 0
                     if index[1] > 0
                     if index[0] < array2d.shape[0] - 1
                     if index[1] < array2d.shape[1] - 1
                     if array2d[index] < array2d[index[0] - 1, index[1] - 1]
                     if array2d[index] < array2d[index[0] - 1, index[1]]
                     if array2d[index] < array2d[index[0] - 1, index[1] + 1]
                     if array2d[index] < array2d[index[0], index[1] - 1]
                     if array2d[index] < array2d[index[0], index[1] + 1]
                     if array2d[index] < array2d[index[0] + 1, index[1] - 1]
                     if array2d[index] < array2d[index[0] + 1, index[1]]
                     if array2d[index] < array2d[index[0] + 1, index[1] + 1]
                   ]
    return local_minima

def hs_corr(x, a, b, dic, prepared_ffts):
    # unpack vector
    r, theta = x
    # generate heaviside filter
    hsfilter = gen_hsfilter(a,b, r, theta)
    
    # assign to clean variables for legibility
    ah = a * hsfilter
    bh = b * hsfilter
    
    # invert selection
    hsfilter_inv = ~hsfilter

    # apply to clean variables
    ch = a * hsfilter_inv
    dh = b * hsfilter_inv

    results = fxcorr(ah, bh, dic, prepared_ffts)
    results_inv = fxcorr(ch, dh, dic, prepared_ffts)

    # return vectorised result containing dx,dy,cc
    return np.concatenate((results,results_inv),axis=None)

def hs_corr_norm(x, a, b, d, prepared_ffts, method='cv.TM_SQDIFF_NORMED'):
    # unpack vector 
    r,theta = x 

    #apply heaviside
    hsfilter = gen_hsfilter(a,b, r, theta)

    # apply to clean variables
    ah = a * hsfilter
    bh = b * hsfilter

    # get inverse
    hsfilter_inv = ~hsfilter

    # apply to clean variables
    ch = a * hsfilter_inv
    dh = b * hsfilter_inv
    # cross correlate using freg
    results = np.array(reg(ah,bh, method))
    results_inv = np.array(reg(ch,dh, method))

    return np.concatenate((results,results_inv),axis=None)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def disc_locate(a, b, mins, d, prepared_ffts):
    """ Attempts to locate discontinuity in subset
        returns r, theta and a fit measure
    """
    # Guess initial discontinuity angle based on shifts
    loc1, loc2 = [[i-a.shape[0] for i in j] for j in mins]
    y1,x1 = loc1
    y2,x2 = loc2
    jump = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    # find vector connecting peaks -> this should indicate initial angle of slip band given that slip band 
    v1 = np.array([x2-x1,y2-y1])
    v2 = np.array([1,0])
    ang = np.degrees(angle_between(v2,v1))
    # Search in r for given angle
    ang_pos = ang+90 % 360
    ang_neg = ang+270 % 360

    n = int(a.shape[0]/2)
    t1 = np.full(n,ang_pos)
    t2 = np.full(n,ang_neg)
    angles = np.concatenate((t1,t2))
    rs = np.linspace(0, n, n, dtype=np.int16)
    distan = np.concatenate((rs,rs))
    s_spac = np.c_[distan,angles]
    r_res = np.array([hs_corr_norm(i,a,b,d, prepared_ffts, 'cv.TM_CCOEFF_NORMED')[2] for i in s_spac])

    # return location of r, theta and a quality of fit parameter

    max_loc = np.argmax(r_res)
    r, theta = s_spac[max_loc,:]

    return r, theta, jump

def correlate_subsets(a, b, d, prepared_ffts):
    """ Returns dx, dy, cc, hson, r, theta, jump
    Locates theta based on vector drawn between peak locations in correlogram
    
    Works well for test data set but is not able to pick up multiple peaks in real datasets

    """ 
    # This library is extremely efficient.

    # pre-allocate output array
    out = np.full((7,),0, dtype=np.float64)

    # Obtain correlogram through square-difference normed
    res = correlogram(a, b,'cv.TM_SQDIFF_NORMED')

    # find local minima
    mins = local_minima(res)
    if len(mins) == 1:
        # No discontinuity -> try freg
        out[:3] = fxcorr(a,b,d,prepared_ffts)
    elif len(mins) == 2:
        # Single discontinuity -> try finding peak
        # mininmize
        x = disc_locate(a, b, mins, d, prepared_ffts)
        r, theta , jump = x
        dx,dy,cc,dxi,dyi,cci = hs_corr([r,theta], a, b, d, prepared_ffts)
        j = np.sqrt((dx-dxi)**2 + (dy-dyi)**2)
        out = np.array([dx,dy,cc,r,theta,True,j])
        # print(j/jump)
    elif len(mins) > 2:
        # Multiple discontinuities - TO DO
        mins = mins[:2]
        # Single discontinuity -> try finding peak
        # mininmize
        x = disc_locate(a, b, mins, d, prepared_ffts)
        r, theta, jump = x
        dx,dy,cc,dxi,dyi,cci = hs_corr([r,theta], a, b, d, prepared_ffts)
        j = np.sqrt((dx-dxi)**2 + (dy-dyi)**2)
        out = np.array([dx,dy,cc,r,theta,True,j])
    else:
        # No correlation found
        out = out

    return out

def disc_locate_angular(a, b, d, prepared_ffts):
    """ Searches in fine theta space to locate a discontinuity
        1. Searches for min CC value in circle near centre (notice sqdiff has good correlation for low values)
        2. Searches in a sector of theta_min +- 10 degrees for a range of r values
        3. Locates best fit of r and theta based on sum of CC and inverted mask CC values

        This method is effectively brute forcing which is optimised by considering the symmetry
        in reg(r,theta) over a subset
    """

    # search space
    t = np.linspace(1,360,36, dtype=np.float16)
    r = np.full(t.shape,a.shape[0]/10, dtype=np.float16)

    c1 = np.c_[r,t]

    res = np.array([hs_corr_norm(i,a,b,d, prepared_ffts, 'cv.TM_SQDIFF_NORMED')  \
        for i in c1]).reshape(36,6)
    dx,dy,cc,dxi,dyi,cci = [res[:,i] for i in range(6)]

    # find min
    cc_sum = cc+cci
    minval = np.min(cc_sum[np.nonzero(cc_sum)])
    minloc = np.where(cc_sum == minval)[0]
    r_sol, t_sol = c1[minloc[0], :]
    # Arc search
    t_min = (t_sol - 10) % 360
    t_max = (t_sol + 10) % 360
    t2 = np.linspace(t_min,t_max, 10, dtype=np.float16)
    r2 = np.linspace(0,a.shape[0]/4, 10, dtype=np.float16)

    r2, t2 = np.meshgrid(r2,t2)

    r2 = r2.flatten()
    t2 = t2.flatten()
    c2 = np.c_[r2,t2]
    res = np.array([hs_corr_norm(i,a,b,d, prepared_ffts, 'cv.TM_SQDIFF_NORMED') for i in c2]).reshape(10*10,6)
    dx,dy,cc,dxi,dyi,cci = [res[:,i] for i in range(6)]
    # find min
    cc_sum = cc+cci
    minval = np.min(cc_sum[np.nonzero(cc_sum)])
    minloc = np.where(cc_sum == minval)[0]
    if len(minloc)>1:
        minloc = minloc[0]
    minloc = int(minloc)
    # minimums
    r_sol, t_sol = c2[minloc, :]

    # find jump
    x0 = dx[minloc]
    x1 = dxi[minloc]
    y0 = dy[minloc]
    y1 = dyi[minloc]
    jump = np.sqrt((x1-x0)**2 + (y1-y0)**2)
    return r_sol, t_sol, jump

def minimise_rt_lstsq(a, b, cc_thresh, d, prepared_ffts):

    # obtain untreated dx,dy,
    dxu, dyu, ccu = reg(a, b, method="cv.TM_SQDIFF_NORMED")

    if ccu > cc_thresh: # Catch which excludes excellent correlations
        r, theta, jump = disc_locate_angular(a,b, d, prepared_ffts)
        x = np.array([r,theta])
        # obtain dx,dy,cc from r,theta
        dx,dy,cc,dxi,dyi,cci = hs_corr_norm(x,a,b,d,prepared_ffts)
        # compare values
        u = dx-dxi
        v = dy-dyi
        j = np.sqrt(u**2+v**2)
        if j >= 1:
            # If there is a jump accept treatment
            dx , dy, _, dxi, dyi, _ = hs_corr(x, a, b, d, prepared_ffts)
            r = float(r)
            theta = float(theta)
            out = [dx, dy, cc, r, theta, True, j]
        elif j < 1:
            # If no jump, reject treatment
            dx, dy, _ = fxcorr(a,b,d,prepared_ffts)
            out = [dx, dy, ccu, 0., 0., False, 0.]

        return out
    else:
        dxu,dyu,cc = fxcorr(a, b, d, prepared_ffts)
        out = [dxu,dyu,ccu,False,False,False,False]
        return out