# Heaviside treatment functions

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy.XCF import fxcorr
import os
import cv2
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

def gen_hsfilter(subsets, r, theta):
    # preallocate arrays
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
    # cols = np.arange(subsets.shape[0])
    # rows = np.arange(subsets.shape[0])
    #print(ys.dtype)
    #filter_arr = np.empty((subsets.shape[0],subsets.shape[0]), dtype=np.bool)
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


#@jit(nopython=False)
def ncc(b, a, mode="full"):
    """
    b is deformed subset
    a is reference subset

    This is a matlab port following Ujash Joshi, University of Toronto, 2017
    """

    # If this happens, it is probably a mistake
    if np.ndim(b) > np.ndim(a) or \
            len([i for i in range(np.ndim(b)) if b.shape[i] > a.shape[i]]) > 0:
        print("normxcorr2: TEMPLATE larger than IMG. Arguments may be swapped.")

    template = b - np.mean(b)
    image = a - np.mean(a)

    a1 = np.ones(template.shape)
    # Faster to flip up down and left right then use fftconvolve instead of scipy's correlate
    ar = np.flipud(np.fliplr(template))
    out = fftconvolve(image, ar.conj(), mode=mode)
    image = fftconvolve(np.square(image), a1, mode=mode) - \
            np.square(fftconvolve(image, a1, mode=mode)) / (np.prod(template.shape))

    # Remove small machine precision errors after subtraction
    image[np.where(image < 0)] = 0

    template = np.sum(np.square(template))
    out = out / np.sqrt(image * template)

    # Remove any divisions by 0 or very close to 0
    out[np.where(np.logical_not(np.isfinite(out)))] = 0

    CCmax=np.abs(np.amax(out))
    loc=np.argmax(out)
    loc1,loc2=np.unravel_index(loc,out.shape)

    row_shift=loc1-a.shape[0]
    col_shift=loc2-a.shape[0]

    return col_shift, row_shift, CCmax



#@jit(nopython=False)
def hs_corr(x, subsets, dic, prepared_ffts):
    # unpack vector
    r, theta = x
    # generate heaviside filter
    hsfilter = gen_hsfilter(subsets, r, theta)
    
    # assign to clean variables for legibility
    a = subsets[:,:,0]#*hsfilter
    b = np.multiply(subsets[:,:,1],hsfilter)
    
    # cross correlate using freg - this is still faster than silly numpy correlate
    results = fxcorr(a, b, dic, prepared_ffts)
    # invert selection
    hsfilter_inv = ~hsfilter

    # apply to clean variables
    c = subsets[:,:,0]# * hsfilter_inv
    d = subsets[:,:,1] * hsfilter_inv

    results_inv = fxcorr(c, d, dic, prepared_ffts)
    
    # return vectorised result containing dx,dy,cc
    return np.concatenate((results,results_inv),axis=None)

#@jit(nopython=False)
def hs_corr_norm(x,subsets,d,prepared_ffts):
    # unpack vector
    r,theta = x

    #apply heaviside
    hsfilter = gen_hsfilter(subsets, r, theta)
    # apply to clean variables
    a = subsets[:,:,0]# * hsfilter
    b = np.multiply(subsets[:,:,1],hsfilter)

    # get inverse
    hsfilter_inv = ~hsfilter

    # apply to clean variables
    c = subsets[:,:,0]# * hsfilter_inv
    d = subsets[:,:,1] * hsfilter_inv

    # cross correlate using freg
    results = np.array(ncc(b,a))
    results_inv = np.array(ncc(d,c))
    
    results = np.concatenate((results,results_inv),axis=None)
    return results

def A_mat(x,y=0,n=2):
    # Function to generate A matrices for polynomial solutions of Ax = b
    if y == 0:
        A = np.zeros((len(x),6))
        A[:,0]=np.squeeze(x**2)
        A[:,1]=np.squeeze(x)
        A[:,2]=np.ones(len(x))
        return A 
    elif n == 2:
        A = np.zeros((len(x),6))
        A[:,0]=np.squeeze(x**2)
        A[:,1]=np.squeeze(x)
        A[:,2]=np.squeeze(y**2)
        A[:,3]=np.squeeze(y)
        A[:,4]=np.squeeze(x*y)
        A[:,5]=np.ones(len(x))
        return A
    elif n == 3:
        A = np.zeros((len(x),10))
        A[:,0]=np.squeeze(x**3)
        A[:,1]=np.squeeze((x**2)*y)
        A[:,2]=np.squeeze(x*(y**2))
        A[:,3]=np.squeeze(y**3)
        A[:,4]=np.squeeze(x**2)
        A[:,5]=np.squeeze(x*y)
        A[:,6]=np.squeeze(y**2)
        A[:,7]=np.squeeze(x)
        A[:,8]=np.squeeze(y)
        A[:,9]=np.ones(len(x))
        return A
    elif n == 4:
        A = np.zeros((len(x),15))
        A[:,0]=np.squeeze(x**4)
        A[:,1]=np.squeeze(y**4)
        A[:,2]=np.squeeze((x**3)*y)
        A[:,3]=np.squeeze((x**2) * (y**2))
        A[:,4]=np.squeeze(x*(y**3))
        A[:,5]=np.squeeze(x**3)
        A[:,6]=np.squeeze(y**3)
        A[:,7]=np.squeeze((x**2)*y)
        A[:,8]=np.squeeze(x*(y**2))
        A[:,9]=np.squeeze(x**2)
        A[:,10]=np.squeeze(y**2)
        A[:,11]=np.squeeze(x*y)
        A[:,12]=np.squeeze(x)
        A[:,13]=np.squeeze(y)
        A[:,14]=np.ones(len(x))
        return A

# @jit(nopython=False)
def A_mat(x):
    A = np.zeros((len(x),6))
    A[:,0]=np.squeeze(x**4)
    A[:,1]=np.squeeze(x**3)
    A[:,2]=np.squeeze(x**2)
    A[:,3]=np.squeeze(x)
    A[:,4]=np.ones(len(x))
    return A

# @jit(nopython=False)
def hsin(x,subsets,d,prepared_ffts):
    _,_,c,_,_,_=hs_corr_norm(x,subsets,d,prepared_ffts)
    return c

def A_mat_sin(x):
    A = np.zeros((len(x),3))
    A[:,0] = np.ones(len(x)) # a0
    A[:,1] = np.squeeze(np.sin(x)) #a1
    A[:,2] = np.squeeze(np.cos(x)) #b0
    return A
#@jit(nopython=False)
def locate_disc(subsets, d, prepared_ffts):
    
    # first find theta
    theta = np.linspace(0, 360, 36)
    theta_mod = np.linspace(0,360,10000)
    r = np.full(len(theta), subsets.shape[1]/100)

    # generate trial Ax=b based on measured response
    atr = A_mat(theta)
    vals = np.c_[r,theta]

    # here we measure XCF
    
    fx = np.array([hsin(val,subsets, d, prepared_ffts) for val in vals])
    
    # fit a sin function
    
    fx = fx-np.mean(fx)

    trial_func = lambda c,x : c[0] + c[1]*np.sin(x) + c[2]*np.cos(x)



    th = np.radians(theta)*2
    A = A_mat_sin(th)

    x, _, _, _ = np.linalg.lstsq(A, fx)

    th = np.radians(theta_mod)*2
    fx_m = trial_func(x,th)

    # evaluate fit
    th_f = np.radians(36)*2
    y = fx
    yhat = trial_func(x,theta)
    ybar = np.sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((y-ybar)**2)

    rsq = ssreg/sstot

    if np.isnan(rsq):
        return False,False,False
    # plot for the sake of it

    # fig, ax = plt.subplots()
    # ax.scatter(theta,fx)
    # ax.plot(np.degrees(th)/2, fx_m)
    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    # plt.show()

    loc = np.where(fx_m.max() == fx_m)
    
    th_mx = th[loc]
  
    # catch dual peaks
    if len(th_mx)>1:
        th_mx = th_mx[0]
    elif not th_mx:
        return False

    theta_est = int(np.degrees(th_mx/2))
    theta_est2 = int((theta_est+180) % 360)

    c1 = hsin(np.array([r[0],theta_est]),subsets, d, prepared_ffts)
    c2 = hsin(np.array([r[0],theta_est2]),subsets, d, prepared_ffts)

    if c1>c2:
        theta =  np.uintc(theta_est)
    else:
        theta =  np.uintc(theta_est2)
    # angle space
    theta_arr = np.full((subsets.shape[1]),theta)

    # create search space in r
    r_max = subsets.shape[1]/2
    r = np.linspace(0,r_max,subsets.shape[1])

    c = np.column_stack((r,theta_arr))

    # search
    fx = np.array([hsin(val,subsets, d, prepared_ffts) for val in c])

    # measure
    loc = np.where(fx.max() == fx)
    r_sol = r[loc]

    return int(r_sol[0]), theta, rsq

def minimise_rt_lstsq(subsets, d, prepared_ffts):

    # obtain peak r and theta
    r, theta, rsq = locate_disc(subsets, d, prepared_ffts)
    if rsq == False:
        return np.array([False,False,False,False,False,False,False,False,False])
    x = np.array([r,theta])
    # obtain dx,dy,cc from r,theta
    dx,dy,cc,dxi,dyi,cci = hs_corr_norm(x,subsets,d,prepared_ffts)
    # obtain untreated dx,dy,
    dxu, dyu, ccu = ncc(subsets[:,:,1], subsets[:,:,0])
    # compare values
    u = dx-dxi
    v = dy-dyi
    if u == 0 and v == 0: # both dx and dy are in the same direction
        # Both in same direction => no discontinuity
        dxu,dyu,ccu = fxcorr(subsets[:,:,0],subsets[:,:,1], d, prepared_ffts)
        return np.array([dxu,dyu,ccu,False,False,False,False,False,False])
    elif (u != 0 or v != 0) and (rsq > 0.7): # dx or dy are not the same - presence of discontinuity
        # One value does not match direction of displacment - indicating a kinematic shift
        dx,dy,cc,dxi,dyi,cci = hs_corr(x, subsets, d, prepared_ffts)
        u = dx - dxi
        v = dy - dyi
        j = np.sqrt(u**2+v**2)
        return np.array([dx,dy,cc,r,theta,True,u,v,j])
    else:
        # fit is bad indicating no symmetry
        dxu,dyu,ccu = fxcorr(subsets[:,:,0],subsets[:,:,1], d, prepared_ffts)
        return np.array([dxu,dyu,ccu,False,False,False,False,False,False])