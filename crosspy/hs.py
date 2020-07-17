# Heaviside treatment functions

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy import DIC, Imset
from crosspy.XCF import fxcorr
import os
import cv2
from numba import double, jit, int32, float32, uint8, jitclass

@jitclass([
    ('x', float32),
    ('y', float32)
])
class Point: 
    # class to definy xy points
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 

@jit(nopython=True)
def line_test(pixel, r, theta, subsets):
    # function to test if a pixel lies above or below a line segment
    # Find start and end points of line segments for different cases of theta and r
    dx = subsets.shape[0]
    # Line segment to interogation point
    p1 = Point(0, 0)
    q1 = Point(pixel[0], pixel[1])
    
    # Vector magnitude cases
    #theta = np.degrees(theta)
    theta = theta % 360
    if r == 0:
        r = 1e-8
        
    # Rotation cases
    if theta == 0. or theta == 360.: # vertical to right
        x1 = r
        x2 = q1.x
        if x2 > x1:
            return False
        else:
            return True
    elif theta == 90.: # horizontal line above
        y1 = r
        y2 = q1.y
        if y2>y1:
            return False
        else:
            return True
    elif theta == 180.: # vertical to left
        x1 = -r
        x2 = q1.x
        if x2 > x1:
            return True
        else:
            return False
    elif theta == 270.: # horizontal below
        y1 = -r
        y2 = q1.y
        if y2 < y1:
            return False
        else:
            return True
    elif theta>0 and theta<180:
        theta = np.radians(theta)
        # Tangent line segment
        t1 = Point(r*np.cos(theta), r*np.sin(theta))
        m = -1*(np.cos(theta)/np.sin(theta))
        c = t1.y - m*t1.x
        y1 = q1.y
        y2 = m*q1.x + c
        if y1>y2:
            return False
        else:
            return True
    elif theta>180 and theta<360:
        theta = np.radians(theta)
        
        # Tangent line segment
        t1 = Point(r*np.cos(theta), r*np.sin(theta))
        m = -1*(np.cos(theta)/np.sin(theta))
        c = t1.y - m*t1.x
        
        y1 = q1.y
        y2 = m*q1.x + c
        if y1<y2:
            return False
        else:
            return True

@jit(uint8[:,:,:](uint8[:,:,:], float32, float32)) 
def subset_hsfilter(subsets, r, theta):
    # preallocate arrays
    hsfilter = np.zeros((subsets.shape[0],subsets.shape[0]))
    xc = hsfilter.shape[0]/2
    yc = hsfilter.shape[1]/2
    xc = int(xc)
    yc = int(yc)
    
    # Create x and y coordinates which are centred
    xs,ys = np.meshgrid(range(-xc, xc), range(-yc,yc))
    
    # iterate pixel by pixel
    for col in range(subsets.shape[0]):
        for row in range(subsets.shape[0]):
            #rasters through columns and rows for a given coordinate in xy
            x = xs[row,col]
            y = ys[row,col]
            # Note that y axis is mirrored
            pixel = [x, (-1*y)]
            # Test if pixel is beyond the discontinuity line
            flag  =line_test(pixel, r, theta, subsets)
            if flag:
                hsfilter[row,col] = True
            else:
                hsfilter[row,col] = False
                
    hs_subset = np.zeros(subsets.shape)
    hs_subset[:,:,0] = np.multiply(hsfilter,subsets[:,:,0])
    hs_subset[:,:,1] = np.multiply(hsfilter,subsets[:,:,1])
    return hs_subset



def hs_corr(x, subsets, d, prepared_ffts):
    # unpack vector
    r, theta = x
    # apply heaviside filter
    filtered_subsets = subset_hsfilter(subsets, r, theta)
    
    # assign to clean variables for legibility
    a = filtered_subsets[:,:,0]
    b = filtered_subsets[:,:,1]
    
    # cross correlate using freg - this is still faster than silly numpy correlate
    result = fxcorr(a, b, d, prepared_ffts)
    
    # normalise via division by active pixels
    act_px = np.count_nonzero(a != False)
    result = np.asarray(result)
    result[2] = result[2]/np.sqrt(act_px)
    
    # return vectorised result containign dx,dy,cc
    return result

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

def peak_angle(resolution, subsets, d, prepared_ffts):

    def hsin(x,subsets,d,prepared_ffts):
        _,_,c=hs_corr(x,subsets,d,prepared_ffts)
        return c
    theta = np.linspace(0, 360, resolution)
    r = np.full(len(theta), subsets.shape[0]/4)

    def A_mat(x):
        A = np.zeros((len(x),6))
        A[:,0]=np.squeeze(x**2)
        A[:,1]=np.squeeze(x)
        A[:,2]=np.ones(len(x))
        return A
    # generate trial Ax=b based on measured response
    atr = A_mat(theta)
    vals = np.c_[r,theta]

    # here we measure XCF
    btr = np.array([hsin(val, subsets, d, prepared_ffts) for val in vals])
    fx_measured = btr
    # x params

    x, _, _, _ = np.linalg.lstsq(atr, btr)

    # model
    theta_mod = np.linspace(0, 360, 1000)

    fx_mod = x[0]*theta_mod**2 + x[1] * theta_mod + x[2]

    fx_measured.max()

    loc = np.where(fx_measured.max() == fx_measured)
    theta_sol = theta[loc]

    return float(theta_sol[0])

def peak_r(theta, resolution, subsets, d, prepared_ffts):

    def hsin(x,subsets,d,prepared_ffts):
        _,_,c=hs_corr(x,subsets,d,prepared_ffts)
        return c
    # angle space
    theta = np.full((resolution),theta)
    

    # create search space in r
    r_max = subsets.shape[1]/3
    r = np.linspace(0,r_max,resolution)

    c = np.column_stack((r,theta))

    # search
    fx = np.array([hsin(val,subsets, d, prepared_ffts) for val in c])

    # create model

    def A_mat(x):
        A = np.zeros((len(x),6))
        A[:,0]=np.squeeze(x**2)
        A[:,1]=np.squeeze(x)
        A[:,2]=np.ones(len(x))
        return A

    # Ax = b

    A = A_mat(c[:,0])
    x, _, _, _ = np.linalg.lstsq(A,fx)
  
    # model

    r_mod = np.linspace(0,r_max,1000)
    fx_m = x[0]*r_mod**2 + x[1]*r_mod + x[2]

    loc = np.where(fx.max() == fx)
    r_sol = r_mod[loc]

    return float(r_sol[0]), np.mean(fx_m), np.max(fx_m) 

def least_squares_fit(res, subsets, d,prepared_ffts):
    theta = peak_angle(res, subsets, d, prepared_ffts)
    r, xcf_mean, xcf_peak = peak_r(theta, res, subsets, d, prepared_ffts)
    return r, theta, xcf_mean, xcf_peak

def minimise_rt_lstsq(subsets, d, prepared_ffts):

    def hsin(x,subsets,d,prepared_ffts):
        dx,dy,c=hs_corr(x,subsets,d,prepared_ffts)
        return dx, dy, c

    r, theta, xcf_mean, xcf_peak = least_squares_fit(15, subsets, d, prepared_ffts)

    dx, dy, cc = fxcorr(subsets[:,:,0], subsets[:,:,1], d, prepared_ffts)
    cc = cc/np.sqrt(subsets.shape[0]**2)
    if xcf_peak > cc:
        x = np.array([float(r), float(theta)])
        dx,dy,cc = hsin(x, subsets, d, prepared_ffts)
        return dx,dy,cc,r,theta,True
    else:
        return dx,dy,cc,False,False,False