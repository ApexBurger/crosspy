## THIS FUNCTION IS UNUSED - THE ACTIVE VERSION LIES IN hs.py


from numba import double, jit, njit, vectorize
from numba import int32, float32, uint8, float64, int64, boolean
import numpy as np
import time

# Apply a line and step function
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