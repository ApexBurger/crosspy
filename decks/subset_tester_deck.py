# %% imports

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy import DIC, Imset
import os
import cv2
from numba import double, jit, int32, float32, uint8, jitclass
import time

# %% Parameters

size = 150
overlap = 80

# %% Import test data

# folder path

folder_path=Path(r'C:\Users\alexb\Dropbox\My PC (DESKTOP-R7R8K2S)\Documents\DIC\crosspy\data\Phantom')
Images = Imset(folder_path,'tif')
Ims = Images.imload([0,1])

fig = plt.figure()
ax = fig.gca()
ax.imshow(Ims[:,:,1])
#ax.set_xticks(np.arange(0, Images.shape[0], Images.shape[0]/10))
#ax.set_yticks(np.arange(0, Images.shape[1], Images.shape[1]/10))
plt.grid()


# %% Test run
if __name__=='__main__':
    roi_1stpass = dict(size_pass = size, overlap_percentage = overlap, xcf_mesh=size)
    filter_settings = [4,2,15,8]
    # first pass

    dic_1stpass = DIC(Images,roi_1stpass,filter_settings)
    print(dic_1stpass.n_subsets)
    dic_1stpass.run_sequential(cores=4,hs=False)
    dic_1stpass.plot_displacements()
    dic_1stpass.calculate_strain()
    dic_1stpass.plot_strains()

# %%

from numba import njit, prange

@njit(parallel=True)
def prange_test(A):
    s = 0
    # Without "parallel=True" in the jit-decorator
    # the prange statement is equivalent to range
    for i in prange(A.shape[0]):
        s += A[i]
    return s

if __name__ == '__main__':
    %timeit prange_test(np.arange(1000000000))
    %timeit prange_test(np.arange(1000000000))


# %% get a specific subset

subsetID = 0
d = dic_disc
subset_a=get_subset(d.ims,d.roi[0],d.ss_locations,subsetID,[0])
subset_b=get_subset(d.ims,d.roi[0],d.ss_locations,subsetID,[1])
subsets = np.zeros((size,size,2))
subsets[:,:,0] = subset_a[:,:,0]
subsets[:,:,1] = subset_b[:,:,0]

# %% Masking function

@jit(nopython=True)
def line_test(pixel, r, theta, subsets):
    # function to test if a pixel lies above or below a line segment
    # Find start and end points of line segments for different cases of theta and r
    dx = subsets.shape[0]
    # Line segment to interogation point
    p1x = 0.
    p1y = 0.

    q1x = pixel[0]
    q1y = pixel[1]
    # Vector magnitude cases
    #theta = np.degrees(theta)
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
        theta = np.radians(theta)
        # Tangent line segment
        t1x = np.array([r*np.cos(theta)])
        t1y = np.array([r*np.sin(theta)])
        m = -1*(np.cos(theta)/np.sin(theta))
        c = t1y - m*t1x
        y1 = q1y
        y2 = m*q1x + c
        if y1>y2:
            return False
        else:
            return True
    elif theta>180 and theta<360:
        theta = np.radians(theta)
        # Tangent line segment
        t1x = np.array([r*np.cos(theta)])
        t1y = np.array([r*np.sin(theta)])
        m = -1*(np.cos(theta)/np.sin(theta))
        c = t1y - m*t1x
        y1 = q1y
        y2 = m*q1x + c
        if y1<y2:
            return False
        else:
            return True

@jit(nopython=False) 
def gen_hsfilter(subsets, r, theta):
    # preallocate arrays
    filter_arr = np.zeros((subsets.shape[0],subsets.shape[0]), dtype="int")
    xc = filter_arr.shape[0]/2
    yc = filter_arr.shape[1]/2
    xc = xc
    yc = yc

    # Create x and y coordinates which are centred
    x_length = np.linspace(-xc, xc,subsets.shape[0])
    y_length = np.linspace(-yc,yc,subsets.shape[0])
    xs,ys = np.meshgrid(x_length,y_length)
    hsfilter = filt_loop(subsets, r, theta, xs, ys, filter_arr)
    return hsfilter


@jit(nopython=True) 
def filt_loop(subsets, r, theta, xs, ys, filter_arr):
    cols = np.arange(subsets.shape[0])
    rows = np.arange(subsets.shape[0])
    # iterate pixel by pixel
    for col in cols:
        for row in rows:
            #rasters through columns and rows for a given coordinate in xy
            # Note that y axis is mirrored
            pixel = np.array([xs[row,col], (-1*ys[row,col])])
            # Test if pixel is beyond the discontinuity line
            filter_arr[row,col]  = line_test(pixel, r, theta, subsets)
            # if flag:
            #     filter_arr[row,col] = True
            # else:
            #     filter_arr[row,col] = False
                
    return filter_arr

# %% Time it

start = time.time()
gen_hsfilter(subsets,r=18, theta=10)
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))

# NOW THE FUNCTION IS COMPILED, RE-TIME IT EXECUTING FROM CACHE
start = time.time()
gen_hsfilter(subsets,r=18, theta=10)
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))

# %% Normalised cross correlation function

@jit(nopython=False)
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
    unout = out
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

@jit(nopython=False)
def hs_corr(x, subsets, dic, prepared_ffts):
    # unpack vector
    r, theta = x
    # generate heaviside filter
    hsfilter = gen_hsfilter(subsets, r, theta)
    
    # assign to clean variables for legibility
    a = subsets[:,:,0]#*hsfilter
    b = subsets[:,:,1]*hsfilter
    
    # cross correlate using freg - this is still faster than silly numpy correlate
    results = XCF.fxcorr(a, b, dic, prepared_ffts)
    # invert selection
    hsfilter_inv = ~hsfilter

    # apply to clean variables
    c = subsets[:,:,0]# * hsfilter_inv
    d = subsets[:,:,1] * hsfilter_inv

    results_inv = XCF.fxcorr(c, d, dic, prepared_ffts)

@jit(nopython=False)
def hs_corr_norm(x,subsets,d,prepared_ffts):
    # unpack vector
    r,theta = x

    #apply heaviside
    hsfilter = gen_hsfilter(subsets, r, theta)

    # apply to clean variables
    a = subsets[:,:,0]# * hsfilter
    b = np.multiply(subsets[:,:,1] , hsfilter)

    # get inverse
    hsfilter_inv = ~hsfilter

    # apply to clean variables
    c = subsets[:,:,0]# * hsfilter_inv
    d = subsets[:,:,1] * hsfilter_inv

    # cross correlate using freg
    results = np.array(ncc(b,a))
    results_inv = np.array(ncc(d,c))
    
    # fig, ax = plt.subplots(1,2)
    # ax[0].imshow(a)
    # ax[1].imshow(d)
    # normalise by active pixels
    results = np.concatenate((results,results_inv),axis=None)
    return results

# %% discontinuity estimate

@jit(nopython=False)
def A_mat(x):
    A = np.zeros((len(x),6))
    A[:,0]=np.squeeze(x**4)
    A[:,1]=np.squeeze(x**3)
    A[:,2]=np.squeeze(x**2)
    A[:,3]=np.squeeze(x)
    A[:,4]=np.ones(len(x))
    return A

def hsin(x,subsets,d,prepared_ffts):
    _,_,c,_,_,_=hs_corr_norm(x,subsets,d,prepared_ffts)
    return c

# @jit(nopython=False)
def locate_disc(resolution, subsets, d, prepared_ffts):
    
    # first find theta
    theta = np.linspace(0, 360, resolution)
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
    th_f = np.radians(resolution)*2
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
    theta_arr = np.full((resolution),theta)

    # create search space in r
    r_max = subsets.shape[1]/3
    r = np.linspace(0,r_max,resolution)

    c = np.column_stack((r,theta_arr))

    # search
    fx = np.array([hsin(val,subsets, d, prepared_ffts) for val in c])

    # measure
    loc = np.where(fx.max() == fx)
    r_sol = r[loc]

    return int(r_sol[0]), theta, rsq

def minimise_rt_lstsq(subsets, d, prepared_ffts):

    # obtain peak r and theta

    r, theta, rsq = locate_disc(20, subsets, d, prepared_ffts)
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
        dxu,dyu,ccu = XCF.fxcorr(subsets[:,:,0],subsets[:,:,1], d, prepared_ffts)
        return np.array([dxu,dyu,ccu,False,False,False,False,False,False])
    elif (u != 0 or v != 0) and rsq > 0.07: # dx or dy are not the same - presence of discontinuity
        # One value does not match direction of displacment - indicating a kinematic shift
        dx,dy,cc,dxi,dyi,cci = hs_corr(x, subsets, d, prepared_ffts)
        u = dx - dxi
        v = dy - dyi
        j = np.sqrt(u**2+v**2)
        return np.array([dx,dy,cc,r,theta,True,u,v,j])
    else:
        # fit is bad indicating no symmetry
        dxu,dyu,ccu = XCF.fxcorr(subsets[:,:,0],subsets[:,:,1], d, prepared_ffts)
        return np.array([dxu,dyu,ccu,False,False,False,False,False,False])

# %% Run dic for heaviside

roi_disc = dict(size_pass = 80, overlap_percentage = 60, xcf_mesh=80)
dic_disc = DIC(ims, roi_disc,filter_settings)

@jit(nopython=False)
def subset_compare(subset_n,d,imnos,prepared_ffts,discontinuity=False):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements 
    if discontinuity == True:
        subsets = np.stack((ref,test),axis=2)
        results = minimise_rt_lstsq(subsets,d,prepared_ffts)
        return results
    else:
        dxs,dys,phs=fxcorr(ref,test,d,prepared_ffts)
        return dxs,dys,phs,None, None, None

prepared_ffts = plan_ffts(dic_disc,ffttype='fftw_numpy')
@jit(nopython=False, parallel=True)
def run_DIC(d,prepared_ffts,imnos=[0,1],discontinuity=False, cores=None, ffttype='fftw_numpy'):
    #fft type can be : fftw_numpy (default), fftw_scipy, else defaults to numpy

    
    #discontinuity enables or disables slip trace tracking via heaviside filtering


    #pre-allocate for this DIC pair
    phs=np.zeros(d.n_subsets)
    dxs=np.zeros(d.n_subsets)
    dys=np.zeros(d.n_subsets)

    #enable the pyfftw cache for speedup
    pyfftw.interfaces.cache.enable()   

    #check for discontinuity tracker
    if discontinuity == True:
        #pre-allocate heavisde assets
        r = np.zeros(d.n_subsets)
        theta = np.zeros(d.n_subsets)
        hson = np.zeros(d.n_subsets)
        u = np.zeros(d.n_subsets)
        v = np.zeros(d.n_subsets)
        j = np.zeros(d.n_subsets)
        results = np.zeros((len(r),9))

        #main DIC loop - iterates through subsets
        i = 0
        for subset_n in prange(0,d.n_subsets):
            results[i,:] = subset_compare(subset_n,d,imnos,prepared_ffts, True)
            i += 1
        dxs = results[:,0]
        dys = results[:,1]
        phs = results[:,2]
        r = results[:,3]
        theta = results[:,4]
        hson = results[:,5]
        u = results[:,6]
        v = results[:,7]
        j = results[:,8]
        #format to maps
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')
        r_map = np.reshape(r,(d.n_rows,d.n_cols),'F')
        theta_map = np.reshape(theta,(d.n_rows,d.n_cols),'F')
        hson_map = np.reshape(hson,(d.n_rows,d.n_cols),'F')
        u_map = np.reshape(u,(d.n_rows,d.n_cols),'F')
        v_map = np.reshape(v,(d.n_rows,d.n_cols),'F')
        j_map = np.reshape(j,(d.n_rows,d.n_cols),'F')

        return dx_map, dy_map, ph_map, r_map, theta_map, hson_map, u_map, v_map, j_map
    else:
        for subset_n in range(0,d.n_subsets):
            dxs[subset_n],dys[subset_n],phs[subset_n]=subset_compare(d,imnos,subset_n,prepared_ffts)

        #translate best_dxs etc back onto image grid
        dx_map=np.reshape(dxs,(d.n_rows,d.n_cols),'F')
        dy_map=np.reshape(dys,(d.n_rows,d.n_cols),'F')
        ph_map=np.reshape(phs,(d.n_rows,d.n_cols),'F')

        return dx_map,dy_map,ph_map

if __name__=='__main__':
    start = time.time()
    dx_map, dy_map, ph_map, r_map, theta_map, hson_map, u_map, v_map, j_map = run_DIC_cuda(dic_disc, discontinuity=True, prepared_ffts=prepared_ffts)
    end = time.time()
    print(end - start)