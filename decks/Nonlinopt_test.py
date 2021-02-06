""" This is a test deck written by Alexander Bergsmo 2021.
The purpose of this deck is to attempt an implementation of the classical DIC formulation relying on first order shape
functions as well as extension by heaviside DIC (or subset splitting methodology)

The general process of the method is outlined below

1. A ROI is selected within the map
2. The ROI is discretised into subsets similar to coarse-fine method
3. The undeformed subset is applied a binary mask the size of the ROI and the deformed image remains the same
4. Instead of generating a correlogram per subset, the correlation coefficient is minimized following a Hessian approximation
5. Once a suitable match is made, the subset is split according to parameters r,theta in polar coordinates from subset centre
6. The external side of the split subset is then minimized according to the Hessian approximation.
7. If there exists a "jump" - aka if the split subset moves - then a discontinious displacement has been detected and treatment is accepted



"""
#%%

def ROI_selector(a, b):
    """ This function selects a rectangular region of interest within the loaded images (a,b)
    """
    pass




def interp_subset(array,factor):
    """ This function up-samples an array with spline biquintic interpolation
    This is done "lazily" by scipy
    """
    x = np.arange(0,array.shape[1],1)
    y = np.arange(0,array.shape[0],1)
    z = array

    # Interpolation is biquintic - can be changed to lesser orders
    f = interpolate.interp2d(x,y,z,kind="quintic")

    # Upsampling is px_old * factor = px_new
    step = 1/factor

    xnew = np.arange(0,array.shape[1],step)
    ynew = np.arange(0,array.shape[0],step)
    # normalise
    image = f(ynew,xnew)
    image *= (255.0/image.max())

    # floor to remove below zero values
    image = np.where(image < 0, 0, image)
    return  image.astype(np.float32)

def interp_subset_cv(array,factor):
    dimx = int(array.shape[1]*factor)
    dimy = int(array.shape[0]*factor)
    image = cv.resize(array.astype(np.float32), ( dimx, dimy ), interpolation = cv.INTER_CUBIC )
    return image

def formmask(a, b, n, d):
    """ generate a mask
    """
    mask = 0
    return mask

def maskunion(a,b,mask):
    pass

def correlation_criterion(x, a, b):
    """ This function returns the correlation coefficient for respective input parameters
    x = [u, v, dudx, dudy, dvdy, dvdx]
    """
    cc = 0
    return cc

def corr_min(a,b):
    """ This function returns the minima/extremum (depending on criteria) of the correlation criterion for subsets f,g
    """

def subset_split():
    """ can take from existing split
    """
    pass

def run_opt(d,imnos=[0,1],hs=False, cc_t=0.):
    """ Runtime function for optimising correlation around first order shape function
    """

    # preallocate for this DIC pair
    phs=np.zeros(d.n_subsets)
    u = np.zeros(d.n_subsets)
    v = np.zeros(d.n_subsets)
    dudx = np.zeros(d.n_subsets)
    dvdy = np.zeros(d.n_subsets)
    dudy = np.zeros(d.n_subsets)
    dvdx = np.zeros(d.n_subsets)

    for i in range(d.n_subsets):
        print(i)

    




class DIC_classical():
    """ incorporate the main runtime into the class later...
    """
    pass





#%%
import numpy as np 
import matplotlib.pyplot as plt 
import crosspy
from pathlib import Path
from crosspy import DIC, Imset
import os
import cv2 as cv

# folder path

folder_path=Path(r'C:\Users\alexb\Desktop\crosspy\data\Disc_tester')
# import images
Images = Imset(folder_path,'tif',[0,1])

#%%
filter_settings=[4,2,15,8]
images = Images.imload([0,1])

roi = dict(size_pass = 32, overlap_percentage = 90, xcf_mesh=200)
dic_1stpass = DIC(Images,roi,filter_settings)

dic_1stpass.roi
#%%

im1 = images[:,:,0]
im2 = images[:,:,1]
d = dic_1stpass
#%%
# get some subsets

n = 0

from crosspy.ImagePreparation import get_subset

ref=get_subset(d.ims,d.roi[0],d.ss_locations,n,0)
test=get_subset(d.ims,d.roi[0],d.ss_locations,n,1)

#%%
# attempt to upsample
from scipy import interpolate
import cv2 as cv

z = ref

f = 2
%timeit a = interp_subset(ref,factor=f)
%timeit b = interp_subset_cv(test,factor=f)
fig,ax = plt.subplots(1,2)

ax[0].imshow(a)
ax[1].imshow(b)
# %%

# form with mask

def formmask(roigrid, f):
    # mask is a regularized grid upsampled
    xs = roigrid.shape[1]*f
    ys = roigrid.shape[0]*f

    # create a mask with same size as roigrid
    mask = np.zeros(shape=(ys,xs))

    return mask

def maskunion(subset,x1,y1,f,mask):
    # pad transformed subest to achieve same size as deformed map
    # pad mask with subset width
    # mask = np.pad(mask,subset.shape[0])
    off = subset.shape[0]
    # indices
    w = int(subset.shape[0])
    print(x1,w,off)
    x1 = (x1*f) #+ off
    x2 = x1 + w
    y1 = (y1*f) #+ off
    y2 = y1 + w
    print(x1,x2)
    mask[y1:y2,x1:x2] = subset
    union = mask
    return union

roigrid = im1
mask = formmask(roigrid,f)
masku = maskunion(a,d.ss_locations[n,0],d.ss_locations[n,1],f,mask)


fig,ax = plt.subplots(1,2)

ax[0].imshow(masku)
ax[1].imshow(interp_subset(im2,f))

# %%
import cv2 as cv
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

ax = masku.astype(np.float32)
bx = interp_subset(im2,f).astype(np.float32)

print(ax.shape,bx.shape)
cc = correlogram(ax, bx, 'cv.TM_CCOEFF_NORMED')


fig, axs = plt.subplots(1,3,figsize=(20,30))
axs[0].imshow(ax)
axs[1].imshow(bx)
axs[2].imshow(cc)
# transform subset using first order transformation

# %%
def reg(a,b,f, method = 'cv.TM_SQDIFF_NORMED'):
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

    # if abs(col_shift) > a.shape[0]/2 or abs(row_shift) > a.shape[1]/2:
    #     col_shift = 0
    #     row_shift = 0
    
    return -col_shift/f, -row_shift/f, cc, res


f = 4
ref_up = interp_subset_cv(ref,f).astype(np.float32)
test_up = interp_subset_cv(test,f).astype(np.float32)

dxf,dyf,ccf,res2 = reg(ref_up,test_up,f, method='cv.TM_SQDIFF_NORMED')

dx1,dy1,cc1,res1 = reg(ref,test,1, method='cv.TM_SQDIFF_NORMED')

fix,ax = plt.subplots(1,2,figsize=(10,20))

ax[0].imshow(res2)
ax[0].set_title("dx = {}, dy = {}, cc={}".format(dxf,dyf,ccf))
ax[1].imshow(res1)
ax[1].set_title("dx = {}, dy = {}, cc={}".format(dx1,dy1,cc1))
# %%

def fxcorr_cv(subset1,subset2,d):
    """ This function uses openCV to upsample and correlate subsets
    It is fast but the subpixel accuracy in inefficient as upsampling
    is done in real space over the entire subset
    
    To improve, look at what ncorr does...
    or, write a custom function which finds the interger shift location, then performs
    a "mini" cross correlation on a smaller subset, upsampled
    """
    # upsampling factor
    f = d.roi[2]/100

    # upsample subsets by bicubic interpolation
    ref_up = interp_subset_cv(ref,f).astype(np.float32)
    test_up = interp_subset_cv(test,f).astype(np.float32)


    col_shift,row_shift,ccmax,_ = reg(ref_up,test_up,f, method='cv.TM_CCOEFF_NORMED')

    return col_shift, row_shift, ccmax

fxcorr_cv(ref,test,d)

# def warp_subset(x, subset):
#     u, v, dudx, dudy, dvdx, dvdy = x

#     c = int(subset.shape[0]/2)
#     w = int(subset.shape[0]/2)
#     xs = np.arange(-w,w,1)
#     ys = xs
#     xx, yy = np.meshgrid(xs,ys)
    
    
#     xcur = xref + u
#     warped_subset = xcur
#     return warped_subset


# c = int(ref.shape[0]/2)
# w = int(ref.shape[0]/2)
# xs = np.arange(-w,w,1)
# ys = xs
# xx, yy = np.meshgrid(xs,ys)

# xcur = lambda xref,xc,yref,yc,u,dudx,dudy:\
#         xref+u+dudx(xref-xc)+dudy(yref-yc)
# ycur = lambda xref,xc,yref,yc,v,dvdy,dvdx:\
#         yref+v+dvdx(xref-xc)+dvdy(yref-yc)

# x_f = xcur(xx,0, )

# plt.imshow(yy)
# %%
