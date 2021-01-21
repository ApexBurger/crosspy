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





def subset_gen(a, b, n, d):
    """ this function generates upsampled subset images and applies a mask to them
    """
    xc, yc = d.ss_locations[n,:]
    w = d.roi[0]

    mask = np.zeros(a.shape)
    mask = np.pad(mask, a.shape[0])
    l = int(w/2)
    xc = xc+a.shape[1]
    yc = yc+a.shape[0]
    mask[yc-l:yc+l,xc-l:xc+l] = np.ones((w,w))

    a_masked = np.pad(a, a.shape[0]) * mask
    
    return a_masked

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
import cv2

# folder path

folder_path=Path(r'C:\Users\alexb\Desktop\crosspy\data\Disc_tester')
# import images
Images = Imset(folder_path,'tif',[0,1])

#%%
filter_settings=[4,2,15,8]
images = Images.imload([0,1])

roi = dict(size_pass = 32, overlap_percentage = 90, xcf_mesh=150)
dic_1stpass = DIC(Images,roi,filter_settings)

dic_1stpass.roi
#%%

a = images[:,:,0]
b = images[:,:,1]
d = dic_1stpass
n = 500


a_masked = subset_gen(a,b,n,d)

fig,ax = plt.subplots(1,2,figsize=(15,10))
ax[0].imshow(a_masked)
ax[1].imshow(b)
#%%
# artifically deform an image

def image_deform(imarray, du,dv):
    """ This function artifically deforms an image array
    according to input parameters


    """


    return imarray_sheared

im_sheared = image_shear(images[:,:,0],-1,-1)

plt.imshow(im_sheared)

#%% Initial guess



#%% Newton - Raphson method for finding deformation maping parameters

def NR_def():
    pass


