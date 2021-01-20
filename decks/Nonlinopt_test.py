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





def subset_gen(a, b, xc, yc, w):
    """ this function generates upsampled subset images and applies a mask to them
    """
    mask = np.zeros(a.shape)

    l = int(w/2)
    mask[yc-l:yc+l,xc-l:xc+l] = np.ones((w,w))
    
    a_masked = a * mask

    return a_masked

def correlation_criterion(x, f, g):
    """ This function returns the correlation coefficient for respective input parameters
    x = [u, v, dudx, dudy, dvdy, dvdx]
    """
    cc = 0
    return cc

def corr_min(f,g):
    """ This function returns the minima/extremum (depending on criteria) of the correlation criterion for subsets f,g
    """

def subset_split():
    """ can take from existing split
    """
    pass

class DIC_classical(DIC):

    




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

images = Images.imload([0,1])
#%%

a = images[:,:,0]
b = images[:,:,1]

a_masked = subset_gen(a,b,32,32,32)

fig,ax = plt.subplots(2)
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


