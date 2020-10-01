#%%
# preamble
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy import DIC, Imset
import os
import cv2

# folder path

folder_path=Path(r'C:\Users\alexb\Dropbox\My PC (DESKTOP-R7R8K2S)\Documents\DIC\crosspy\data\tester_hor')
Images = Imset(folder_path,'tif')

# Image crop

crop_on = False

if crop_on:
    img = Images.imload([0,1])
    y = 175
    h = 1000
    x = 525
    w = 1000
    Images = img[y:y+h, x:x+w, :]
else:
    Images = Images.imload([0,1])

# Image crop

roi_1stpass = dict(size_pass = 40, overlap_percentage = 80, xcf_mesh=40)
filter_settings = [4,2,15,8]
# first pass

dic_1stpass = DIC(Images,roi_1stpass,filter_settings)

print(dic_1stpass.n_subsets)
dic_1stpass.run_sequential(cores=4,hs=False)

# dic_1stpass.plot_displacements()
dic_1stpass.calculate_strain()
dic_1stpass.plot_strains()
# %%
