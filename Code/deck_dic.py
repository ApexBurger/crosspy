#%%
%load_ext autoreload
%autoreload 2

#%%

import os as o
o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
#o.chdir('D:/DIC/crosspy/Code')
from Classes import *
from imprep_functions import *
from XCF_functions import *
from pathlib import Path
import time

t0=time.time()

folder_path = Path(r"/Users/tom/Documents/GitHub/crosspy/Data/Siyang")
Images = Imset(folder_path,'tif')

fig = plt.figure()
plt.imshow(Images.imload([1]))

# %% Instantiate and run the DIC

#fft filter settings: high pass, high pass width, low pass, low pass width
filter_settings=[4,2,15,8]
roi_1stpass = dict(size_pass = 200, overlap_percentage = 70, xcf_mesh=250)

#build the dic class (but don't run it yet):
dic_1stpass=DIC(Images,roi_1stpass,filter_settings)
#run the dic on specified images within the stack, and get displacements:
dic_1stpass.run_sequential(par=True)
dic_1stpass.plot_displacements()

print(time.time()-t0)


#%%
Images_cor =im_correct(Images, dic_1stpass.dx_maps, dic_1stpass.dy_maps, dic_1stpass.x_pos, dic_1stpass.y_pos)