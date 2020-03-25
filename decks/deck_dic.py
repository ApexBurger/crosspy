#%%
#%load_ext autoreload
#%autoreload 2

#%%
import matplotlib.pyplot as plt
from pathlib import Path
from crosspy import DIC, Imset

if __name__=='__main__':

    #folder_path = Path(r"C:\Users\tpm416\Documents\GitHub\crosspy\data\Siyang")
    folder_path=Path(r'/Users/tom/Documents/GitHub/crosspy/data/Siyang/')
    Images = Imset(folder_path,'tif')

    # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = 200, overlap_percentage = 70, xcf_mesh=250)

    # build the dic class (but don't run it yet):
    dic_1stpass=DIC(Images,roi_1stpass,filter_settings)
    # run the dic on specified images within the stack, and get displacements:
    dic_1stpass.run_sequential()
    dic_1stpass.plot_displacements()

    dic_1stpass.strain_sequential()
    dic_1stpass.plot_strains()

    # # correct the images and instantiate a new DIC class
    # roi_2ndpass = dict(size_pass = 200, overlap_percentage = 80, xcf_mesh=250)
    # dic_2ndpass = DIC(dic_1stpass.correct(),roi_2ndpass,filter_settings)

    # # run the second pass
    # dic_2ndpass.run_sequential()
    # dic_2ndpass.plot_displacements()
    
    # # strain calc
    # dic_2ndpass.strain_sequential(strain_method='l2')

# %%
