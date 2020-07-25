# Heaviside deck

load_ext autoreload
%autoreload 2


# %%

# %%

if __name__ == "__main__":
    import os as o

    from pathlib import Path
    import time

    import crosspy as xpy
    
    t0=time.time()

    folder_path=Path(r'D:\DIC\crosspy\data\Phantom')
    Images = xpy.Imset(folder_path,'tif',[0,1])

    ss_size_final = 70

    # # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = ss_size_final*2, overlap_percentage = 80, xcf_mesh=ss_size_final*2)

    # # build the dic class (but don't run it yet):
    dic_1stpass=xpy.DIC(Images[0,1],roi_1stpass,filter_settings)

    # # run the dic on specified images within the stack, and get displacements:
    dic_1stpass.run_sequential(cores=4)
    dic_1stpass.plot_displacements()

    # # correct the images and instantiate a new DIC class
    corrected_images=dic_1stpass.correct(method='polynomial',printing=1)

    roi_2ndpass = dict(size_pass = ss_size_final, overlap_percentage = 80, xcf_mesh=ss_size_final)
    dic_2ndpass = xpy.DIC(corrected_images,roi_2ndpass,filter_settings)

    # # run the second pass
    dic_2ndpass.run_sequential(cores=4, hs=True)
    dic_2ndpass.plot_displacements()

    dic_2ndpass.calculate_strain()
    dic_2ndpass.plot_strains()

    dic_2ndpass.save_data()


# %%

print(dic_2ndpass.hs_maps[:,:,0])

import matplotlib.pyplot as plt 

plt.hist(dic_2ndpass.th_maps[:,:,0])

# %%
