#%%
%load_ext autoreload
%autoreload 2

#%%
if __name__=='__main__':
    import os as o
    from pathlib import Path
    import time
    import numexpr
    import crosspy as xpy
    
    t0=time.time()

    folder_path=Path(r'../data/Tom') 
    Images = xpy.Imset(folder_path,'tif',[0,1])

    # # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = 300, overlap_percentage = 75, xcf_mesh=400)

    # # build the dic class (but don't run it yet):
    dic_1stpass=xpy.DIC(Images[0,1],roi_1stpass,filter_settings)

    # # run the dic on specified images within the stack, and get displacements:
    dic_1stpass.run_sequential()
    dic_1stpass.plot_displacements()
    
    # # correct the images and instantiate a new DIC class
    corrected_images=dic_1stpass.correct(method='polynomial',printing=1)

    roi_2ndpass = dict(size_pass = 250, overlap_percentage = 80, xcf_mesh=500)
    dic_2ndpass = xpy.DIC(corrected_images,roi_2ndpass,filter_settings,savingfolder=dic_1stpass.folder)

    # # run the second pass
    dic_2ndpass.run_sequential()
    dic_2ndpass.plot_displacements()

    dic_2ndpass.calculate_strain()
    dic_2ndpass.plot_strains()

    dic_2ndpass.save_data()


# %%
