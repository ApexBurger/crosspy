#%%
%load_ext autoreload
%autoreload 2

#%%
if __name__=='__main__':
    import os as o
    o.chdir('/Users/tom/Documents/GitHub/crosspy/Code_dev')
    #o.chdir('D:/DIC/crosspy/Code_dev')

    from pathlib import Path
    import time
    import numexpr

    import crosspy
    
    t0=time.time()

    folder_path=Path(r'/Users/tom/Documents/GitHub/crosspy/data/Tom')
    Images = Imset(folder_path,'tif',[0,1])

    # %% Instantiate and run the DIC

    # # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = 300, overlap_percentage = 75, xcf_mesh=400)

    # t1=time.time()
    # # build the dic class (but don't run it yet):
    dic_1stpass=DIC(Images,roi_1stpass,filter_settings)
    # # run the dic on specified images within the stack, and get displacements:
    dic_1stpass.run_sequential()
    dic_1stpass.plot_displacements()

    # # correct the images and instantiate a new DIC class
    corrected_images=dic_1stpass.correct(method='polynomial',printing=1)

    roi_2ndpass = dict(size_pass = 100, overlap_percentage = 80, xcf_mesh=500)
    dic_2ndpass = DIC(corrected_images,roi_2ndpass,filter_settings)

    # # run the second pass
    dic_2ndpass.run_sequential()
    dic_2ndpass.plot_displacements()
    # # strain calc
    # dic_2ndpass.strain_sequential()
    # print(time.time()-t0)

    #%%
    from StrainCalc import *
    dic_2ndpass.calculate_strain()
    dic_2ndpass.plot_strains()


#%%
dic_2ndpass.save_data()


# # %%
# from pathlib import Path
# import h5py
# file = ('/results'+str(dic_1stpass.roi[0])+'_'+str(dic_1stpass.roi[1]))

# with h5py.File(file, 'w') as f:
#     f.create_dataset('dx maps', data=dic_1stpass.dx_maps)
#     f.create_dataset('dy maps', data=dic_1stpass.dy_maps)
#     f.create_dataset('Peak heights', data=dic_1stpass.ph_maps)
#     f.create_dataset('Strain 11', data=dic_1stpass.strain_11)
#     f.create_dataset('Strain 22', data=dic_1stpass.strain_22)
#     f.create_dataset('Strain 12', data=dic_1stpass.strain_12)
#     f.create_dataset('Effective strain', data=dic_1stpass.strain_eff)

# %%
