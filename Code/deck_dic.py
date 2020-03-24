#%%
# %load_ext autoreload
# %autoreload 2

#%%
if __name__=='__main__':
    import os as o
    # o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
    o.chdir(r'C:\Users\tpm416\Documents\GitHub\crosspy\data\Siyang')
    from Classes import *
    from imprep_functions import *
    from XCF_functions import *
    from pathlib import Path
    import time
    from ImageCorrection_functions import *

    t0=time.time()

    folder_path = Path(r"C:\Users\tpm416\Documents\GitHub\crosspy\data\Siyang")
    Images = Imset(folder_path,'tif')

    fig = plt.figure()
    plt.imshow(Images.imload([1]))

    # %% Instantiate and run the DIC

    # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = 200, overlap_percentage = 70, xcf_mesh=250)

    t1=time.time()
    # build the dic class (but don't run it yet):
    dic_1stpass=DIC(Images,roi_1stpass,filter_settings)
    # run the dic on specified images within the stack, and get displacements:
    dic_1stpass.run_sequential()
    # dic_1stpass.plot_displacements()
    print(time.time()-t1)

    # correct the images and instantiate a new DIC class
    roi_2ndpass = dict(size_pass = 100, overlap_percentage = 80, xcf_mesh=250)
    dic_2ndpass = DIC(dic_1stpass.correct(),roi_2ndpass,filter_settings)

    # run the second pass
    dic_2ndpass.run_sequential()
    dic_2ndpass.plot_displacements()
    # strain calc

    dic_2ndpass.strain_sequential(strain_method='l2')
    print(time.time()-t0)


# %%
