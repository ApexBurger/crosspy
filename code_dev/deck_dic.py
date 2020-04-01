#%%
%load_ext autoreload
%autoreload 2

#%%
if __name__=='__main__':
    import os as o
    #o.chdir('/Users/tom/Documents/GitHub/crosspy/Code')
    o.chdir('D:/DIC/crosspy/Code_dev')
    from Classes import *
    from ImagePreparation import *
    from XCF import *
    from pathlib import Path
    import time
    import numexpr
    t0=time.time()

    folder_path = Path(r"D:\DIC\crosspy\data\Siyang")
    Images = Imset(folder_path,'tif')

    # %% Instantiate and run the DIC

    # # fft filter settings: high pass, high pass width, low pass, low pass width
    filter_settings=[4,2,15,8]
    roi_1stpass = dict(size_pass = 200, overlap_percentage = 70, xcf_mesh=250)

    # t1=time.time()
    # # build the dic class (but don't run it yet):
    dic_1stpass=DIC(Images,roi_1stpass,filter_settings)
    # # run the dic on specified images within the stack, and get displacements:
    dic_1stpass.run_sequential()
    dic_1stpass.plot_displacements()
    print(time.time()-t0)

    # # correct the images and instantiate a new DIC class
    # roi_2ndpass = dict(size_pass = 100, overlap_percentage = 80, xcf_mesh=250)
    # dic_2ndpass = DIC(dic_1stpass.correct(),roi_2ndpass,filter_settings)

    # # run the second pass
    # dic_2ndpass.run_sequential()
    # dic_2ndpass.plot_displacements()
    # # strain calc
    # dic_2ndpass.strain_sequential()
    # print(time.time()-t0)
    


#%%
from StrainCalc import *
dic_1stpass.calculate_strain()

dic_1stpass.plot_strains()
#%%

#dic_1stpass.plot_strain_meta()
import matplotlib.pyplot as plt
from plotting import *

plot_4(dic_1stpass,0,cmap="RdBu")



#%% save data

dic_1stpass.save_data()



# %%
from pathlib import Path
import h5py
file = ('/results'+str(dic_1stpass.roi[0])+'_'+str(dic_1stpass.roi[1]))

with h5py.File(file, 'w') as f:
    f.create_dataset('dx maps', data=dic_1stpass.dx_maps)
    f.create_dataset('dy maps', data=dic_1stpass.dy_maps)
    f.create_dataset('Peak heights', data=dic_1stpass.ph_maps)
    f.create_dataset('Strain 11', data=dic_1stpass.strain_11)
    f.create_dataset('Strain 22', data=dic_1stpass.strain_22)
    f.create_dataset('Strain 12', data=dic_1stpass.strain_12)
    f.create_dataset('Effective strain', data=dic_1stpass.strain_eff)

# %%
