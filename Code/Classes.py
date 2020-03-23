# Bergsmo & McAuliffe 2020

from PIL import Image
import numpy as np
from imprep_functions import *
from XCF_functions import *
import matplotlib.pyplot as plt

class Imset:

    # This instantiates a DIC image class, holding metadata. The image is not loaded until the method .imload() is called.
    # Inputs are folder path and image file extension - will identify all images of that same extension within the folder.

    def __init__(self, folder_path,extension):

        self.folder = folder_path
        self.extension = extension
        self.paths = sorted(self.folder.glob('*.'+extension))
        self.names=[path.name for i,path in enumerate(self.paths)]
        self.n_ims=len(self.names)

        print('Found '+str(self.n_ims)+' images')

    def __len__(self):
        return self.n_ims
        # this is redundant but nice to define

    def imload(self,numbers):
        # This loads the images identifed in the Imset enumerated by 'numbers'

        for i, n in enumerate(numbers):

            path=self.paths[n]
            #open as grayscale
            image=Image.open(str(path)).convert('LA')
            imarray=np.array(image)[:,:,0]

            imarray=np.squeeze(imarray)

            #loads a 2D array then adds on a new axis
            imarray=np.expand_dims(imarray,-1)
            #concatenate all the arrays
            if i==0:
                imarray_stack=np.array(imarray)
            else:
                imarray_stack=np.concatenate((imarray_stack,imarray),axis=2)

        if len(numbers)==1:
            imarray_stack=np.squeeze(imarray_stack)
        
        return imarray_stack

class DIC:
# A class that holds information relating to a DIC run on an imageset
# Call .run(filter_settings...) to map x, y displacements
# At the moment it runs on the top left portion of the image (ie. if there are pixels to the right and down that...
# can't be accommodated by square subsets of given roi_size, they will be ignored).

    def __init__(self,imageset,roi,filter_settings):

        if imageset.n_ims<1:
            raise Exception('No images found!')

        self.ims=imageset.imload(range(0,imageset.n_ims))
        self.n_ims=imageset.n_ims
        self.n_rows,self.n_cols,self.ss_locations,self.ss_spacing=gen_ROIs(self.ims.shape[0:2],roi)
        self.n_subsets=self.ss_locations.shape[0]
        self.roi=roi
        self.filter_settings=filter_settings
        self.x_pos = self.ss_locations[:,0].reshape(self.n_rows,self.n_cols)+roi['size_pass']/2
        self.y_pos = self.ss_locations[:,1].reshape(self.n_rows,self.n_cols)+roi['size_pass']/2

    def run(self,imnos=[0,1]):
            
        phs=np.zeros(self.n_subsets)
        dxs=np.zeros(self.n_subsets)
        dys=np.zeros(self.n_subsets)

        for subset_n in range(0,self.n_subsets):
            #grab the reference and test subsets, and get subpixel registration
            ref=get_subset(self.ims,self.roi,self.ss_locations,subset_n,imnos[0])
            test=get_subset(self.ims,self.roi,self.ss_locations,subset_n,imnos[1])
            dxs[subset_n],dys[subset_n],phs[subset_n]=fxcorr(ref,test,self.roi,self.filter_settings)

            #translate best_dxs etc back onto image grid
            dx_map=np.reshape(dxs,(self.n_rows,self.n_cols),'F')
            dy_map=np.reshape(dys,(self.n_rows,self.n_cols),'F')
            ph_map=np.reshape(phs,(self.n_rows,self.n_cols),'F')

        return dx_map,dy_map,ph_map

    def run_sequential(self):
        #Perform DIC on consecutive images, using the previous as a reference.

        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        for i in range(0,self.n_ims-1):
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=self.run([i,i+1])

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps

        #return dx_maps, dy_maps, ph_maps

    def run_cumulative(self):
        #Perform DIC on sequential images, using the first as a reference.

        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        for i in range(0,self.n_ims-1):
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=self.run([0,i+1])

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps
        
        #return dx_maps, dy_maps, ph_maps

    def plot_results(self,colmap='plasma'):

        if self.ph_maps.any()==False:
            raise Exception('No DIC results to plot!')

        for i in range(0,self.n_ims-1):
            fig,(ax11,ax12,ax13)=plt.subplots(1,3,figsize=(10,10)) 
            ax11.imshow(self.dx_maps[:,:,i],cmap=colmap)
            ax11.set_title('X-displacements, map '+str(i+1))
            ax12.imshow(self.dy_maps[:,:,i],cmap=colmap)
            ax12.set_title('Y-displacements, map '+str(i+1))
            ax13.imshow(self.ph_maps[:,:,i],cmap=colmap)
            ax13.set_title('CC peak heights, map '+str(i+1))


    def strain_calc(self):
        # strain calc based on deformation map
        pass

class Im(Imset):

    def __init__(self):
        pass
