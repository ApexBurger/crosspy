# Bergsmo & McAuliffe 2020

from PIL import Image
import numpy as np
from imprep_functions import *
from XCF_functions import *
from runDIC_functions import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools

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
        self.roi=list([roi['size_pass'],roi['overlap_percentage'],roi['xcf_mesh']])

        self.n_rows,self.n_cols,self.ss_locations,self.ss_spacing=gen_ROIs(self.ims.shape[0:2],self.roi)
        self.n_subsets=self.ss_locations.shape[0]

        self.fftfilter,self.hfilter=gen_filters(self.roi,filter_settings)
        self.filter_settings=filter_settings

    def run_sequential(self,par=False,chunks=10,cores=None):
        #Perform DIC on consecutive images, using the previous as a reference.
        #chunks and cores only apply if par=True ; if cores=None looks for maximum for your system.

        #preallocate for all DIC pairs
        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        suffix=''

        for i in range(0,self.n_ims-1):
            if par: suffix=' (parallel) '
            print('Running sequential DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=run_DIC(self,[i,i+1],par,chunks,cores)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps

        #return dx_maps, dy_maps, ph_maps

    def run_cumulative(self,par=False,chunks=10,cores=None):
        #Perform DIC on sequential images, using the first as a reference.
        #chunks and cores only apply if par=True ; if cores=None looks for maximum for your system.

        #preallocate for all DIC pairs
        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        suffix=''

        for i in range(0,self.n_ims-1):
            if par: suffix=' (parallel) '

            print('Running cumulative DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=run_DIC(self,[0,i+1],par,chunks,cores)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps
        
        #return dx_maps, dy_maps, ph_maps

    def plot_displacements(self,colmap='RdBu'):

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
            
            plt.show()

    def strain_sequential(self,par=False,chunks=10,cores=None):
        #Perform strain calculation on consecutive images, using the previous as a reference.

        #preallocate for all DIC pairs
        strain_11 = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        strain_22 = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        strain_12 = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        rotation = 
        strain_eff = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        suffix=''

        for i in range(0,self.n_ims-1):
            if par: suffix=' (parallel) '
            print('Running sequential DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
            strain_11[:,:,i],strain_22[:,:,i],strain_12[:,:,i], strain_eff[:,:,i]=run_DIC(self)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps

        #return dx_maps, dy_maps, ph_maps

    def strain_cumulative(self):
        #Perform strain calculation on sequential images, using the first as a reference.



class Im(Imset):

    def __init__(self):
        pass
