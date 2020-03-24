# Bergsmo & McAuliffe 2020

from PIL import Image
import numpy as np
from imprep_functions import *
from ImageCorrection_functions import *
from XCF_functions import *
from runDIC_functions import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools
from StrainCalc_functions import *

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

    def __init__(self,images,roi,filter_settings):

        #if fed an Imset class
        if isinstance(images,Imset):
            if images.n_ims<1:
                raise Exception('No images found!')

            self.imageset=images
            self.ims=images.imload(range(0,images.n_ims))
            self.n_ims=images.n_ims

        if isinstance(images,np.ndarray):
            self.ims=images
            self.n_ims=images.shape[2]

        self.roi=list([roi['size_pass'],roi['overlap_percentage'],roi['xcf_mesh']])

        self.n_rows,self.n_cols,self.ss_locations,self.ss_spacing=gen_ROIs(self.ims.shape[0:2],self.roi)
        self.n_subsets=self.ss_locations.shape[0]

        self.fftfilter,self.hfilter=gen_filters(self.roi,filter_settings)
        self.filter_settings=filter_settings

        self.x_pos = self.ss_locations[:,0].reshape(self.n_rows,self.n_cols)+roi['size_pass']/2
        self.y_pos = self.ss_locations[:,1].reshape(self.n_rows,self.n_cols)+roi['size_pass']/2

    def run_sequential(self,par=False,cores=None,chunk_length=50):
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
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=run_DIC(self,[i,i+1],par,cores,chunk_length=50)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps

        #return dx_maps, dy_maps, ph_maps

    def run_cumulative(self,par=False,cores=None,chunk_length=50):
        #Perform DIC on sequential images, using the first as a reference.
        #chunks and cores only apply if par=True ; if cores=None looks for maximum for your system.

        #preallocate for all DIC pairs
        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        suffix=''

        for i in range(0,self.n_ims):
            if par: suffix=' (parallel) '

            print('Running cumulative DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=run_DIC(self,[0,i+1],par,cores,chunk_length)

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

    def strain_sequential(self, strain_method='l2'):
        if self.dx_maps.any()==False:
            raise Exception('No displacements available for strain calculation!')
        #Perform strain calculation on consecutive images, using the previous as a reference.
        self.mapnos = np.size(self.dx_maps, 2)

        # preallocate arrays - strain, rotation and deformation gradient
        # are stored in tensors for each subset for each map
        strain = np.zeros((self.n_rows, self.n_cols, 3, 3, self.mapnos))
        strain_eff = np.zeros((self.n_rows,self.n_cols, 1, self.mapnos))
        rotation = np.zeros((self.n_rows, self.n_cols, 3, 3, self.mapnos))
        F = np.zeros((self.n_rows, self.n_cols, 3, 3, self.mapnos))
        suffix=''

        for i in range(0,self.mapnos):
            print('Calculating sequential strain on map ' +str(i+1)+' of '+str(self.mapnos)+suffix)
            strain[:,:,:,:,i], strain_eff[:,:,:,i], rotation[:,:,:,:,i], F = strain_calc(self,
                mapnos=i, strain_method=strain_method)

        self.strain_11 = strain[:,:,0,0,:]
        self.strain_22 = strain[:,:,1,1,:]
        self.strain_12 = strain[:,:,1,0,:]
        self.strain_eff = strain_eff
        self.rotation = rotation
        self.deformation_gradient = F

    def strain_cumulative(self):
        #Perform strain calculation on sequential images, using the first as a reference.
        pass
    def correct(self):
        images_corrected=im_correct(self.imageset,self)
        return images_corrected
    
    def plot_strains(self,colmap='RdBu'):
        
        if self.strain_eff.any()==False:
            raise Exception('No strain results to plot!')

        for i in range(0,self.mapnos):
            fig,((ax11,ax12),(ax21,ax22))=plt.subplots(2,2,figsize=(10,10)) 
            ax11.imshow(self.strain_11[:,:,i],cmap=colmap)
            ax11.set_title('XX strains, map '+str(i+1))
            ax12.imshow(self.strain_22[:,:,i],cmap=colmap)
            ax12.set_title('YY strains, map '+str(i+1))
            ax21.imshow(self.strain_12[:,:,i],cmap=colmap)
            ax21.set_title('Shear strains, map '+str(i+1))
            ax22.imshow(self.strain_eff[:,:,i],cmap=colmap)
            ax22.set_title('Effective strain, map '+str(i+1))
            plt.show()

class Im(Imset):

    def __init__(self):
        pass
