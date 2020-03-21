# Bergsmo & McAuliffe 2020

from PIL import Image
import numpy as np
from imprep_functions import *
from XCF_functions import *
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

def subset_compare(d,imnos,subset_n):
    #grab the reference and test subsets, and get subpixel registration
    ref=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[0])
    test=get_subset(d.ims,d.roi[0],d.ss_locations,subset_n,imnos[1])
    #get the displacements
    dxs,dys,phs=fxcorr(ref,test,d)
    return dxs,dys,phs

def subset_compare_par(d,imnos,chunks,cores):
    subset_compare_partial=functools.partial(subset_compare,d,imnos)

    with Pool(cores) as p:
        output = p.map(subset_compare_partial,range(0,len(d.ss_locations)),chunks)
    return output

class dic:
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

    def run(self,imnos=[0,1],par=False,chunks=10,cores=None):
        
        #preallocate for this DIC pair
        phs=np.zeros(self.n_subsets)
        dxs=np.zeros(self.n_subsets)
        dys=np.zeros(self.n_subsets)

        #serial version of this function
        if par==False:

            for subset_n in range(0,self.n_subsets):
                dxs[subset_n],dys[subset_n],phs[subset_n]=subset_compare(self,imnos,subset_n)

        #parallel version
        if par==True:
            output = subset_compare_par(self,imnos,chunks,cores)

            for subset_n in range(0,len(output)):
                dxs[subset_n]=output[subset_n][0]
                dys[subset_n]=output[subset_n][1]
                phs[subset_n]=output[subset_n][2]

        #translate best_dxs etc back onto image grid
        dx_map=np.reshape(dxs,(self.n_rows,self.n_cols),'F')
        dy_map=np.reshape(dys,(self.n_rows,self.n_cols),'F')
        ph_map=np.reshape(phs,(self.n_rows,self.n_cols),'F')

        return dx_map,dy_map,ph_map

    def run_sequential(self,par=False,chunks=10,cores=None):
        #Perform DIC on consecutive images, using the previous as a reference.
        #chunks and cores only apply if par=True ; if cores=None looks for maximum for your system.

        #preallocate for all DIC pairs
        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        for i in range(0,self.n_ims-1):
            if par: suffix=' (parallel) '
            print('Running sequential DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=self.run([i,i+1],par,chunks,cores)

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

        for i in range(0,self.n_ims-1):
            if par: suffix=' (parallel) '
            print('Running cumulative DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
            dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=self.run([0,i+1],par,chunks,cores)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps
        
        #return dx_maps, dy_maps, ph_maps

    def plot_results(self,colmap='RdBu'):

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

    def strain_calc(self):
        # strain calc based on deformation map
        pass

class Im(Imset):

    def __init__(self):
        pass
