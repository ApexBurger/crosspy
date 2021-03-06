# Bergsmo & McAuliffe 2020

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import functools
import time
import h5py
import numpy as np
import os as o
import crosspy


class Imset:
    """
    This instantiates a DIC image class, holding metadata. 
    The image is not loaded until the method .imload() is called.

    Inputs:
        folder_path (pathlib.Path object) - path to the folder of images.
        extension (str) - image file extension to search for.
        indices (list or int) - specify which images to take.
    
    Methods:
        self.imload(numbers) - Load and return corresponding images identiifed with metadata.
    
    Attributes:
        self.folder - file path to saving folder.
        self.extension - extension of images.
        self.paths - file paths of images.
        self.names - file names of images.
        self.n_ims - number of images in the set.

    """
    def __init__(self, folder_path,extension,indices=None):
        self.folder = folder_path
        self.extension = extension
        foundpaths = sorted(self.folder.glob('*.'+extension))

        if indices==None:
            self.paths = foundpaths
        
        else:
            if type(indices)==int:
                indices=[indices]
            self.paths=[foundpaths[pathno] for _,pathno in enumerate(indices)]

        self.names=[path.name for i,path in enumerate(self.paths)]
        self.n_ims=len(self.names)

        #print('Found '+str(self.n_ims)+' images')


    def __len__(self):
        return self.n_ims
        # this is redundant but nice to define


    def __getitem__(self,nslice):
        Imset2=Imset(self.folder,self.extension,nslice)
        return Imset2


    def imload(self,numbers):
        """
        Loads the images identified in the parent Imset which are enumerated by the input.

        Inputs
            numbers (list of int): Which images in the imset to load.

        Outputs:
            imarray_stack: numpy array of images of dimensions [H, W, n_images]

        """
        # for massive images
        Image.MAX_IMAGE_PIXELS = None
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
    """
    A class that holds information relating to a DIC run on an imageset
    At the moment it runs on the top left portion of the image (ie. if there are pixels to the right and down 
    that can't be accommodated by square subsets of given roi_size, they will be ignored).

    Inputs:
        images (Imset) - Imset containing metadata for loading image stacks.
        roi (dict) - ROI settings for DIC. 
                     This should be dictionary with keys "size_pass", "overlap_percentage", "XCF_mesh".
        filter_settings (list) - High / low pass Fourier filter settings. 
                     This should be [high_pass_filter, high_pass_filter_width, low_pass_filter, low_pass_filter_width] in px.
        savingfolder (str or pathlib.Path) - location to write results to if requested.

    Methods:
        self.run_sequential(cores, ffttype, hs, cc_t) - run DIC on images with displacements/strains relative to previous image in stack. Saves results to attributes.
        self.run_cumulative(cores, ffttype, hs, cc_t) - run DIC on images with displacements/strains relative to first image in stack. Saves results to attributes.
        
        ## These can only be run after .run_sequential() or .run_cumulative() ##
        self.calculate_strain(strain_method) - infer the strain using calculated displacements. Saves results to attributes.
        self.correct(method, printing, fn) - correct the images by removing the calculated displacements. Returns an updated stack of images.
        self.plot_displacements(colmap) - plot the calculated displacements if available in attributes.

        ## These can only be run after .calculate_strains() ##
        self.plot_strains(colmap) - plot the calculated strains if available in attributes.
        self.plot_strain_meta(bins) - plots histograms of strains if available in attributes.
        self.save_data(output_folder) - save results to hdf5.
    
    Attributes:
        ## Available after .run_sequential() or .run_cumulative() ##
        self.dx_maps - x displacement maps
        self.dy_maps - y displacement maps
        self.ph_maps - cross-correlation peak height maps

        ## and if heaviside being used ##
        self.rd_maps
        self.th_maps
        self.hs_maps
        self.js_maps

        ## Available after .calculate_strains() ##
        self.strain_11 - 11 (xx) strain
        self.strain_22 - 22 (yy) strain
        self.strain_12 - 12 (xy) strain
        self.strain_eff - von Mises effective strain
        self.rotation - 12 (xy) rotation
        self.deformation_gradient - full deformation gradient tensor

    """

    # Call .run(filter_settings...) to map x, y displacements


    def __init__(self,images,roi,filter_settings,savingfolder=None):
        #if fed an Imset class
        if isinstance(images,(crosspy.Imset,Imset)):
            if images.n_ims<1:
                raise Exception('No images found!')

            self.imageset=images
            self.ims=images.imload(range(0,images.n_ims))
            self.n_ims=images.n_ims
            self.folder = images.folder

        else:
            self.ims=images
            self.n_ims=images.shape[2]
            self.folder = savingfolder
        
        self.roi=list([roi['size_pass'],roi['overlap_percentage'],roi['xcf_mesh']])
        self.n_rows,self.n_cols,self.ss_locations,self.ss_spacing=crosspy.gen_ROIs(self.ims.shape[0:2],self.roi)
        self.n_subsets=self.ss_locations.shape[0]

        self.fftfilter,self.hfilter = crosspy.gen_filters(self.roi,filter_settings)
        self.filter_settings=filter_settings

        self.x_pos = self.ss_locations[:,0].reshape(self.n_rows,self.n_cols)+roi['size_pass']/2
        self.y_pos = self.ss_locations[:,1].reshape(self.n_rows,self.n_cols)+roi['size_pass']/2

    def run_sequential(self,cores=1, ffttype='fftw_numpy', hs=False, cc_t=0., px_size=None, cormeth="efficient"):
        """
        Run DIC with reference images as previous in the stack. Generates attributes for this class.

        Inputs
            cores (int) - how many CPU cores to run subsets in parallel on.
            ffttype (str) - one of 'fftw_numpy', 'fftw_scipy', or 'numpy' - which FFT implementation to use.
            hs (bool) - run Heaviside or not.
            cc_t (float) - cross-correlation peak height threshold.
            px_size (float) = pixel size in length used to plot j maps
            cormeth (str) = method of correlation in disc cont, either "efficient" or "cv" -> cv allows for NSQ and subpixel


        Outputs
            self.rd_maps
            self.th_maps
            self.hs_maps
            self.js_maps
            self.ph_maps
            self.dx_maps
            self.dy_maps

        """

        #preallocate for all DIC pairs
        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        suffix=' ...'
        t0=time.time()

        if hs == True:
            rd_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
            th_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
            hs_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
            js_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

            for i in range(0,self.n_ims-1):
                print('Running sequential DIC on image pair '+str(i+1)+' of '+str(self.n_ims-1)+suffix+', total subsets per image: '+str(self.n_subsets)) 
                dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i],rd_maps[:,:,i],th_maps[:,:,i],hs_maps[:,:,i],js_maps[:,:,i] \
                = crosspy.run_DIC(d=self, imnos=[i,i+1],cores=cores, hs=True, cc_t=cc_t,cormeth=cormeth)

            self.rd_maps = rd_maps
            self.th_maps = th_maps
            self.hs_maps = hs_maps

            if px_size == None:
                print("No defined pixel size, J-map is returned in pixels")
                self.js_maps = js_maps
            else:
                self.js_maps = js_maps * px_size

        else:
            for i in range(0,self.n_ims-1):
                print('Running sequential DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
                dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i] = crosspy.run_DIC(self, imnos=[i,i+1], cores=cores, cormeth=cormeth)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps
        print('... Completed in (s) '+str(time.time()-t0))

        #return dx_maps, dy_maps, ph_maps

    def run_cumulative(self,cores=1,ffttype='fftw_numpy', hs=False, cc_t=0., px_size=None,cormeth="efficient"):
        """
        Run DIC with first image as reference. Generates attributes for this class.

        Inputs
            cores (int) - how many CPU cores to run subsets in parallel on.
            ffttype (str) - one of 'fftw_numpy', 'fftw_scipy', or 'numpy' - which FFT implementation to use.
            hs (bool) - run Heaviside or not.
            cc_t (float) - cross-correlation peak height threshold.
            px_size = pixel size in length used to plot j maps

        Outputs
            self.rd_maps
            self.th_maps
            self.hs_maps
            self.js_maps
            self.ph_maps
            self.dx_maps
            self.dy_maps
 
        """

        #preallocate for all DIC pairs
        ph_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dx_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
        dy_maps=np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

        suffix=' ...'
        t0=time.time()

        if hs == True:
            rd_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
            th_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
            hs_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))
            js_maps = np.zeros((self.n_rows,self.n_cols,self.n_ims-1))

            for i in range(0,self.n_ims-1):
                print('Running sequential DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix +', total subsets per image: ' + str(self.n_subsets))
                dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i],rd_maps[:,:,i],th_maps[:,:,i],hs_maps[:,:,i],js_maps[:,:,i] \
                = crosspy.run_DIC(d=self, imnos=[0,i+1],cores=cores, hs=True, cc_t=cc_t,cormeth=cormeth)

            self.rd_maps = rd_maps
            self.th_maps = th_maps
            self.hs_maps = hs_maps
            if px_size == None:
                print("No defined pixel size, J-map is returned in pixels")
                self.js_maps = js_maps
            else:
                self.js_maps = js_maps * px_size

        else:
            for i in range(0,self.n_ims-1):
                print('Running cumulative DIC on image pair ' +str(i+1)+' of '+str(self.n_ims-1)+suffix)
                dx_maps[:,:,i],dy_maps[:,:,i],ph_maps[:,:,i]=crosspy.run_DIC(self,[0,i+1], cores=cores, cormeth=cormeth)

        self.ph_maps=ph_maps
        self.dx_maps=dx_maps
        self.dy_maps=dy_maps
        print('... Completed in (s) '+str(time.time()-t0))
        #return dx_maps, dy_maps, ph_maps

    def plot_displacements(self,colmap='RdBu'):

        if self.ph_maps.any()==False:
            raise Exception('No DIC results to plot!')
        if hasattr(self, "rd_maps") == False:
            for i in range(0,self.n_ims-1):
                fig,(ax11,ax12,ax13)=plt.subplots(1,3,figsize=(10,10)) 
                ax11.imshow(self.dx_maps[:,:,i],cmap=colmap)
                ax11.set_title('X-displacements, map '+str(i+1))
                ax12.imshow(self.dy_maps[:,:,i],cmap=colmap)
                ax12.set_title('Y-displacements, map '+str(i+1))
                ax13.imshow(self.ph_maps[:,:,i],cmap=colmap)
                ax13.set_title('CC peak heights, map '+str(i+1))
                
                plt.show()
        else:
            for i in range(0,self.n_ims-1):
                fig,((ax11,ax12,ax13),(ax21,ax22,ax23))=plt.subplots(2,3,figsize=(10,10)) 
                ax11.imshow(self.dx_maps[:,:,i],cmap=colmap)
                ax11.set_title('X-displacements, map '+str(i+1))
                ax12.imshow(self.dy_maps[:,:,i],cmap=colmap)
                ax12.set_title('Y-displacements, map '+str(i+1))
                ax13.imshow(self.ph_maps[:,:,i],cmap=colmap)
                ax13.set_title('CC peak heights, map '+str(i+1))
                ax21.imshow(self.rd_maps[:,:,i],cmap=colmap)
                ax21.set_title('R , map '+str(i+1))
                ax22.imshow(self.th_maps[:,:,i],cmap=colmap)
                ax22.set_title('Theta , map '+str(i+1))
                ax23.imshow(self.hs_maps[:,:,i],cmap=colmap)
                ax23.set_title('Heavisided, map '+str(i+1))
                
                plt.show()

    def calculate_strain(self, strain_method='l2'):
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
        suffix=' ...'
        t0=time.time()

        for i in range(0,self.mapnos):
            print('Calculating strain on map ' +str(i+1)+' of '+str(self.mapnos)+suffix)
            strain[:,:,:,:,i], strain_eff[:,:,:,i], rotation[:,:,:,:,i], F[:,:,:,:,i] = crosspy.strain_calc(self,
                mapnos=i, strain_method=strain_method)

        self.strain_11 = strain[:,:,0,0,:]
        self.strain_22 = strain[:,:,1,1,:]
        self.strain_12 = strain[:,:,1,0,:]
        self.strain_eff = strain_eff[:,:,0,:]
        self.rotation = rotation[:,:,:,:,:]
        self.deformation_gradient = F[:,:,:,:,:]
        print('... Completed in (s) '+str(time.time()-t0))

    def correct(self,method='polynomial',printing=0,fn=None):
        """
        Remove displacements saved to DIC attributes from the input Imset by fitting and subtracting a polynomial. Returns corrected images.

        Inputs:
            method (str) - one of 'polynomial' or 'rigid'  displacement correction, latter is affine.
            printing (bool) - whether to print to console.
            fn - optional functional form to solve for BG correction. If none, quadratic assumed. See polynom_im_correct documentation for more.

        Outputs:
            images_corrected - numpy array of corrected images. 

        """

        print('Correcting images based on DIC results ...')
        t0=time.time()

        #choose one of affine or polynomial methods
        if method=='rigid':
            images_corrected=crosspy.im_correct(self,printing)
        elif method=='polynomial':
            images_corrected=crosspy.polynom_im_correct(self,printing,fn)
        else:
            raise Exception('Method not recognised!')
            

        print('... Completed in (s) '+str(time.time()-t0))
        return images_corrected
    
    def plot_strains(self,colmap='RdBu',vmin=0, vmax=0.2):
        print('Quick plotting strains')
        if self.strain_eff.any()==False:
            raise Exception('No strain results to plot!')

        for i in range(0,self.mapnos):
            crosspy.plot_4(self, i, colmap)


    def plot_strain_meta(self, bins=5):
        for i in range(0,self.mapnos):
            fig,((ax11,ax12),(ax21,ax22))=plt.subplots(2,2,figsize=(10,10)) 
            ax11.hist(self.strain_11[:,:,i].flatten(), density=True, bins=bins)
            ax11.set_title('XX strains histogram, map '+str(i+1))
            ax12.hist(self.strain_22[:,:,i].flatten(), density=True, bins=bins)
            ax12.set_title('YY strains histogram, map '+str(i+1))
            ax21.hist(self.strain_12[:,:,i].flatten(), density=True, bins=bins)
            ax21.set_title('Shear strains histogram, map '+str(i+1))
            ax22.hist(self.strain_eff[:,:,i].flatten(), density=True, bins=bins)
            ax22.set_title('Effective strain histogram, map '+str(i+1))

    def save_data(self, output_folder = None):
        currentfolder=o.getcwd()
        o.chdir(self.folder)

        if output_folder == None:
            output_folder = self.folder
        file = ('results'+str(self.roi[0])+'_'+str(self.roi[1]))
        if self.dx_maps.any() == False:
            raise Exception('No displacement data to save!')
        elif self.strain_11.any() == False:
            print('No strain maps available, only displacement data will be saved!')
        else:
            print('Saving displacement and strain data')
        
        # e11 = pd.DataFrame(data=self.strain_11)

        with h5py.File(file, 'w') as f:
            f.create_dataset('dx maps', data=self.dx_maps)
            f.create_dataset('dy maps', data=self.dy_maps)
            f.create_dataset('Peak heights', data=self.ph_maps)
            f.create_dataset('Strain 11', data=self.strain_11)
            f.create_dataset('Strain 22', data=self.strain_22)
            f.create_dataset('Strain 12', data=self.strain_12)
            f.create_dataset('Effective strain', data=self.strain_eff)

        o.chdir(currentfolder)

