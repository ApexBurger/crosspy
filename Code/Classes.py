# Classes used for the DPyC library package. 
# Written by Alexander Bergsmo and associated legends at Imperial College London 2020

from PIL import Image
import numpy as np

class Imset:

    # This instantiates a DIC image class which performs all necessary
    # operations on the dic image and gives it parameters
    # This class will allow images to be read and converted to DICable images
    # DIm classes have certain traits which will simplify DIC operations

    def __init__(self, folder_path,extension):

        self.folder = folder_path
        self.extension = extension
        self.paths = sorted(self.folder.glob('*.'+extension))
        self.names=[path.name for i,path in enumerate(self.paths)]

    def imload(self,numbers):
        # This loads the images provided by 'numbers'
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

    def imfilter(self, size_pass_1, overlap_pass_1, size_pass_2, overlap_pass_2):
        print('to-do')
        # this method sets up the filter windows and boundaries for ROIs
        #return fft_filter, h_filter

class im:
    # image sub class - will contain information of individual images
    pass

