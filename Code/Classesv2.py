# Classes used for the DPyC library package. 
# Written by Alexander Bergsmo at Imperial College London 2020
from PIL import Image
import numpy as np
import os

class Imset:
    # This instantiates a DIC image class which performs all necessary
    # operations on the dic image and gives it parameters
    # This class will allow images to be read and converted to DICable images
    # DIm classes have certain traits which will simplify DIC operations
    def __init__(self, filename, filepath, imageorder):
        self.filename = filename # this is confusing
        self.path = filepath
        self.ftype = filename.split('.')[-1]
        self.files = [i for i in os.dir(self.path)]


    def imload(self):
        # This loads the images
        im_list = []
         #needs to get the files
            

        if self.ftype == 'tiff' or self.ftype == 'tif':
            imarray = Image.

        return imarray

    def conv_greyscale(self):
        # This converts a loaded image to greyscale
        imarray =
        imarray_conv = 

        return imarray_conv

    def imfilter(self, size_pass_1, overlap_pass_1, size_pass_2, overlap_pass_2):
        # this method sets up the filter windows and boundaries for ROIs

        return fft_filter, h_filter

        

class im(imset):
    # image sub class - will contain information of individual images
    pass