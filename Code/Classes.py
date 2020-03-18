# Classes used for the DPyC library package. 

# Written by Alexander Bergsmo at Imperial College London 2020

from PIL import Image
import numpy as np

class Imset:

    # This instantiates a DIC image class which performs all necessary
    # operations on the dic image and gives it parameters
    # This class will allow images to be read and converted to DICable images
    # DIm classes have certain traits which will simplify DIC operations

    def __init__(self, folder_path,extension):

        self.path = folder_path
        self.extension = extension
        self.imagepaths = sorted(self.path.glob('*.'+extension))

    def imload(self,numbers):
        # This loads the images provided by 'numbers'
        for i, n in enumerate(numbers):

            path=self.imagepaths[n]
            #open as grayscale
            image=Image.open(str(path))
            imarray=np.array(image)
            imarray=np.squeeze(imarray)

            #loads a 2D array then adds on a new axis
            imarray=np.expand_dims(imarray,-1)
                        #concatenate all the arrays
            if i==0:
                imarray_stack=np.array(imarray)
            else:
                imarray_stack=np.concatenate((imarray_stack,imarray),axis=2)
        return imarray_stack


    def conv_greyscale(self):

        # This converts a loaded image to greyscale

        imarray =

        imarray_conv = 



        return imarray_conv

    # Filter settings are currently only definable here:
    filter_settings = [4,2,16,32]
    # TO DO - GUI to help user chose the filter

    def gen_filters(self, size_pass, overlap_pass):
        # Inputs:
        #   roisize = subset size - (128, 256 etc.)
        #   fpasset = [high pass cut off, high pass width, low pass cut off, low pass width]
        # outputs:
        #   fftfilter = non idea (gaussian) band pass filter
        #   hfilter = hanning filter

        pi = np.pi
        cos = np.cos
        dot = np.dot
        sqrt = np.sqrt
        exp = np.exp

        lcutoff = fpassset[2]
        lwidth = fpassset[3]/2
        hcutoff = fpassset[0]
        hwidth = fpassset[1]/2

        if lcutoff < hcutoff:
            print('low pass filter smaller than high pass filter')
            return

        # generate square grid
        u = range(0, roisize)
        # meshgrid function
        meshv, meshu = np.meshgrid(u,u)
        meshvf = meshv-roisize/2-0.5
        meshuf = meshu-roisize/2-0.5

        # create Hann window
        hfilter = (cos(((pi*(meshuf)/roisize)))*(cos((pi*meshvf/roisize))))

        # create fft filter
        distf = sqrt((meshvf*meshvf)+(meshuf*meshuf))

        # lowpass
        lfftfilter = exp(-((distf-lcutoff)/sqrt(2)*lwidth/2)**2)
        lfftfilter[distf > (lcutoff+2*lwidth)] = 0
        lfftfilter[distf < (lcutoff+2*lwidth)] = 1


        # highpass
        hfftfilter = exp(-((hcutoff-distf)/sqrt(2)*hwidth/2))**2
        hfftfilter[distf < (hcutoff-2*hwidth)] = 0
        hfftfilter[distf > hcutoff] = 1

        # combine
        fftfilter = hfftfilter*lfftfilter
        fftfilter = np.fft.fftshift(fftfilter)
        
        return fftfilter, hfilter



        








