# Classes
from PIL import Image
import numpy as np
class DIm:
    # This instantiates a DIC image class which performs all necessary
    # operations on the dic image and gives it parameters
    # This class will allow images to be read and converted to DICable images
    # DIm classes have certain traits which will simplify DIC operations
    def __init__(self, filename, filepath):
        self.array = np.array(Image.open(filepath))
        self.imnum = [i for i in filename]
        self.imtype = filename.split('.')[-1]
    # 
    def convert_grayscale():
        pass

    

    
