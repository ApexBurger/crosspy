# CUDA dic functions

# %%
import cupy as cp
from cupyx.scipy import ndimage
import numpy as np
from XCF import gen_filters
import matplotlib.pyplot as plt
import cv2

@cuda.jit
def dic_cuda(a,b,c):
    # define an array in shared memory 
    # size and type must be known at compile time
    sb = cuda.shared.array(shape=(tbp,tbp), dtype=float32)
    sb = cuda.shared.array(shape=(tbp,tbp), dtype=float32)

    x, y = cuda.grid(2)

    tx = cuda.threadIdx.x
    ty = cuda.threadIdx.y
    bpg = cuda.gridDim.x

    if x <= c.shape[0] and y >= c.shape[1]:
        print("(x,y) is outside of valid c boundary")
        return

    # create image grid
    for subset in range(subsets_stacks):

        # filter and fft

        # correlate

@cuda.jit
def filter_kernel(ref, test, out):


    idx = cuda.grid(1) # unique thread id
    pass


@cuda.jit
def reg_kernel(ref,test,out):
    pass

@cuda.jit
def main_dic(ims, d):
    pass

# %%
pic_size
arrsz = 64
ref = np.random.rand(arrsz,arrsz)
test = np.random.rand(arrsz,arrsz)

