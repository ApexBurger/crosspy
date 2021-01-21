# Bergsmo & McAuliffe 2020
import numpy as np
import crosspy

def gen_ROIs(imshape,roi):
    #generate where the subsets are for an arbitrary image size

    size_pass=roi[0]
    overlap=roi[1]/100

    rows=imshape[0]
    cols=imshape[1]

    spacing = round(size_pass*(1-overlap))
    
    #size_pass = subset size = 128,256 etc, overlap = proportion, eg 0.5, image1 and image2 must be same shape
    col_remainder=cols%spacing
    row_remainder=rows%spacing

    n_col_sets=int((cols-col_remainder)/spacing) #number of spacings in horizontal direction
    
    #ensure there is space for all the subsets
    while n_col_sets*spacing+size_pass>cols:
        n_col_sets=n_col_sets-1

    n_row_sets=int((rows-row_remainder)/spacing) #number of spacings in vertical direction

    while n_row_sets*spacing+size_pass>rows:
        n_row_sets=n_row_sets-1

    #generate an array of subset locations [top L corner row, top L corner col] - can get bot R row and col from known subset size
    ss_locations=np.zeros((n_row_sets*n_col_sets,2),dtype=np.int32)
    for c in range(0,n_col_sets):
        for r in range(0,n_row_sets):
            i = (c*n_row_sets)+r #index in this list
            ss_locations[i,0]=r*spacing
            ss_locations[i,1]=c*spacing
    
    return n_row_sets, n_col_sets, ss_locations, spacing

def get_subset(ims,size_pass,ss_locations,n,imno):
    #grab a subset from a stacked numpy array of images already loaded in RAM (imno is number in stack, n is subset number in ss_locations)
    
    subset=ims[ss_locations[n,0]:(ss_locations[n,0]+int(size_pass)),ss_locations[n,1]:(ss_locations[n,1]+int(size_pass)),imno]

    # if norm==True:
    #     #normalise
    #subset_norm=(subset-np.ones_like(subset)*np.mean(subset))/np.std(subset)
    # else:
    #     subset_norm=np.array(subset) 

    return subset #subset_norm


