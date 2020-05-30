# Image correction functions
import numpy as np
from scipy.optimize import least_squares
from scipy import interpolate
import cv2 as cv
import matplotlib.pyplot as plt

def im_correct(Images, d, method = "remapping" ):

    # Function to translate and rotate images to correct for 
    # rigid body rotation and translation

    shift_x=d.dx_maps
    shift_y=d.dy_maps
    pos_x=d.x_pos
    pos_y=d.y_pos

    # Load the first image and generate pixel arrays
    image_ref = Images.imload([0])
    x, y = np.meshgrid(range(0,np.size(image_ref,1)),range(0,np.size(image_ref,0)))
    image_c = np.zeros([np.size(image_ref,0),np.size(image_ref,1),len(Images)])
    image_c[:,:,0] = image_ref

    for i in range(0,len(Images)-1):
        # load second image and the shifts of a XCF pass
        im = Images.imload([i+1])
        x_shifts = shift_x[:, :, i]
        y_shifts = shift_y[:, :, i]

        # Correct rigid translation in x and y 
        # by subtracting the mean shift to every coordinate
        x_shifts_new = (x_shifts - np.mean(x_shifts)).flatten('F')
        y_shifts_new = (y_shifts - np.mean(y_shifts)).flatten('F')

        ## Guess rotation scheme
        # - rotation is confined between 90 degrees 
        # - with moving rotation axis
        params0 = np.array([np.size(im,1)/8, np.size(im,0)/1, -2.]) #initial [x_centre, y_centre, theta]
        ub = np.array([np.size(im,1), np.size(im,0), 45]) #upper bound [x_centre, y_centre, theta]
        lb = np.array([1, 1, -45]) #lower bound [x_centre, y_centre, theta]

        # A least squares regression scheme is utilised to approximate the 
        # centre of rotation and rotation angle
        l2_fit = rot_calc(x_shifts_new, y_shifts_new, pos_x, pos_y, params0, ub, lb)
        params1 = l2_fit.x
        
        # correct image rotation
        xc = params1[0] # Centre x
        yc = params1[1] # Centre y
        theta = params1[2] * np.pi/180 # Rotation theta
        rotation = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])     

        if method == 'remapping':
            # This method uses opencv remap - very fast and accurate
            x_cor, y_cor = x-np.mean(x_shifts), y-np.mean(y_shifts) #shifts
            points_cor = np.einsum('ji, mni -> jmn', rotation, np.dstack([x-xc, y-yc]))
            x_map = points_cor[0,:,:]+xc
            y_map = points_cor[1,:,:]+yc
            x_map = x_map.astype(np.float32)
            y_map = y_map.astype(np.float32)
            image_c[:,:,i+1] = cv2.remap(src=im, map1=x_map, map2=y_map, interpolation=cv2.INTER_CUBIC)
            
        elif method == 'affine':
            # Applies an affine translation and rotation - very fast but can cause blurry images
            rows,cols = im.shape
            # Translation
            x_shift = -np.mean(x_shifts)
            y_shift = -np.mean(y_shifts)
            M = np.float32([[1,0,x_shift],[0,1,y_shift]])
            dst = cv2.warpAffine(im,M,(cols,rows))
            # Rotation
            M = cv2.getRotationMatrix2D((xc,yc),np.degrees(theta),1)
            image_c[:,:,i+1] = cv2.warpAffine(dst,M,(cols,rows))

        # Below interpolates corrected values on grid following surface f(x,y) = z
        image_c[:,:,1] = interpolate.griddata(points=points_corr, values=im.flatten('F'), xi=(x,y), method='linear')
       
    return image_c

def rot_calc(x_shift, y_shift, x_pos, y_pos, params0, ub, lb):
    params1 = least_squares(rotation_fun, params0, jac='3-point', bounds=(lb, ub), args=(x_shift, y_shift, x_pos.flatten('F'), y_pos.flatten('F')))
    
    return params1

def rotation_fun(params, x_shift, y_shift, x_pos, y_pos):
    # WEEEEEEEEEEEEE
    xc = params[0]
    yc = params[1]
    theta = params[2] * np.pi/180

    # calculate rotation matrix
    rotm = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    
    # rotate existing grid about (xc, yc)
    temp = np.column_stack((x_pos-xc,y_pos-yc))
    new_points = np.matmul(temp,rotm)

    # calculate shifts
    shifts = np.column_stack((x_pos[:]-xc-new_points[:,0], y_pos[:]-yc-new_points[:,1]))

    # Find residuals
    a = shifts[:,0]-x_shift[:]
    b = shifts[:,0]-y_shift[:]
    resids = np.column_stack((a,b))
    resids = resids.flatten('F')

    return resids

def polynom_im_correct(d,printing,fn):
    #fits a least squares image to plane for a function to generate A (x,y) for Ax=b, with x parameters to fit.

    if fn==None:

        #create a function to generate the input parameters to the Ax=b mode
        def gen_A(x,y):
            A=np.zeros([len(x),6])
            A[:,0]=np.squeeze(x**2)
            A[:,1]=np.squeeze(x)
            A[:,2]=np.squeeze(y**2)
            A[:,3]=np.squeeze(y)
            A[:,4]=np.squeeze(x*y)
            A[:,5]=np.ones(len(x))
            return A

    else:
        gen_A=fn

    colmin=[]
    rowmin=[]
    colmax=[]
    rowmax=[]
    deformedimages=[]

    for imno in range(0,d.n_ims-1):

        #grab displacmenets from a dic pass
        dxs=np.squeeze(d.dx_maps[:,:,imno])
        dys=np.squeeze(d.dy_maps[:,:,imno])
        dims=dxs.shape

        #build grids of x, y
        xs=(np.linspace(0,dims[1]-1,dims[1])+0.5)*d.ss_spacing
        ys=(np.linspace(0,dims[0]-1,dims[0])+0.5)*d.ss_spacing
        xgrid,ygrid=np.meshgrid(xs,ys)

        #get overall co ordinate positions: x + dx, y + dy
        x_measured=xgrid+dxs
        y_measured=ygrid+dys
        #x_measured=dxs
        #y_measured=dys

        #unravel into 1D
        x_measured_1d=np.reshape(x_measured,(-1)) #rows
        y_measured_1d=np.reshape(y_measured,(-1)) #rows
        x=np.reshape(xgrid,(1,-1)).T #columns
        y=np.reshape(ygrid,(1,-1)).T #columns

        A=gen_A(x,y)

        #get least squares sol to (quadratic) plane of best fit
        x_params,_,_,_=np.linalg.lstsq(A,x_measured_1d)
        y_params,_,_,_=np.linalg.lstsq(A,y_measured_1d)

        #model values of the x, y planes
        xhat_1d=np.array(A@x_params)
        yhat_1d=np.array(A@y_params)

        #vectors taking original to LS deformed plane
        dxhat_1d=xhat_1d-np.squeeze(x)
        dyhat_1d=yhat_1d-np.squeeze(y)

        #now caculate the x,y positions for the full image with the parameters calculated above
        deformedimage=d.ims[:,:,imno+1].astype('float32')
        deformedimage=(deformedimage-np.mean(deformedimage))/np.std(deformedimage)

        dims=np.shape(deformedimage)
        xs=np.linspace(0,dims[1]-1,dims[1])
        ys=np.linspace(0,dims[0]-1,dims[0])
        xgrid,ygrid=np.meshgrid(xs,ys)

        x=np.reshape(xgrid,(1,-1)).T 
        y=np.reshape(ygrid,(1,-1)).T 

        A_full=gen_A(x,y)

        #model values of the x, y planes for full image
        xhat_1d=np.array(A_full@x_params)
        yhat_1d=np.array(A_full@y_params)
        xhat=np.reshape(xhat_1d,dims).astype('float32')
        yhat=np.reshape(yhat_1d,dims).astype('float32')

        #dxhat=xhat-xgrid.astype('float32')
        #dyhat=yhat-ygrid.astype('float32')

        deformedimage_corrected=cv.remap(deformedimage,xhat,yhat,cv.INTER_CUBIC,borderValue=0)

        #%% extract the centres of the corrected and undeformed images

        subset = deformedimage_corrected==0
        subset = ~subset #return true where we want to keep values

        #get the bottom left corner of the image
        inds=np.argwhere(subset)
        colmin+=[np.amin(inds[:,1])]
        rowmin+=[np.amin(inds[:,0])]
        colmax+=[np.amax(inds[:,1])]
        rowmax+=[np.amax(inds[:,0])]

        deformedimage_corrected[~subset]=0

        #remap to 0 - 255 dynamic range
        im_min=np.amin(deformedimage_corrected)
        im_max=np.amax(deformedimage_corrected)
        deformedimage_corrected-=im_min
        deformedimage_corrected=255*deformedimage_corrected/(im_max-im_min)

        #append the deformed image to a list
        deformedimages+=[deformedimage_corrected]

        if printing==1:
            deformedimage_corrected_toplot=deformedimage_corrected[rowmin[imno]:rowmax[imno],colmin[imno]:colmax[imno]]
            originalimage=d.ims[:,:,0].astype('float32')

            fig,([ax11,ax12,ax13],[ax21,ax22,ax23])=plt.subplots(nrows=2,ncols=3) 

            ax11.imshow(subset)
            ax11.set_title('remapping')
            ax11.axis('off')

            ax12.imshow(xhat)
            ax12.set_title('x-shifts')
            ax12.axis('off')

            ax13.imshow(yhat)
            ax13.set_title('y-shifts')
            ax13.axis('off')

            ax21.imshow(deformedimage,cmap='gray')
            ax21.set_title('deformed_original')
            ax21.axis('off')

            ax22.imshow(deformedimage_corrected_toplot,cmap='gray')
            ax22.set_title('deformed_corrected')
            ax22.axis('off')

            ax23.imshow(originalimage,cmap='gray')
            ax23.set_title('undeformed')
            ax23.axis('off')

            plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
            plt.show()

        new_rowmin=np.amin(rowmin)
        new_colmin=np.amin(colmin)
        new_rowmax=np.amax(rowmax)
        new_colmax=np.amax(colmax)

    #get the first image to begin the output stack (needs to be row x col x 1)
    imscorrected=np.expand_dims(d.ims[new_rowmin:new_rowmax,new_colmin:new_colmax,0].astype('float32'),-1)

    for i in range(0,d.n_ims-1):
        imscorrected=np.concatenate((imscorrected,np.expand_dims(deformedimages[i][new_rowmin:new_rowmax,new_colmin:new_colmax],-1)),axis=2)

    return imscorrected