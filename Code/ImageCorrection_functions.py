# Image correction functions
import numpy as np
from scipy.optimize import least_squares
from scipy import interpolate

def im_correct(Images, shift_x, shift_y, pos_x, pos_y):
    # Function to translate and rotate images to correct for 
    # rigid body rotation and translation

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
        # Apply corrections by creating the an array containing the
        # corrected coordinates and then mapping those coordinates
        # to existing pixel grid
        im_shifts = np.column_stack((x.flatten('F')-np.mean(x_shifts), y.flatten('F')-np.mean(y_shifts)))
        im_rot = np.dot(np.column_stack((im_shifts[:,0]-xc,im_shifts[:,1]-yc)),rotation)
        im_correct = np.column_stack((im_rot[:,[0]]+xc, im_rot[:,[1]]+yc))
        points_corr = [(row[0],row[1]) for row in im_correct]

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