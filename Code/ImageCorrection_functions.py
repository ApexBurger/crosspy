# Image correction functions
import numpy as np
from scipy.optimize import least_squares

def im_correct(Images, shift_x, shift_y, pos_x, pos_y):
    # Function to translate and rotate images to correct for 
    # rigid body rotation and translation

    # Load the first image
    image_ref = Images.imload()
    x, y = np.meshgrid(np.arange(0,np.size(image_ref,1),np.arange(0,np.size(image_ref,0))))
    image_c = np.zeros(np.size(image_ref,0),np.size(image_ref,1),len(Images))
    image_c[:,:,1] = Image_ref

    for i in len(Images)-1:

        im = Images.imload()
        x_shifts = shift_x[:, :, i]
        y_shifts = shift_y[:, :, i]

        # Correct rigid translation in x and y

        x_shifts_new = x_shifts[:] - np.mean(x_shifts[:])
        y_shifts_new = y_shifts[:] - np.mean(y_shifts[:])

        # Guess rotation - rotation is confined 90 degrees with moving rotation axis
        params0 = np.array([np.size(im,1)/8, np.size(im,0)/1, -2]) #initial [x_centre, y_centre, theta]
        up = np.array([np.size(im,1), np.size(im,0), 45]) #upper bound [x_centre, y_centre, theta]
        lb = np.array([1, 1,-45]) #lower bound [x_centre, y_centre, theta]
        
        params1 = rot_calc(x_shifts_new, y_shifts_new, pos_x, pos_y, params0, ub, lb)

        # Now we fit a model to determine the minimum rotation for the shifts
        # free parameters are centre of rotation and rotation angle (xc, yc, theta)

        # correct image rotation

        xc = params1[0] # Centre x
        yc = params1[1] # Centre y
        theta = params1[2] # Rotation theta
        rotation = np.array([np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)])

        # Apply corrections

        im_shift = np.array([x-np.mean(x_shifts), y-np.mean(y_shifts)])
        im_rot = np.array([(im_shift[:,0]-xc)])
    return image_c

def rot_calc(x_shift, y_shift, x_pos, y_pos, params0, ub, lb):
    
    fun = lambda params0 : rotation_fun(params0, x_shift, y_shift, x_pos, y_pos)

    params1 = least_squares(fun, params0, bounds=(ub, lb))
    
    return params1

def rotation_fun(params, x_shift, y_shift, x_pos, y_pos):
    # WEEEEEEEEEEEEE
    xc = params[0]
    yc = params[1]
    theta = params[2]

    # calculate rotation matrix

    rotm = np.array([np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)])
    rotm = np.degrees(rotm)

    # rotate existing grid about (xc, yc)

    new_points = np.array([x_pos[:]-xc, y_pos[:]-yc])*rotm

    # calculate shifts
    shifts = np.array([x_pos[:]-xc-new_points[:,0], y_pos[:]-yc-new_points[:,1]])

    # Find residuals

    resids = np.array([shifts[:,0]-x_shift],[shifts[:,1]-y_shift])
    resids = np.abs(resids)

    return resids