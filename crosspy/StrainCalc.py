# Strain calculation functions
import numpy as np
from numpy.polynomial.polynomial import polyder, polyval, polyfit
import crosspy

def strain_calc(d, mapnos = 0, strain_method = 'l2'):
    # This function calculates the strain between two images
    #   Images
    #   nodes_x = roi centres x
    #   nodes_y = roi centres y
    #   disp_x = dx_maps
    #   disp_y = dy_maps

    # get displacements and positions
    x = d.x_pos
    y = d.y_pos
    dx = d.dx_maps[:,:,mapnos]
    dy = d.dy_maps[:,:,mapnos]
    # determine strain method
    if strain_method == '9nodes':
        e, e_eff, r, f = strain_n9()
    elif strain_method == '8nodes':
        e, e_eff, r, f = strain_n8()
    elif strain_method == '4nodes':
        e, e_eff, r, f = strain_n4()
    elif strain_method == 'l2':
        e, e_eff, r, f = strain_l2(dx, dy, x, y)
    else:
        raise Exception('Invalid strain method!')


    return e, e_eff, r, f
    
def strain_n4():
    raise Exception('U A KEEN BEAN - WORK IN PROGRESS')

def strain_n8():
    raise Exception('U A MEAN TEEN - WORK IN PROGRESS')

def strain_n9():
    raise Exception('U A LEAN MACHINE - WORK IN PROGRESS')

def strain_l2(dx, dy, x, y):
    # Function to determine strain via polynomial fitting and l2 min
    # Preallocate arrays
    rows = np.size(dx,0)
    cols = np.size(dx,1)
    e_temp = np.zeros(shape=(rows,cols,3,3))
    eeff_temp = np.zeros(shape=(rows,cols,1))
    rotation_temp = np.zeros(shape=(rows,cols,3,3))
    e11_temp = np.zeros(shape=(rows,cols))
    e22_temp = np.zeros(shape=(rows,cols))
    e12_temp = np.zeros(shape=(rows,cols))
    
    # First obtain strain values for corners and edges of the map
    e11_temp, e22_temp, e12_temp = strain_l2_corners(dx, dy, x, y, e11_temp, e22_temp, e12_temp)
    e11_temp, e22_temp, e12_temp = strain_l2_edges(dx, dy, x, y, e11_temp, e22_temp, e12_temp)

    # Obtain bulk values in loop below - corners and edges are already determined
    for i in range(0,rows-1):
        for j in range(0,cols-1):
            dhor_3pt_x = np.array([dx[i,j-1], dx[i,j], dx[i,j+1]])
            dhor_3pt_y = np.array([dy[i,j-1], dy[i,j], dy[i,j+1]])
            dver_3pt_x = np.array([dx[i-1,j], dx[i,j], dx[i+1,j]])
            dver_3pt_y = np.array([dy[i-1,j], dy[i,j], dy[i+1,j]])
            pos_x3 = np.array([x[i,j-1], x[i,j], x[i,j+1]])
            pos_y3 = np.array([y[i-1,j], y[i,j], y[i+1,j]])

            # Determine second order poly fit and derivative
            coef_x = polyder(polyfit(pos_x3, dhor_3pt_x,2))
            coef_y = polyder(polyfit(pos_y3, dver_3pt_y,2))
            coef_xy = polyder(polyfit(pos_x3, dhor_3pt_y,2))
            coef_yx = polyder(polyfit(pos_y3, dhor_3pt_x,2))

            # Obtain values of polynomial fit at the centre of object pixel
        
            du_dx = polyval(coef_x, x[i,j]) # du from dx map
            dv_dy = polyval(coef_y, y[i,j]) # dv from dy map
            du_dy = polyval(coef_xy, x[i,j]) # du from dy map
            dv_dx = polyval(coef_yx, y[i,j]) # dv from dx map


            # Create the deformation gradient F from displacements u and v
            F = np.array([[du_dx, du_dy, 0], [dv_dx, dv_dy, 0], [0, 0, -(du_dx+dv_dy)]])+np.eye(3)
            Ft = F.transpose()
            C_test = np.dot(F,Ft)
            C = np.matmul(F.transpose(), F) # Green-Lagrange tensor
            V, Q = np.linalg.eig(C) # eigenvalues V and vector Q
            V = np.diag(V)
            Ut = np.sqrt(V)
            U = np.matmul(Q.transpose(), np.matmul(Ut, Q))
            U_1 = np.linalg.inv(U)
            R = np.dot(F, U_1) # rotations

            # Determine green strain tensor and rotation tensor from F by symm and anti symm parts
            e_temp[i,j,:,:] = 0.5*(np.matmul(F.transpose(),F-np.eye(3)))
            rotation_temp[i,j,:,:] = R
            # Determine the effective strain
            xs = np.dot(e_temp[i,j],e_temp[i,j])
            eeff_temp[i,j,:] = np.sqrt((2/3)*np.tensordot(e_temp[i,j],e_temp[i,j]))

    # Form outputs
    strain = e_temp
    strain_effective = eeff_temp
    rotation = R
    
    return strain, strain_effective, rotation, F

    
def strain_l2_corners(dx, dy, x, y, e11, e22, e12):
    # Use first order polynomial fitting for the 4 corners of the map
    rows = np.size(dx,0)-1
    cols = np.size(dx,1)-1
    
    ## first corner - top left
    dx_x1 = [dx[0,0], dx[0,1]]
    dx_y1 = [dy[0,0], dy[0,0]]
    dy_x1 = [dx[0,0], dx[1,1]]
    dy_y1 = [dy[0,0], dy[1,0]]
    pos_x1 = [x[0,0], x[0,1]]
    pos_y1 = [y[0,0], y[1,0]]
    # determine first order of poly fit
    coef_x = polyfit(pos_x1, dx_x1,1)
    coef_y = polyfit(pos_y1, dy_y1,1)
    coef_xy = polyfit(pos_x1, dx_y1,1)
    coef_yx = polyfit(pos_y1, dy_x1,1)
    # equal to strain
    e11[0,0] = coef_x[0]
    e22[0,0] = coef_y[0]
    e12[0,0] = 0.5*(coef_xy[0]+coef_yx[0])

    ## second corner - top right
    dx_x2 = [dx[0,cols-1], dx[0,cols]]
    dx_y2 = [dy[0,cols], dy[1,cols]]
    dy_x2 = [dx[0,cols-1], dx[0,cols]]
    dy_y2 = [dy[0,cols], dy[1,cols]]
    pos_x2 = [x[0,cols-1], x[1,cols]]
    pos_y2 = [y[0,cols], y[1,cols]]
    # determine first order of poly fit
    coef_2 = polyfit(pos_x2, dx_x2,1)
    coef_2 = polyfit(pos_y2, dy_y2,1)
    coef_xy = polyfit(pos_x2, dx_y2,1)
    coef_yx = polyfit(pos_y2, dy_x2,1)
    # equal to strain
    e11[0,cols] = coef_x[0]
    e22[0,cols] = coef_y[0]
    e12[0,cols] = 0.5*(coef_xy[0]+coef_yx[0])

    ## third corner - bottom left
    dx_x3 = [dx[rows,0], dx[rows,1]]
    dx_y3 = [dy[rows-1,0], dy[rows,0]]
    dy_x3 = [dx[rows,0], dx[rows,1]]
    dy_y3 = [dy[rows-1,0], dy[rows,0]]
    pos_x3 = [x[rows,0], x[rows,1]]
    pos_y3 = [y[rows-1,0], y[rows,0]]
    # determine first order of poly fit
    coef_x = polyfit(pos_x3, dx_x3,1)
    coef_y = polyfit(pos_y3, dy_y3,1)
    coef_xy = polyfit(pos_x3, dx_y3,1)
    coef_yx = polyfit(pos_y3, dy_x3,1)
    # equal to strain
    e11[rows,0] = coef_x[0]
    e22[rows,0] = coef_y[0]
    e12[rows,0] = 0.5*(coef_xy[0]+coef_yx[0])
    ## fourth corner - bottom right
    dx_x4 = [dx[rows,cols-1], dx[rows,cols]]
    dx_y4 = [dy[rows-1,cols], dy[rows,cols]]
    dy_x4 = [dx[rows,cols-1], dx[rows,cols]]
    dy_y4 = [dy[rows-1,cols], dy[rows,cols]]
    pos_x4 = [x[rows,cols-1], x[rows,cols]]
    pos_y4 = [y[rows-1,cols], y[rows,cols]]
    # determine first order of poly fit
    coef_x = polyfit(pos_x4, dx_x4,1)
    coef_y = polyfit(pos_y4, dy_y4,1)
    coef_xy = polyfit(pos_x4, dx_y4,1)
    coef_yx = polyfit(pos_y4, dy_x4,1)
    # equal to strain
    e11[rows,cols] = coef_x[0]
    e22[rows,cols] = coef_y[0]
    e12[rows,cols] = 0.5*(coef_xy[0]+coef_yx[0])

    return e11, e22, e12
    
def strain_l2_edges(dx, dy, x, y, e11, e22, e12):
    # Use polynomial fit to find strain on edges of map
    rows = np.size(dx,0)-1
    cols = np.size(dx,1)-1

    # Top edge
    for i in [0]:
        for j in range(2,cols):
            dx_3pt_x = [dx[i,j-1], dx[i,j], dx[i,j+1]]
            dy_2pt_y = [dy[i,j], dy[i+1,j]]
            dx_3pt_y = [dy[i,j-1], dy[i,j], dy[i,j+1]]
            dy_2pt_x = [dx[i,j], dx[i+1,j]]
            x_3pt = [x[i,j-1], x[i,j], x[i,j+1]]
            y_2pt = [y[i,j], y[i+1,j]]

            coef_x = polyder(polyfit(x_3pt, dx_3pt_x, 2))
            coef_y = polyfit(y_2pt, dy_2pt_y, 1)
            coef_xy = polyder(polyfit(x_3pt, dx_3pt_y, 2))
            coef_yx = polyfit(y_2pt, dy_2pt_x, 1)
            e11[i,j] = polyval(coef_x,x[i,j])
            e22[i,j] = coef_y[0]
            e12[i,j] = 0.5*(polyval(coef_xy, x[i,j])+coef_yx[0])

    # Bottom edge

    for i in [rows]:
        for j in range(2,cols):
            dx_3pt_x = [dx[i,j-1], dx[i,j], dx[i,j+1]]
            dy_2pt_y = [dy[i-1,j], dy[i,j]]
            dx_3pt_y = [dy[i,j-1], dy[i,j], dy[i,j+1]]
            dy_2pt_x = [dx[i-1,j], dx[i,j]]
            x_3pt = [x[i,j-1], x[i,j], x[i,j+1]]
            y_2pt = [y[i-1,j], y[i,j]]

            coef_x = polyder(polyfit(x_3pt, dx_3pt_x, 2))
            coef_y = polyfit(y_2pt, dy_2pt_y, 1)
            coef_xy = polyder(polyfit(x_3pt, dx_3pt_y, 2))
            coef_yx = polyfit(y_2pt, dy_2pt_x, 1)
            e11[i,j] = polyval(coef_x,x[i,j])
            e22[i,j] = coef_y[0]
            e12[i,j] = 0.5*(polyval(coef_xy, x[i,j])+coef_yx[0])

    # Left edge

    for i in range(2,rows):
        for j in [0]:
            dx_2pt_x = [dx[i,j], dx[i,j+1]]
            dy_3pt_y = [dy[i-1,j], dy[i,j], dy[i+1,j]]
            dx_2pt_y = [dy[i,j], dy[i,j+1]]
            dy_3pt_x = [dx[i-1,j], dx[i,j], dx[i+1,j]]
            x_2pt = [x[i,j], x[i,j+1]]
            y_3pt = [y[i-1,j], y[i,j], y[i+1,j]]

            coef_x = polyfit(x_2pt, dx_2pt_x,1)
            coef_y = polyder(polyfit(y_3pt, dy_3pt_y, 2))
            coef_xy = polyfit(x_2pt, dx_2pt_y, 1)
            coef_yx = polyder(polyfit(y_3pt, dy_3pt_x, 2))
            e11[i,j] = coef_x[0]
            e22[i,j] = polyval(coef_y,y[i,j])
            e12[i,j] = 0.5*(coef_xy[0] + polyval(coef_yx, y[i,j]))

    # Right edge 

    for i in range(2,rows):
        for j in [cols]:
            dx_2pt_x = [dx[i,j-1], dx[i,j]]
            dy_3pt_y = [dy[i-1,j], dy[i,j], dy[i+1,j]]
            dx_2pt_y = [dy[i,j-1], dy[i,j]]
            dy_3pt_x = [dx[i-1,j], dx[i,j], dx[i+1,j]]
            x_2pt = [x[i,j-1], x[i,j]]
            y_3pt = [y[i-1,j], y[i,j], y[i+1,j]]

            coef_x = polyfit(x_2pt, dx_2pt_x,1)
            coef_y = polyder(polyfit(y_3pt, dy_3pt_y, 2))
            coef_xy = polyfit(x_2pt, dx_2pt_y, 1)
            coef_yx = polyder(polyfit(y_3pt, dy_3pt_x, 2))
            e11[i,j] = coef_x[0]
            e22[i,j] = polyval(coef_y,y[i,j])
            e12[i,j] = 0.5*(coef_xy[0] + polyval(coef_yx, y[i,j]))

        return e11, e22, e12