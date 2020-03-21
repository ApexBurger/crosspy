# Strain calculation functions

def strain_l2(shift_x, shift_y, pos_x, pos_y, im_size, im_num, im_interval):
    pass

def strain_l2_corners(displacements_x, displacements_y, pos_x, pos_y, e_11=0, e_22=0, e_12=0):
    if e_11 == 0 and e22 == 0 and e12 == 0:
        e11 = np.zeros(displacements_x.shape)
        e22 = e11
        e12 = e11
    
    