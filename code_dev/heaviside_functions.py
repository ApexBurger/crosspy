# Heaviside functions

class Point: 
    def __init__(self, x, y): 
        self.x = x 
        self.y = y 
    
def line_test(pixel, r, theta, imarray):
    # Find start and end points of line segments for different cases of theta and r
    dx = imarray.shape[0]
    # Line segment to interogation point
    p1 = Point(0, 0)
    q1 = Point(pixel[0], pixel[1])
    
    # Vector magnitude cases
    
    if r == 0:
        r = 1e-8
        
    # Rotation cases
    if theta == 0. or theta == 360.: # vertical to right
        x1 = r
        x2 = q1.x
        if x2 > x1:
            return False
        else:
            return True
    elif theta == 90.: # horizontal line above
        y1 = r
        y2 = q1.y
        if y2>y1:
            return False
        else:
            return True
    elif theta == 180.: # vertical to left
        x1 = -r
        x2 = q1.x
        if x2 > x1:
            return True
        else:
            return False
    elif theta == 270.: # horizontal below
        y1 = -r
        y2 = q1.y
        if y2 < y1:
            return False
        else:
            return True
    elif theta>0 and theta<180:
        theta = np.radians(theta)
        # Tangent line segment
        t1 = Point(r*np.cos(theta), r*np.sin(theta))
        m = -1*(np.cos(theta)/np.sin(theta))
        c = t1.y - m*t1.x
        y1 = q1.y
        y2 = m*q1.x + c
        if y1>y2:
            return False
        else:
            return True
    elif theta>180 and theta<360:
        theta = np.radians(theta)
        # Tangent line segment
        t1 = Point(r*np.cos(theta), r*np.sin(theta))
        m = -1*(np.cos(theta)/np.sin(theta))
        c = t1.y - m*t1.x
        
        y1 = q1.y
        y2 = m*q1.x + c
        if y1<y2:
            return False
        else:
            return True
        
    
#@jit(nopython=True)
def subset_hsfilter(imarray, r, theta):
    if type(imarray) == list:
        imarray = imarray[0]
    # preallocate
    hsfilter = np.zeros((imarray.shape[0],imarray.shape[0]))
    xc = hsfilter.shape[0]/2
    yc = hsfilter.shape[1]/2
    xc = int(xc)
    yc = int(yc)
    
    # Create x and y coordinates which are centred
    xs,ys = np.meshgrid(range(-xc, xc), range(-yc,yc))
    
    
    for col in range(imarray.shape[0]):
        for row in range(imarray.shape[0]):
            #rasters through columns and rows for a given coordinate in xy
            x = xs[row,col]
            y = ys[row,col]
            # Note that y axis is mirrored
            pixel = [x, (-1*y)]
            
            # Test if pixel is beyond the discontinuity line
            if line_test(pixel, r, theta, imarray):
                hsfilter[row,col] = 1
            else:
                hsfilter[row,col] = False
                
            
                
    hs_subset = np.zeros(imarray.shape)
    hs_subset[:,:,0] = np.multiply(hsfilter,imarray[:,:,0])
    hs_subset[:,:,1] = np.multiply(hsfilter,imarray[:,:,1])
    return hs_subset

def hs_corr(x, imarrays):
    
    r, theta = x
    
    filtered_subsets = subset_hsfilter(imarrays, r, theta)
    
    a = filtered_subsets[:,:,0]
    b = filtered_subsets[:,:,1]
    
    result = XCF.fxcorr(a, b, dic_disc, XCF.plan_ffts(dic_disc))
    
    act_px = np.count_nonzero(a != False)
    result = np.asarray(result)
    result[2] = result[2]/act_px
    
    return result

def fit_samples(ia, n=5):
    
    def mapspace(r_range,theta_range, ia):
        dxs = np.zeros((len(r_range), len(theta_range)))
        dys = np.zeros((len(r_range), len(theta_range)))
        ccs = np.zeros((len(r_range), len(theta_range)))


        for i, r in enumerate(r_range):
            for j, theta in enumerate(theta_range):
                x = [r,theta]
                c = hs_corr(x, ia)
                ccs[i,j] = c[2]
        return dxs, dys, ccs
    
    # step 1 do a coarse search for peak height
    # set range for r and theta
    rs = np.linspace(0,ia.shape[1]/2, n) # r
    thetas = np.linspace(0,360, n) # theta
    
    # coarse search
    _,_,z = mapspace(rs,thetas,ia)
    
    # find peak
    loc = np.where(z == z.max())

    # finer search space
    
    ind1 = loc[0]-1 % len(r_range)
    ind2 = loc[0]+1 % len(r_range)
    ind3 = loc[0]-1 % len(theta_range)
    ind
    r_fine = np.linspace(r_range[ind1],r_range[loc[0]+1],5)
    theta_fine = np.linspace(theta_range[loc[1]-1], theta_range[loc[1]+1],5)

    # new points

    _,_,zf = mapspace(r_fine,theta_fine,ia)

    # find peak - best fit of r and theta
    loc = np.where(zf == zf.max())
    r = r_fine[loc[1]]
    theta = theta_fine[loc[0]]
    return r, theta