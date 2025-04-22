#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile

# for contours
import scipy.optimize as so
#phoutils isophote
from Snap import dvectors, deriv, graph


# In[ ]:


"""
This takes the density contour function from lab 7 to apply to our contour graphs

"""

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, 
         colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), density=True)
    # NOTE : if you are using the latest version of python, in the above: 
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    fmt = {}
    
    ### Adjust Here #### 
    
    # Contour Levels Definitions
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.25))
    new_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.50))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.66))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.85))
    
    
    # You might need to add a few levels


    # Array of Contour levels. Adjust according to the above
    levels = [one_sigma, new_sigma, two_sigma, three_sigma][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.25','0.50','0.66', '0.85'][::-1]

    
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour
    



# In[ ]:


"""
Here I am putting the stuff from lab 7 into one function to make things a little cleaner later on
In this case, I am creating a graph that determines the position of particles relavite to the center of mass of the galaxy
"""
def create_graph(filename,ptype):
    """
    This function creates a contour graph of a specific type of particle in the snapshot of a galaxy
    Inputs:
        filename: the name of the file that we are taking data from
        ptype: the type of particle that we are looking at
    """
    COMD = CenterOfMass(filename,ptype)
    COMP = COMD.COM_P(0.1)
    #COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
    # Determine positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 

    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)

    # Determine velocities of disk particles relatiev to COM motion
    #vxD = COMD.vx - COMV[0].value 
    #vyD = COMD.vy - COMV[1].value 
    #vzD = COMD.vz - COMV[2].value 

    # total velocity 
    #vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    #v = np.array([vxD,vyD,vzD]).T
    #graph(xD,yD)
        
    # M31 Disk Density 
    fig, ax= plt.subplots(figsize=(12, 10))

    # ADD HERE
    # plot the particle density for M31 using a 2D historgram
    # plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
    # cmap options: 
    # https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
    #   e.g. 'magma', 'viridis'
    # can modify bin number to make the plot smoother
    plt.hist2d(xD,yD,bins = 1500, norm = LogNorm(), cmap='viridis')
               #,range=[[-40, 40], [-40, 40]])
    #print(xD)
    cbar = plt.colorbar()
    cbar.set_label("Number of disk particle per bin", fontsize=15)

    # ADD HERE
    # make the contour plot
    # x pos, y pos, contour res, contour res, axis, colors for contours.
    # remember to adjust this if there are other contours added
    # density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])

    density_contour(xD,yD,80,80,ax = ax, colors=['yellow','red','green','cyan'])

    # Add axis labels
    plt.xlabel('x (kpc)', fontsize=22)
    plt.ylabel('y (kpc)', fontsize=22)

    #set axis limits
    plt.ylim(-400,400)
    plt.xlim(-400,400)

    #adjust tick label font size
    label_size = 22
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size



    # Save to a file
    #plt.savefig('Lab7_M31Disk.png')
    plt.show()
    


# In[ ]:


def create_graph2(filename,ptype):
    """
    This function is meant to return values defined in create_graph, without actually returning a graph
    This function returns the position and velocity data of a galaxy snap for a specific particle type.
    Inputs:
        filename: the name of the file that we are taking data from
        ptype: the type of particle that we are looking at
    Outputs:
        r: a transposed array containing the distance vector of all particles of a specific type 
        relative to a galaxies center of mass
        v: a transposed array containing the velocity vector of all 
        particles of a specific type relative to a galaxies center of mass velocity
    """
    COMD = CenterOfMass(filename,ptype)
    COMP = COMD.COM_P(0.1)
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
    # Determine positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 

    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)

    # Determine velocities of disk particles relatiev to COM motion
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 

    # total velocity 
    #vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T
    #graph(xD,yD)
    return r,v


# In[ ]:


def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) 
            such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) 
            such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    
    # normalize the angular momentum vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to
    # z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

