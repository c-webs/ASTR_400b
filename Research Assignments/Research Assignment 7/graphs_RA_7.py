#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt

# my modules
#from ReadFile import Read
from CenterOfMass import CenterOfMass
#from MassProfile import MassProfile


import scipy.optimize as so

# In[9]:


def get_COM(filename,ptype):
    """
    This function returns the position and velocity data of a galaxy snap for a specific particle type
    after redifining the coordinate system to contain the center of mass of the galaxy at the origin point.
    Inputs:
        filename ('string'): the name of the file that we are taking data from
        ptype ('int'): the type of particle that we are looking at
    Outputs:
        r ('np.ndarray'): a transposed array containing the distance vector of all particles of a specific type 
        relative to a galaxies center of mass
        v ('np.ndarray'): a transposed array containing the velocity vector of all 
        particles of a specific type relative to a galaxies center of mass velocity
    """
    # Computes the center of mass and center of mass moton of the galaxy
    COMD = CenterOfMass(filename,ptype)
    COMP = COMD.COM_P(0.1)
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
    # Determines the positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 

    # Determines the velocities of disk particles relatiev to COM motion
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 


    # gets transposed arrays for the vector location of each particle relative to the center of mass
    r = np.array([xD,yD,zD]).T # transposed 
    # gets transposed arrays for the vector velocity of each particle relative to the center of mass velocity
    v = np.array([vxD,vyD,vzD]).T
    return r,v


# In[11]:


def RotateFrame(posI,velI):
    """a function that rotates the position and velocity vectors
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
    
    # computes the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    
    # normalizes the angular momentum vector
    L_norm = L/np.sqrt(np.sum(L**2))
    
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

