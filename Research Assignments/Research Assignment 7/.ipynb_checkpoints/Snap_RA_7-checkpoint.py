#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import astropy.units as u
#from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm

# my modules
#from ReadFile import Read
#from CenterOfMass import CenterOfMass
#from MassProfile import MassProfile

# for contours
#import scipy.optimize as so
#phoutils isophote


# In[199]:


dataM31 = np.genfromtxt("orbit_M31.txt",dtype=None,names=True)
dataM33 = np.genfromtxt("orbit_M33.txt",dtype=None,names=True)


# In[201]:



# In[203]:


"""
This function is from HW6, it is used to determine the distance between the center of masses of two galaxies as a function of time
in order to determine which points I want to look at the structure of M33, im going to need to have this data available first
"""

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit 

def dvectors(data1,data2):
    """
    This function loops over the datafiles for position and velocity vectors of two galaxies and creates an
    array that represents how the magnitude of the distance between the two galaxies changes over time
    Inputs:
        file1 ('np.ndarray'): Contains the time, position, and velocity data for one of the galaxies
        file2 ('np.ndarray'): Contains the time, position, and velocity data for the other galaxy

    Outputs:
        dis ('np.ndarray'): Contains the distance between two galaxies over time
        time ('np.ndarray'): Contains the time of the simulation
    """
    # gets the coordinates of the two galaxies and assigns them to their own arrays
    coord1 = data1['x'],data1['y'],data1['z']
    coord2 = data2['x'],data2['y'],data2['z']
    distance1 = np.array(coord1)
    distance2 = np.array(coord2)
    # makes an array containing the distances between each galaxy for each point in time
    distvec = ((distance1-distance2)**2).T
    distmag = np.sqrt(distvec[:,0]+distvec[:,1]+distvec[:,2])
    time = data1['t'] * 1e-3
    return distmag,time


# In[205]:


def deriv(dist,time):
    """
    This function takes the time and relative distances two galaxies to compute
    the snapshots and times where they are at pericenter and apocenter
    inputs:
        dist ('np.ndarray'): contains the distances between the galaxies during the simulation
        time ('np.ndarray'): contains the time of the simulation
    outputs:
        apo_snaps ('np.ndarray'): contains the snapshots where the two galaxies are at apocenter
        peri_snaps ('np.ndarray'): contains the snapshots where the two galaxies are at pericenter
        apo_time ('np.ndarray'): contains the times in Gyrs where the two galaxies are at apocenter
        per_time ('np.ndarray'): contains the times in Gyrs where the two galaxies are at pericenter
    """
    an = np.zeros(len(dist)) # initializes an array that will contain the snaps where M33 is at apocenter
    pn = np.zeros(len(dist))  # initializes an array that will contain the snaps where M33 is at pericenter
    apo_t = np.zeros(len(dist)) # initializes an array that will contain the times where M33 is at apocenter
    per_t = np.zeros(len(dist)) # initializes an array that will contain the times where M33 is at pericenter
    n = 0
    nstep = 'min' # tells the code whether it is looking for a perienter or apocenter point
    while n < len(time)-1:
        if nstep == 'max': # locates a maximum point on the graph and adds that snap and time to an and apo_t
            if dist[n] > dist[n+1]:
                apo_t[n] = (time[n])
                an[n] = n
                nstep = 'min'
        if nstep == 'min': # locates a minimum point on the graph and adds that snap and time to pn and peri_t
            if dist[n] < dist[n+1]:
                per_t[n]=time[n]
                pn[n] = n
                nstep = 'max'
        n = n+1
    # eliminates the zeros in an, pn, apo_t, and peri_t to make the arrays contain the data that we need
    apo_index = np.where(an != 0)
    peri_index = np.where(pn != 0)
    # multiplies the an and pn values by 5 to make them match the snaps. 
    # this is because the files being used contain intervals of 5 (row 100 actually represents snap 500)
    apo_snaps = (an[apo_index])*5
    peri_snaps = (pn[peri_index])*5
    apo_time = apo_t[apo_index]
    per_time = per_t[peri_index]
    return apo_snaps, peri_snaps, apo_time, per_time


# In[207]:


dist, time = dvectors(dataM31,dataM33)


# In[209]:


apo_snaps, peri_snaps, at, pt = deriv(dist,time)


# In[301]:


"""
Here I am graphing my pericenter and apocenter points over the graph from Homework 6
to confirm that my code works
"""
def graph(distances,times, peri_times, apo_times):
    """
    This function graphs the distance between two galaxies over time, and then inserts vertical lines
    at peri_times and apo_times to show that deriv is working properly
    inputs:
        distances ('np.ndarray'): Contains the distance between two galaxies over time
        times ('np.ndarray'): Contains the time of the simulation
        peri_times ('np.ndarray'): contains the times in Gyrs where the two galaxies are at apocenter
        apo_times: ('np.ndarray'): contains the times in Gyrs where the two galaxies are at pericenter
    """
    for t in range(len(apo_times)): # graphs the apocenter and pericenter points
        if t == 0: # makes sure that apocenter and pericenter don't show up on the legend more than once
            plt.axvline(x = apo_times[t], linestyle = '-', color = 'blue', label = 'apocenter')
            plt.axvline(x = peri_times[t], linestyle = '-', color = 'green', label = 'pericenter')
        if t > 0:
            plt.axvline(x = apo_times[t], linestyle = '-', color = 'blue')
            plt.axvline(x = peri_times[t], linestyle = '-', color = 'green')
    plt.plot(times,distances, linestyle='-', color='red', label = 'simulation') # graphs the distances between the two galaxies as a function of time
    plt.xlabel("Time (Gyrs)")
    plt.ylabel("Distance (kpc)")
    plt.title("Distance vs. Time")
    plt.grid(True)
    plt.legend()
    plt.show()


# In[303]:


#graph(M3133[0],M3133[2],pt,at)


# In[ ]:




