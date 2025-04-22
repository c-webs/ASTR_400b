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


# In[ ]:


"""
This function is from HW6, it is used to determine the distance between the center of masses of two galaxies as a function of time
in order to determine which points I want to look at the structure of M33, im going to need to have this data available first
"""

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit 

def dvectors(file1,file2):
    """
    This function loops over the datafiles for position and velocity vectors of two galaxies and creates 2
    lists that represents how the magnitudes of the distance and velocity between the two galaxies change over time
    Inputs:
        file1: the file containing the time, position, and velocity data for one of the galaxies
        file2: the file containing the time, position, and velocity data for the other galaxy

    Outputs:
        dis: A list containing the distance between two galaxies over time
        vel: A list containing the velcoity between two galaxies over time
    """
    dist = np.zeros(len(file1))
    vdist = np.zeros(len(file1))
    time = np.zeros(len(file1))
    for i in range(len(file1)):
        # Gathers position and velocity data for each point in time 
        x1 = file1[i][1]
        y1 = file1[i][2]
        z1 = file1[i][3]
        vx1 = file1[i][4]
        vy1 = file1[i][5]
        vz1 = file1[i][6]
        x2 = file2[i][1]
        y2 = file2[i][2]
        z2 = file2[i][3]
        vx2 = file2[i][4]
        vy2 = file2[i][5]
        vz2 = file2[i][6]
        t1 = file1[i][0] *1e-3
        # Computes the magnitudes of distance and velocity for each galaxy and stores them
        # as a data point on an array
        dist[i] = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        vdist[i] = np.sqrt((vx1-vx2)**2+(vy1-vy2)**2+(vz1-vz2)**2)
        time[i] = t1
        #converts arrays to lists
        tim = time.tolist()
        dis = dist.tolist()
        vel = vdist.tolist()
    return dis,vel,tim


# In[ ]:


def deriv(dist,vel,time):
    """
    This function takes the file created in distance and computes areas where
    the derivative of the distance is roughly zero. 
    (It's worth noting that this is an extremely crude approximation, but as shown later, it works in this case)
    It also computes the midpoints in time between these derivatives
    for now, I am using the data from the .txt files created in homework 6, where we stored the center of mass of the galaxies
    at several different intervals. 
    After finding the indexes that I am looking for, and separating them into lists for apocenter, pericenter, and midpoint,
    I multiply the index values by 5 to get the snapshot
    note: 5 is the interval in which I wrote certain snapshots into the .txt file from homework 6
    """
    #dist, vel, time = 
    apo = []
    an = []
    pn = []
    per = []
    apo_t = []
    per_t = []
    index = []
    points = []
    points_t = []
    #distance = []
    #time = []
    n = 0
    nstep = 'min'
    while n < len(time)-1:
        if nstep == 'max':
            if dist[n] > dist[n+1]:
                apo.append(dist[n])
                points.append(dist[n])
                apo_t.append(time[n])
                points_t.append(time[n])
                index.append(n)
                an.append (n)
                nstep = 'min'
        if nstep == 'min':
            if dist[n] < dist[n+1]:
                per.append(dist[n])
                points.append(dist[n])
                per_t.append(time[n])
                points_t.append(time[n])
                index.append(n)
                pn.append(n)
                nstep = 'max'
        n = n+1
    midpoints = []
    m = 0
    while m < len(index)-1:
        midpoints.append(int(np.round((index[m]+(index[m+1]))/2,0)))
        m = m+1
    apo_snaps = np.array(an)*5
    peri_snaps = np.array(pn)*5
    mid_snaps = np.array(midpoints)*5
    return apo_snaps, peri_snaps, mid_snaps, apo_t, per_t
    #return apo,apo_t,per,per_t,index,points,points_t,an,pn


# In[1]:


"""
Here I am graphing my pericenter and apocenter points over the graph from Homework 6
to confirm that my code works
"""
def graph(distances,times, peri_times, apo_times):
    distances2 = distances
    #M3133[0] # gets the distances from M3133
    time_2 = times
    #M3133[2] # gets the times from M3133
    for i in apo_times:
        j = [i]
        k = j * len(distances2)
        plt.plot(k,distances2, linestyle='-', color='blue', label = 'apocenter')
    for i in peri_times:
        j = [i]
        k = j * len(distances2)
        plt.plot(k,distances2, linestyle='-', color='green', label = 'pericenter')
    plt.plot(time_2,distances2, linestyle='-', color='red', label = 'simulation')
    #plt.plot(tesst2,distances2, linestyle='-', color='blue', label = 'test') 
    plt.xlabel("Time (Gyrs)")
    plt.ylabel("Distance (kpc)")
    plt.title("Distance vs. Time")
    #plt.ylim(0,800)
    plt.grid(True)
    #plt.legend()
    plt.show()


# In[ ]:




