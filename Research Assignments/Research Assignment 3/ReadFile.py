#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import astropy.units as u
# Imports numpy and Astropy


# In[101]:


def Read(filename):
    file = open(filename,'r')
    line1 = file.readline()
    label,value = line1.split()
    # Reads the first line of the file and assigns it to the variable, value
    time = float(value)*u.Myr
    # Gives value units of Myr
    line2 = file.readline()
    label2,value2 = line2.split()
    particles = int(value2)
    # Gives the particle
    file.close()
    # Closes the file
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    # assigns the data from the 4th row on to a data array, with the headers being taken from the third row
    return time, particles, data



# In[103]:


value,particles , data = Read("MW_000.txt")
# Allows us to return the different variables of Read as a tuple
#print(value)
#print(particles)
#print(data)


# In[75]:





# In[ ]:




