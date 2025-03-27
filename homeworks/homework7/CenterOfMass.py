#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Homework 4
# Center of Mass Position and Velocity
# Colin Weber


# In[ ]:


# remember this is just a template,
# you don't need to follow every step.
# If you have your own method to solve the homework,
# it is totally fine


# In[5]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[ ]:





# In[82]:


class CenterOfMass:
# Class to define COM position and velocity properties 
# of a given galaxy and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        
    
    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        # x component Center of mass
        a_com = np.sum(a*m) / np.sum(m)
        # y component Center of mass
        b_com = np.sum(b*m) / np.sum(m)
        # z component Center of mass
        c_com = np.sum(c*m) / np.sum(m)
        # returns the 3 components separately
        return a_com,b_com,c_com
    
    def COM_P(self, delta):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     
        # calls COMdefine to try a first guess at COM                                         
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # computes the magnitude of the COM position vector.
        r_COM = np.sqrt(x_COM**2+y_COM**2+z_COM**2)

         # finds the max 3D distance of all particles from the guessed COM   
    
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2.0+y_new**2.0+z_new**2.0)

        # finds the max 3D distance of all particles from the guessed COM                                               
        # re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/2.0
        
        # initial value for the change in COM position                                                      
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # selects all particles within the reduced radius (starting from original x,y,z, m)
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]                                                              
            # computes the center of mass position using the particles in the reduced radius
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2,y2,z2,m2)
            # computes the new 3D COM position with x_COM2,y_COM2,z_COM2
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)                                 
            # Calculates the difference in position between the previous center of mass and the current one                                                                                        
            change = np.abs(r_COM - r_COM2)                    
            # reduces the volume by a factor of 2 again                                                                 
            r_max /= 2.0                                                                                    
            # changes the frame of reference to the newly computed COM by subtracting the new COM
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
            # sets the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # creates an array (np.array) to store the COM position and rounds all values                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])
        return np.round(p_COM, 2)*u.kpc
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # max distance from the center that is used to determine the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determines the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        xV = self.x[:]*u.kpc - x_COM
        yV = self.y[:]*u.kpc - y_COM
        zV = self.z[:]*u.kpc - z_COM
        rV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determines the index for those particles within the max radius
        indexV = np.where(rV < rv_max)
        # determines the velocity and mass of those particles within the mas radius
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # computes the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new,vy_new,vz_new, m_new)
        
        v_COM = np.array([vx_COM, vy_COM, vz_COM])                                                                                       
        return np.round(v_COM, 2)*u.km/u.s


# In[84]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MW_COM = CenterOfMass("MW_000.txt", 2)
M31_COM = CenterOfMass("M31_000.txt", 2)
M33_COM = CenterOfMass("M33_000.txt", 2)


# In[85]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
MW_COM_p = MW_COM.COM_P(0.1)
#print(MW_COM_p)
MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
#print(MW_COM_v)


# In[88]:


# now write your own code to answer questions
M31_COM_p = M31_COM.COM_P(0.1)
M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
#print(M31_COM_p, M31_COM_v)
#print(M31_COM_v)


# In[90]:


M33_COM_p = M33_COM.COM_P(0.1)
#print(M33_COM_p)
M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
#print(M33_COM_v)
#print(M33_COM_v[0])


# In[45]:


# Q.1: 
#MW COM position vector [-1,2,-1]
#MW COM velocity vector [-1,2,-1]
#M31 COM position vector [-378,611,-285]
#M31 COM velocity vector [74,-72,49]
#MW COM position vector [-476,491,-412]
#MW COM velocity vector [44,102,142]


# In[77]:


# Q.2
# Position separation
COMPD = np.sqrt((M31_COM_p[0]-MW_COM_p[0])**2+(M31_COM_p[1]-MW_COM_p[1])**2+(M31_COM_p[2]-MW_COM_p[2])**2)
#print (np.round(COMPD,3)*u.kpc)
# Velocity separation
COMVD = np.sqrt((M31_COM_v[0]-MW_COM_v[0])**2+(M31_COM_v[1]-MW_COM_v[1])**2+(M31_COM_v[2]-MW_COM_v[2])**2)
#print (np.round(COMVD,3)*u.km/u.s)


# In[79]:


# Q.3
# Position separation
COMPD33 = np.sqrt((M31_COM_p[0]-M33_COM_p[0])**2+(M31_COM_p[1]-M33_COM_p[1])**2+(M31_COM_p[2]-M33_COM_p[2])**2)
#print (np.round(COMPD33,3)*u.kpc)
# Velocity separation
COMVD33 = np.sqrt((M31_COM_v[0]-M33_COM_v[0])**2+(M31_COM_v[1]-M33_COM_v[1])**2+(M31_COM_v[2]-M33_COM_v[2])**2)
#print (np.round(COMVD33,3)*u.km/u.s)


# In[81]:


# Q.4
#This is important to helping us understand how two large masses interect, 
#it also allows for us to ignore the dark matter halo, 
#which is much less predictable


# In[60]:


M31_COM_p = M31_COM.COM_P(0.1)
M31_COM_v = M31_COM.COM_V(M31_COM_p[0],M31_COM_p[1],M31_COM_p[2])
#print('M31 COM xyz position:', M31_COM_p, 'and xyz velocity:', MW_COM_v)


# In[ ]:




