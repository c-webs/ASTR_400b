#!/usr/bin/env python
# coding: utf-8

# In[503]:


import numpy as np
import astropy.units as u
from ReadFile import Read


# In[565]:


def ComponentMass(PType,filename):
    """
    This function sums togeter the total value of a particle type in a galaxy
    Inputs:
        PType: The type of particle that this is
        filename: The file that we are extracting the data from
    Outputs:
        M_tot: the total mass of a type of particle in a galaxy
    """
    M_tot = 0
    time, particles, data = Read(filename)
    index = np.where(data['type'] == PType)
    mnew = data['m'][index]*1e10*u.Msun
    for Mass in mnew:
        M_tot = M_tot + Mass
    if M_tot != 0:
        M_tot = M_tot.to(1e12*u.Msun)
    return np.round(M_tot,3)
    # Use a for loop to add to each one as we go through the particles in each galaxy


# In[567]:





# In[569]:


def mtable(tablet):
    files = ["MW_000.txt","M31_000.txt","M33_000.txt"]
    for i in range(len(tablet)):
        for k in range(len(tablet[i])):
            if ComponentMass(i,files[k-1]) != 0:
                tablet[i][k] = ComponentMass(i,files[k-1]).value
            if ComponentMass(i,files[k-1]) == 0:
                tablet[i][k] = ComponentMass(i,files[k-1])
    return tablet
    


# In[571]:


# This sets up the table used for the mtable function
table2 = [["","","",""],
        ["","","",""],
        ["","","",""],
          ["","","",""]]
#Prints the final function
#Takes awhile to run (about a minute), but it does work
final_table = mtable(table2)
#print(final_table)


# In[572]:


del final_table[0] # gets ride of first row that is just zeros
for i in range(len(final_table)): # gets ride of the extra column
    if len(final_table[i]) >=3:
        del final_table[i][0]
# List of what needs to be the first column
galaxy_mass_type = ["Halo Mass (1e12 Msun)","Disk Mass (1e12 Msun)","Bulge Mass (1e12 Msun)"] 
galaxy_names = ["Galaxy Name","MW","M31","M33"] # List of galaxt names that will be the top row
for i in range(len(final_table)):
    final_table[i-1].insert(0,galaxy_mass_type[i-1]) # inserts first column
final_table.insert(0,galaxy_names) # inserts first row


# In[573]:


# Next 4 lines create the row of the total mass for each galaxy
totals = ["Total (1e12 Msun)","","",""]
totals[1] = final_table[1][1] + final_table[2][1] + final_table[3][1]
totals[2] = final_table[1][2] + final_table[2][2] + final_table[3][2]
totals[3] = final_table[1][3] + final_table[2][3] + final_table[3][3]
final_table.append(totals) # adds this row to the table
Local = ["Local Mass (1e12 Msun)",sum(totals[1:3])] # Next row is total local mass
final_table.append(Local) # Adds the row with total local mass



# In[574]:


# calculates F_bar values for each galaxy
f_bar_MW = (final_table[2][1] + final_table[3][1])/totals[1]
f_bar_M31 = (final_table[2][2] + final_table[3][2])/totals[2]
f_bar_M33 = (final_table[2][3] + final_table[3][3])/totals[3]
# puts the different f_bar values into a list
F_bar = ["f_bar",np.round(f_bar_MW,3),np.round(f_bar_M31,3),np.round(f_bar_M33,3)]
# adds F_bar to the table
final_table.append(F_bar)
# Print final_table
#print(final_table)


# In[575]:


# download final table
import pandas as pd
f = final_table


# In[576]:


# This makes the data much nicer looking 
df=pd.DataFrame(f[1:], columns = f[0])
#print (df)


# In[583]:




# In[585]:


import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
from reportlab.lib import colors


# In[587]:


pdf_filename = "table.pdf"
pdf = SimpleDocTemplate(pdf_filename, pagesize=letter)
table = Table([df.columns.tolist()] + df.values.tolist())
pdf.build([table])


# In[ ]:




