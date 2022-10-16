#!/usr/bin/env python
# coding: utf-8

# In[2]:


#import libraries
import numpy as np
import pandas as pd
import math
import scipy.linalg as la
import matplotlib.pyplot as plt
from scipy.stats import gamma


# In[3]:


print("Enter Depth in Feet")
d=float(input('Depth: '))


# In[4]:


print("Enter in-situ Stress Values in Psi:")
sigmah= float(input("sigmah: "))
sigmaH= float(input("sigmaH: "))
sigmav= float(input("sigmav: "))


# In[5]:


print("Enter in-situ stress Direction:")
alphas1= float(input("Azimuth: "))
betas1= float(input("Deviation: "))
alphas= alphas1*0.0175
betas= betas1*0.0175



# In[6]:


print("Enter Borehole Azimuth and Deviation:")
alphab1= float(input("Azimuth: "))
betab1= float(input("Direction: "))
alphab= alphab1*0.0175
betab= betab1*0.0175


# In[7]:


print("Enter Direction of Plane of Weakness")
alphaw1= float(input("Azimuth: "))
betaw1= float(input("Deviation: "))
alphaw= alphaw1*0.0175
betaw= betaw1*0.0175


# In[8]:


print("Enter Pore Pressure in psi:")
pp= float(input(" "))
print("Enter Biot's Coefficient:")
alpha= float(input(" "))
print("Enter Poission's Ratio:")
mu= float(input(" "))


# In[9]:


print("Enter the Rock Strength parameters:")
si= float(input("Cohesion of Rock Matrix in Psi: "))
mi= float(input("Coefﬁcient of Friction of Rock Matrix: "))
sw= float(input("Cohesion of Weak Plane in Psi: "))
mw= float(input("Coefﬁcient of Friction of Weak Plane: "))


# In[10]:


#Defining values for ICS to GCS
cosa= math.cos(alphas)
cosb= math.cos(betas)
sina= math.sin(alphas)
sinb= math.sin(betas)
sigmaics= np.array([[sigmah,0,0], [0,sigmaH,0], [0,0,sigmav]])


# In[11]:


#rotation Matrix for ICS to GCS
E= np.array([[cosa*cosb, sina*cosb, sinb], [-sina, cosa, 0], [-cosa*sinb, -sina*sinb, cosb]])
ET= np.transpose(E)


# In[12]:


#stresstransformationmatrix
sigmaics2ecs= np.dot(np.dot(ET, sigmaics), E)


# In[13]:


#Rotation matrix
B= np.array([[math.cos(alphab)*math.cos(betab), math.sin(alphab)*math.cos(betab), math.sin(betab)], 
            [-math.sin(alphab), math.cos(alphab), 0],
            [-math.cos(alphab)*math.sin(betab), -math.sin(alphab)*math.sin(betab), math.cos(betab)]])

BT= np.transpose(B)


# In[14]:


#Stress transformation matrix
sigmaecs2bcs= np.dot(np.dot(B, sigmaics2ecs), BT)


# In[15]:


#Rotation matrix
W= np.array([[math.cos(alphaw)*math.sin(betaw), math.sin(alphaw)*math.sin(betaw), math.cos(betaw)], 
            [-math.sin(alphaw), math.cos(alphaw), 0],
            [-math.cos(alphaw)*math.cos(betaw), -math.sin(alphaw)*math.cos(betaw), math.sin(betaw)]])

WT= np.transpose(W)


# In[16]:


x= (sigmaecs2bcs[0][0]+sigmaecs2bcs[1][1])
y= (sigmaecs2bcs[0][0]-sigmaecs2bcs[1][1])


# In[20]:


pm1=float(input("Mud Weight in Psi: "))


# CHECKING FAILURE FOR PLANE OF WEAKNESS:

# In[21]:



theta=0
pm=pm1
lists1=[]
lists2=[]
lists3=[]
while  (theta<180):
    thet=theta*0.0175
    sigmarr= pm-pp
    sigma00= x-2*y*math.cos(2*thet)-4*sigmaecs2bcs[0][1]*math.sin(2*thet)-pm-pp
    sigmazz=sigmaecs2bcs[2][2]-2*mu*y*math.cos(2*thet)-4*mu*sigmaecs2bcs[0][1]*math.sin(2*thet)-pp
    taur0= 0
    taurz= 0
    tau0z= 2*(-sigmaecs2bcs[0][2]*math.sin(2*thet)+sigmaecs2bcs[1][2]*math.cos(2*thet))
    
    sigmaccs= [[sigmarr, taur0, taurz], [taur0, sigma00, tau0z], [taurz, tau0z, sigmazz]]

   
    
    eigenvalue= la.eigvals(sigmaccs)
    
    sigma1= np.amax(abs(eigenvalue))
    sigma3= np.amin(abs(eigenvalue))
    C= [[math.cos(thet), math.sin(thet), 0], [-math.sin(thet), math.cos(thet), 0], [0,0,1]]
    
    CT= np.transpose(C)
        
    i= (np.dot(np.dot(W,BT), CT))
    j= np.dot(np.dot(i, sigmaccs), C)
    k= np.dot(np.dot(j,B), WT)
    sigmaecs2wcs= k

    
    tauw= math.sqrt((np.power((k[0][1]),2)+np.power((k[0][2]),2)))
    sigmaw= k[0][0]
    v= tauw/(mw*sigmaw+sw)
    if v>=1:
        lists1.append(theta)
        lists2.append(180+theta)
        pm=pm+0.0052*d
        continue
       
       
        
    else:
        lists3.append(pm)
        pm=pm1
        theta= theta+1


        


# CHECKING FAILURE FOR ROCK MATRIX:
# 

# In[22]:


theta=0
pm=pm1
lists4=[]
lists5=[]
lists6=[]
while  (theta<180):

    thet=theta*0.0175
    
    sigmarr= pm-pp
    sigma00= x-2*y*math.cos(2*thet)-4*sigmaecs2bcs[0][1]*math.sin(2*thet)-pm-pp
    sigmazz=sigmaecs2bcs[2][2]-2*mu*y*math.cos(2*thet)-4*mu*sigmaecs2bcs[0][1]*math.sin(2*thet)-pp
    taur0= 0
    taurz= 0
    tau0z= 2*(-sigmaecs2bcs[0][2]*math.sin(2*thet)+sigmaecs2bcs[1][2]*math.cos(2*thet))
    
    sigmaccs= [[sigmarr, taur0, taurz], [taur0, sigma00, tau0z], [taurz, tau0z, sigmazz]]
   
    
    eigenvalue= la.eigvals(sigmaccs)
    
    sigma1= np.amax(abs(eigenvalue))
    sigma3= np.amin(abs(eigenvalue))
    
 
    u= (sigma1)/abs(((sigma3)+2*(si+mi*sigma3)*(np.sqrt(1+mi*mi) +mi)))
    
    if u>=1:
        lists4.append(theta)
        lists5.append(theta+180)
        pm=pm+0.0052*d
        continue
    
    else:
        lists6.append(pm)
        pm=pm1
        theta= theta+1


# PLOT FOR FAILURE PATTERN:

# In[27]:


print("Estimated Failure Patterns : Green: Rock Matrix & Red: Weak Plane")
r=50
fig=plt.figure()
ax=fig.gca()
ax=plt.subplot(111,projection='polar')
ax.set_theta_direction(-1)
ax.set_theta_offset(math.pi/2)
for j in range(len(lists4)):
    
    (plt.polar(lists4[j]*np.pi/180, r,'g.'))
    (plt.polar(lists5[j]*np.pi/180,r,'g.'))
plt.show()    

for i in range(len(lists1)):
    (plt.polar(lists1[i]*np.pi/180, r,'r.'))
    (plt.polar(lists2[i]*np.pi/180,r,'r.'))
plt.show()
    


# In[24]:


print("Lower Critical Mud Weight for Rock Matrix:")
print(max(lists6)/d/0.052)


# In[25]:


print("Lower Critical Mud Weight for Weak Plane:")
print(max(lists3)/d/0.052)


# In[ ]:





# In[ ]:





# In[ ]:




