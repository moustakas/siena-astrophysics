import numpy as np
import math
r4585=np.loadtxt('cmass_dr11_north_ir4585.v11.1.release.txt')
ra=r4585[:,0]
dec=r4585[:,1]
z=r4585[:,2]
peculiarv=r4585[:,3]
weight_cboss=r4585[:,4]
weight_cp=r4585[:,5]
weight_z=r4585[:,6]
veto=r4585[:,7]
#ztrue=data4560[:,8]
#flag=data4560[:,9]
'''
##### Weight Cut ######

print 'Cutting Weights'
wc=weight_cboss>0
wcnew=weight_cboss[wc]

wp=weight_cp>0
wpnew=weight_cp[wp]

wz=weight_z>0
wznew=weight_z[wz]

##### Veto Cut #######

print 'Veto Cut'
vc=veto==1
vcnew=veto[vc]

##### Final Cut #######

print 'Finalizing Random Cut'
tot=wc*wp*wz*vc
totrand=r4560[tot]

########### Distances ################

a=np.arange(0,len(totrand))
np.random.shuffle(a)
sample=totrand[a[0:5000]]
H=70
c=2.99792458*10**5
Zshift=sample[:,2] #Redshift Values
Vr=Zshift*c #Recession Velocity of Galaxy in km/s
rho=Vr/H  #Distance from Earth to Galaxy in Mpc

# Converting RA and Dec to Radians

RArad=(sample[:,0])*((math.pi)/180)
Decrad=((sample[:,1])*((math.pi)/180))

# Converting to Cartesian Coordinates

x=rho*np.sin(Decrad)*np.cos(RArad)
y=rho*np.sin(Decrad)*np.cos(RArad)
z=rho*np.cos(Decrad)

# Distances
print 'Calculating Distances'
d = []

for i in range(0,len(sample)):     #Cycling through each data point
    stepx=x[i]
    stepy=y[i]
    stepz=z[i]
    print i
    for n in range(i+1,len(sample)):   #Finding distances between a data point and every other data point
        val=[np.sqrt(((stepx-x[n])**2)+((stepy-y[n])**2)+((stepz-z[n])**2))]
        d.append(val[0])
d = np.array(d)

#####################################################

# Histogram

import matplotlib.pylab as plt
hist=plt.hist(d,bins=100)
plt.xlabel('Galactic Distances (Mpc)')
plt.ylabel('Frequency')
plt.title('Galactic Distance Distribution of 5000 Random Mock CMASS Galaxies')
plt.show()
frequ=hist[0]
dist=hist[1]
m=dist[0:-1]
diff=np.diff(dist)/2
mid=m+diff
vals=np.column_stack((mid,frequ))
'''

####### RA and Dec Scatter Plot ########

import matplotlib.pylab as plt
plt.plot(ra,dec,'o',markersize=0.2)
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Declination (Degrees)')
plt.title('Mock Data (North 4585)')
plt.xlim(90,280)
plt.ylim(-10,80)
plt.show()
