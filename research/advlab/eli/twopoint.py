import numpy as np
from operator import add
import matplotlib.pylab as plt
import math
from matplotlib.colors import LogNorm
import matplotlib as mpl
DD=np.loadtxt('onedddcalctest.dat',dtype='float')
DR=np.loadtxt('oneddrcalctest.dat',dtype='float')
RR=np.loadtxt('onedrrcalctest.dat',dtype='float')


DD = DD.transpose()
RR = RR.transpose()
DR = DR.transpose()

#DD+=np.flipud(DD)
#DR+=np.flipud(DR)
#RR+=np.flipud(RR)

#DDnew1=np.rot90(DD)
#DRnew1=np.rot90(DR)
#RRnew1=np.rot90(RR)
#DDnew2=np.rot90(DDnew1)
#DRnew2=np.rot90(DRnew1)
#RRnew2=np.rot90(RRnew1)

#DD+=DDnew2
#DR+=DRnew2
#RR+=RRnew2


'''
########## Bin Reduction ##########
def binred(array):
    import numpy as np
    binsnow=200
    binsneeded=100
    global newbins1
    global newbins
    newbins1=[]
    newbins=[]
    for i in range(0,200,2):
        #Adding columns
        global val1
        val1=(array[:,i]+array[:,i+1])
        newbins1.append(val1)
    newbins1=np.array(newbins1)
    return newbins1
    for j in range(0,200,2):
        global val
        val=(newbins1[j,:]+newbins1[j+1,:])
        newbins.append(val)
        return newbins
    newbins=np.array(newbins)    
                    
binred(DD)      
#print newbins        
#print val    
''' 




ndata=100000
nrand=100000

#print DD.shape

#print sum(sum(DD))
#print sum(sum(DR))
#print sum(sum(RR))
'''
# Rebin
DDnew = np.zeros((100,100))
DRnew = np.zeros((100,100))
RRnew = np.zeros((100,100))


for i in range(0,200,2):
    for j in range(0,200,2):
        DDnew[i/2][j/2] = DD[i][j] + DD[i+1][j] + DD[i][j+1] + DD[i+1][j+1]
        DRnew[i/2][j/2] = DR[i][j] + DR[i+1][j] + DR[i][j+1] + DR[i+1][j+1]
        RRnew[i/2][j/2] = RR[i][j] + RR[i+1][j] + RR[i][j+1] + RR[i+1][j+1]
DD = DDnew
DR = DRnew
RR = RRnew

'''

DD /=(ndata**2-ndata)/2.
DR /=(nrand*ndata)/1.
RR /=(nrand**2-nrand)/2.
theta = (DD - 2*DR + RR)/RR

#R^2 WEIGHTING

nbins=400
rangeval=300

# Correct for little h
rangeval *= 0.7

#R Values
'''
for i in range(nbins):
    for j in range(nbins):
        r2=((nbins/2)-i)**2 + (j-(nbins/2))**2
        theta[i][j] *= r2
'''

plt.figure(figsize=(8,8))


#extent=
#plot=plt.imshow(theta)
extent=[-rangeval,rangeval,-rangeval,rangeval]

plt.subplot(2,2,1)
a=plt.imshow(DD,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('DD')

plt.subplot(2,2,2)
b=plt.imshow(RR,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('RR')

plt.subplot(2,2,3)
c=plt.imshow(DR,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('DR')

### Mirror Over the X-Axis #### 
#nbins = 100
newtheta= np.zeros((nbins,nbins))
newtheta += theta
for i in range(0,nbins):
    newtheta[i] += theta[(nbins-1)-i]

plt.subplot(2,2,4)
d=plt.imshow(newtheta,extent=extent,norm=mpl.colors.LogNorm(vmin=0.03,vmax=0.2))
plt.colorbar(d)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title(r'$\xi$')

#plt.ylim(0,1)
#plt.xlim(0,30)
#plt.xlabel('Distance (Mpc)')
#plt.ylabel('Theta')
#plt.title('DR10 Correlation Estimator with 20000 Galaxies')

plt.tight_layout()

plt.show()

if 1:
    w1D = theta.transpose()[200]
    w1D2= w1D[200:]
    xvals = np.linspace(0,200,nbins/2)
    print w1D,xvals
    plt.figure()
    plt.plot(xvals,w1D2,'o')
    #plt.xlim(0,200)
    #plt.ylim(0,0.04)
    plt.xlabel('r')
    plt.ylabel('Xi')
    plt.title('1D Correlation Function (Grid)')
plt.show()

#np.savetxt('ourxi.dat',newtheta)

