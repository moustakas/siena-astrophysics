#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astroML.utils import completeness_contamination

#reading in OII fitsfile and the values for rz and gr
deep2dir = '/home/desi2/data/'
oii = fits.getdata(deep2dir + 'deep2egs-oii.fits.gz',1)
rz_oii = oii['CFHTLS_R']-oii['CFHTLS_Z']
gr_oii = oii['CFHTLS_G']-oii['CFHTLS_R']

#finding the galaxies that have a redshift higher than 1 and a concentration of
#OII larger than stated value
zcut = (oii['z']>1.0)*1
oiicut = (oii['oii_3727']>8E-17)*1

#creating variables for the parameters needed for plt.scatter
colors_oii = np.vstack((rz_oii,gr_oii)).T
labels_oii = zcut & oiicut



#stars = fits.getdata('deep2egs-stars.fits.gz',1)
#rz_stars = stars['CFHTLS_R']-stars['CFHTLS_Z']
#gr_stars = stars['CFHTLS_G']-stars['CFHTLS_R']

#gnb = GaussianNB()
#dofit = gnb.fit(colors_oii,labels_oii)
hh = 0.02

def plot_boundary(clf,xlim,ylim,ax,useproba=False):
    sns.set(palette='Reds',style='white')
    hh = 0.03  # step size in the mesh
    xx, yy = np.meshgrid(np.arange(xlim[0],xlim[1],hh),
                         np.arange(ylim[0],ylim[1],hh))
    if useproba:
        zz = clf.predict_proba(np.c_[xx.ravel(),yy.ravel()])
        zz = zz[:,1]
    else:
        zz = clf.predict(np.c_[xx.ravel(),yy.ravel()])
    zz = zz.reshape(xx.shape)
    #ax.pcolormesh(xx,yy,zz,cmap=plt.cm.Paired)
    ax.contour(xx,yy,zz,[0.5],linewidths=5)


def kneighbor(colors_oii,labels_oii):
    from sklearn.neighbors import KNeighborsClassifier
    kvals = 10 #comment me
    clf = KNeighborsClassifier(n_neighbors=kvals)
    clf.fit(colors_oii, labels_oii)
    pred = clf.predict(colors_oii)
    compl, contam = completeness_contamination(pred,labels_oii)
    return clf, compl, contam

# K Neighbors
print('Working on Kneighbors classifier...')
knc_clf, knc_compl, knc_contam = kneighbor(colors_oii,labels_oii)
print('Completeness = {}, contamination = {}'.format(knc_compl,knc_contam))

#setting the min and max values
minoii = 8E-17 # erg/s/cm2
rzmin = -0.5
rzmax = 2.0
grmin = -0.5
grmax = 1.5
hh=0.2
#xx,yy = np.meshgrid(np.arange(rzmin,rzmax,hh),np.arange(grmin,grmax,hh))
#zz = clf.predict(np.c_[xx.ravel(),yy.ravel()])
#zz = zz.reshape(xx.shape)



fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.contour(xx, yy, zz, [0.5], lw = 2,colors='g')
#ax.scatter(colors_oii[:,0],colors_oii[:,1],c=labels_oii,marker='o', alpha=0.6)
#plt.axis('tight')
#plt.xlim([rzmin,rzmax])
#plt.ylim([grmin,grmax])

    #Plotting the Completeness and the contamination
classifier = [1]
compl = [knc_compl]
contam = [knc_contam]
markers = ['d']

fig, ax = plt.subplots(1,1)
ax.plot(compl,'bo',label='Completeness')
ax.set_xlabel('classifier')
ax.set_ylabel(r'Completeness')
ax.set_xlim(-0.25,3.25)
ax.set_ylim(0,1)
ax.legend(loc='upper left',frameon=True)
myxticks = (['K-Neighbors'])
plt.xticks(classifier,myxticks,rotation=45)
ax2 = ax.twinx()
ax2.plot(contam,'rd',label='Contamination')
ax2.set_ylabel(r'Contamination')
ax2.set_xlim(-0.25,3.25)
ax2.set_ylim(0,1)
ax2.margins(0.2)
ax2.legend(frameon=True)

#setting the min and max values
plt.savefig('/home/desi2/deep2_gauss_ex.png')

plt.show()
# 1. make a plot of g-r vs r-z using cfhtls_g, cfhtls_r, ... for both
# the galaxies and the stars - code the stars, and the strong
# [OII]-emitting galaxies differently
