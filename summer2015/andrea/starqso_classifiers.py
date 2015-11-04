#!/usr/bin/env python

import os
import pandas
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import svm
from astroML.utils import completeness_contamination

def plot_boundary(clf,xlim,ylim,ax,useproba=False):
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
    ax.contour(xx,yy,zz,[0.5],linewidths=2.5,colors='orange')

def naive_bayes(colors,labels):
    from sklearn.naive_bayes import GaussianNB
    clf = GaussianNB()
    clf.fit(colors,labels)
    pred = clf.predict(colors)
    compl, contam = completeness_contamination(pred,labels)
    return clf, compl, contam

def gaussmix_bayes(colors,labels):
    from astroML.classification import GMMBayes
    ncomp = 5 # number of completions, can over compensate
    clf = GMMBayes(ncomp)#,min_covar=1E-5,covariance_type='full')
    clf.fit(colors,labels)
    pred = clf.predict(colors)
    compl, contam = completeness_contamination(pred,labels)
    return clf, compl, contam

def kneighbor(colors,labels):
    from sklearn.neighbors import KNeighborsClassifier
    kvals = 10 #comment me
    clf = KNeighborsClassifier(n_neighbors=kvals)
    clf.fit(colors, labels)
    pred = clf.predict(colors)
    compl, contam = completeness_contamination(pred,labels)
    return clf, compl, contam

def kernel_svm(colors, labels):
    from sklearn.svm import SVC
    clf = SVC(kernel='rbf')
    clf.fit(colors, labels)
    pred = clf.predict(colors)
    compl, contam = completeness_contamination(pred,labels)
    return clf, compl, contam

def get_starqso(usesim=False):
    from astropy.io import fits
    """Document me"""

    if usesim:
        # Simulate some data
        np.random.seed(0)
        mu1 = [1, 1]
        cov1 = 0.3 * np.eye(2)
        mu2 = [5, 3]
        cov2 = np.eye(2) * np.array([0.4, 0.1])
        ugr = np.concatenate([np.random.multivariate_normal(mu1, cov1, 100),
                              np.random.multivariate_normal(mu2, cov2, 100)])
        objtype = np.zeros(200)
        objtype[100:] = 1
        uglim = [-1,5]
        grlim = [-1,8]
    else: 
        filename = os.getenv('HOME')+'/data/stars_qsos_sdss.fits'
        fitsdata = fits.getdata(filename,1)
        ug = fitsdata['u'] - fitsdata['g']
        gr = fitsdata['g'] - fitsdata['r']
        objtype = (fitsdata['class'] == 'QSO')*1
            
        ugr = np.vstack((gr,ug)).T
        uglim = [-0.5,3.0]
        grlim = [-0.5,1.5]

    return ugr, objtype, uglim, grlim

def main():
    """Document me"""

    # Comment
    ugr, objtype, uglim, grlim = get_starqso(usesim=False)

    # Gaussian NB
    print('Working on Gaussian NB classifier...')
    gnb_clf, gnb_compl, gnb_contam = naive_bayes(ugr,objtype)
    print('Completeness = {}, contamination = {}'.format(gnb_compl,gnb_contam))

    # Gaussian Mixture
    print('Working on Gaussian mixture Gauss classifier...')
    gmm_clf, gmm_compl, gmm_contam = gaussmix_bayes(ugr,objtype)
    print('Completeness = {}, contamination = {}'.format(gmm_compl,gmm_contam))

    # K Neighbors
    print('Working on Kneighbors classifier...')
    knc_clf, knc_compl, knc_contam = kneighbor(ugr,objtype)
    print('Completeness = {}, contamination = {}'.format(knc_compl,knc_contam))

    # Kernel SVM
    print('Working on Kernel SVM classifier...')
    svm_clf, svm_compl, svm_contam = kernel_svm(ugr,objtype)
    print('Completeness = {}, contamination = {}'.format(svm_compl,svm_contam))

    # Make the plot!
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
    plt.xlim(grlim)
    plt.ylim(uglim)
    #plt.xticks(fontsize=14)
    #plt.yticks(fontsize=14)
    ax3.set_xlabel('g-r')
    ax4.set_xlabel('g-r')
    ax3.set_ylabel('u-g')
    ax1.set_ylabel('u-g')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
  
    #removing the plot frame lines
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax4.spines['top'].set_visible(False)

    ax1.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax4.spines['right'].set_visible(False)


    #Axis ticks only on the bottom and left
    ax1.get_xaxis().tick_bottom()
    ax2.get_xaxis().tick_bottom()
    ax3.get_xaxis().tick_bottom()
    ax4.get_xaxis().tick_bottom()
    
    ax1.get_yaxis().tick_left() 
    ax2.get_yaxis().tick_left() 
    ax3.get_yaxis().tick_left() 
    ax4.get_yaxis().tick_left() 
    
    # Gaussian NB
    ax1.scatter(ugr[:,0],ugr[:,1],c=objtype)
    plot_boundary(gnb_clf,grlim,uglim,ax1)
    
    # Gaussian Mixed
    ax2.scatter(ugr[:,0],ugr[:,1],c=objtype)
    plot_boundary(gmm_clf,grlim,uglim,ax2,useproba=True)

    #KNeighbors
    ax3.scatter(ugr[:,0],ugr[:,1],c=objtype)
    plot_boundary(knc_clf,grlim,uglim,ax3)

    #Kernel SVM
    ax4.scatter(ugr[:,0],ugr[:,1],c=objtype)
    plot_boundary(svm_clf,grlim,uglim,ax4)

    plt.savefig(os.getenv('HOME')+'/classifiers/starqso_classifiers.png')

    #Plotting the Completeness and the contamination
    classifier = [0,1,2,3]
    compl = [gnb_compl,gmm_compl,knc_compl,svm_compl]
    contam = [gnb_contam,gmm_contam,knc_contam,svm_contam]

    fig, ax = plt.subplots(1,1)
    fig.subplots_adjust(bottom=0.25)
    ax.plot(compl,'o',label='Completeness',markersize=20)
    ax.set_xlabel('classifier')
    ax.set_ylabel(r'Completeness')
    ax.set_xlim(-0.25,3.25)
    ax.set_ylim(0,1)
   # ax.legend()#loc='upper left',frameon=True)
    myxticks = (['Niave Bayes Gaussian','Gaussian Mixed','K-Neighbors','Kernel SVM'])
    plt.xticks(classifier,myxticks,rotation=45)
    ax2 = ax.twinx()
    ax2.plot(contam,'rd',label='Contamination',markersize=20)
    ax2.set_ylabel(r'Contamination')
    ax2.set_xlim(-0.25,3.25)
    ax2.set_ylim(0,1)
    ax2.margins(0.2)

    
    plt.savefig(os.getenv('HOME')+'/classifiers/starqso_contam.png')
    plt.show()

if __name__ == '__main__':
    main()
    
