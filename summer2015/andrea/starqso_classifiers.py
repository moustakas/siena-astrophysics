#!/usr/bin/env python

import os
import pandas
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import svm
from astroML.utils import completeness_contamination

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
            
        data = np.vstack((gr,ug,objtype)).T
        pdata = pandas.DataFrame(data,columns=['g-r','u-g','type'])
        uglim = [-0.5,3.0]
        grlim = [-0.5,1.5]

    return pdata, uglim, grlim
    #return ugr, objtype, uglim, grlim

def main():
    """Document me"""

    # Comment
    data, uglim, grlim = get_starqso(usesim=False)
    #ugr, objtype, uglim, grlim = get_starqso(usesim=False)

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

    # Temporary hack!
    sns.set(context='talk',style='dark',font_scale=1.5)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
    sns.lmplot('g-r','u-g',data,hue='type',markers=['s','o'],
               fit_reg=False,ax=ax1)
    sns.lmplot('g-r','u-g',data,hue='type',markers=['s','o'],
               fit_reg=False,ax=ax2)
    sns.lmplot('g-r','u-g',data,hue='type',markers=['s','o'],
               fit_reg=False,ax=ax3)
    sns.lmplot('g-r','u-g',data,hue='type',markers=['s','o'],
               fit_reg=False,ax=ax4)
    fig.tight_layout(pad=0.9, w_pad=0.5, h_pad=1.0)
    plt.show()

    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
    plt.xlim(grlim)
    plt.ylim(uglim)
    ax1.scatter(data['g-r'].values,data['u-g'].values,c=data['type'].values,
                cmap='Reds')
    plt.show()
    


    plt.xlim(grlim)
    plt.ylim(uglim)

    ax3.set_xlabel('g-r')
    ax4.set_xlabel('g-r')
    ax3.set_ylabel('u-g')
    ax1.set_ylabel('u-g')
    sns.lmplot('g-r','u-g',data,hue='type',markers=['s','o'],
               fit_reg=False,ax=ax1,sharex=True,sharey=True)
    plt.show()



    ax2.scatter(ugr[:,0],ugr[:,1],c=objtype)
    ax3.scatter(ugr[:,0],ugr[:,1],c=objtype)
    ax4.scatter(ugr[:,0],ugr[:,1],c=objtype)
    plt.show()
    
    
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
    plt.xlim(grlim)
    plt.ylim(uglim)
    plt.tight_layout(pad=0.9, w_pad=0.5, h_pad=1.0)
    ax3.set_xlabel('g-r')
    ax4.set_xlabel('g-r')
    ax3.set_ylabel('u-g')
    ax1.set_ylabel('u-g')
   
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

    #Plotting the Completeness and the contamination
    classifier = [0,1,2,3]
    compl = [gnb_compl,gmm_compl,knc_compl,svm_compl]
    contam = [gnb_contam,gmm_contam,knc_contam,svm_contam]
    markers = ['x','o','v','d']

    fig, ax = plt.subplots(1,1)
    ax.plot(compl,'bo',label='Completeness')
    ax.set_xlabel('classifier')
    ax.set_ylabel(r'Completeness')
    ax.set_xlim(-0.25,3.25)
    ax.set_ylim(0,1)
    ax.legend(loc='upper left',frameon=True)
    myxticks = (['Niave Bayes Gaussian','Gaussian Mixed','K-Neighbors','Kernel SVM'])
    plt.xticks(classifier,myxticks,rotation=45)
    ax2 = ax.twinx()
    ax2.plot(contam,'rd',label='Contamination')
    ax2.set_ylabel(r'Contamination')
    ax2.set_xlim(-0.25,3.25)
    ax2.set_ylim(0,1)
    ax2.margins(0.2)
    ax2.legend(frameon=True)


    plt.savefig(os.getenv('HOME')+'/classifiers/starqso_classifiers.png')
    plt.show()

if __name__ == '__main__':
    main()
    
