#!/usr/bin/env python

import os
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
    ax.contour(xx,yy,zz,[0.5],colors = 'y',linewidths=2)

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
    kvals = 3 #comment me
    clf = KNeighborsClassifier(n_neighbors=kvals,p=1)
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

def get_deep2(usesim=False):
    from astropy.io import fits
    """Document me"""

    if usesim:
        # Simulate some data
        np.random.seed(0)
        mu1 = [1, 1]
        cov1 = 0.3 * np.eye(2)
        mu2 = [5, 3]
        cov2 = np.eye(2) * np.array([0.4, 0.1])
        rzg = np.concatenate([np.random.multivariate_normal(mu1, cov1, 100),
                              np.random.multivariate_normal(mu2, cov2, 100)])
        objtype = np.zeros(200)
        objtype[100:] = 1
        rzlim = [-0.5,2.0]
        grlim = [-0.5,1.5]
    else: 
        filename = '/home/desi2/data'+'/deep2-oii.fits.gz'
        filestars = '/home/desi2/data'+'/deep2-stars.fits.gz'
        oii = fits.getdata(filename,1)
        stars = fits.getdata(filestars,1)
        rmagcut = (oii['CFHTLS_R']<23.4)*1
        oii = oii[:][np.where(rmagcut==1)]
        #print(oii['CFHTLS_R'].max())
        rz = oii['CFHTLS_R']-oii['CFHTLS_Z']
        rz_stars = stars['CFHTLS_R']-stars['CFHTLS_Z']
        #rz = 10**(-0.4*(oii['CFHTLS_R']-oii['CFHTLS_Z'])) # flux ratio
        gr = oii['CFHTLS_G']-oii['CFHTLS_R']
        gr_stars = stars['CFHTLS_G']-stars['CFHTLS_R']
        #gr = 10**(-0.4*(oii['CFHTLS_G']-oii['CFHTLS_R']))
        zcut = (oii['z']>0.6)*1
        oiicut = (oii['oii_3727']>8E-17)*1
        rzg_oii = np.vstack((rz,gr)).T
        rzg_stars = np.vstack((rz_stars,gr_stars)).T
        objtype = zcut & oiicut
        objtype_stars = np.zeros(len(stars))

        rzlim = [-0.5,2.0]
        grlim = [-0.5,1.5]

    return rzg_oii, objtype, rzlim, grlim, rzg_stars, objtype_stars

def main():
    """Document me"""

    # Comment
    rzg_oii, objtype, rzlim, grlim, rzg_stars, objtype_stars = get_deep2(usesim=False)

    # Gaussian NB
    print('Working on Gaussian NB classifier...')
    gnb_clf, gnb_compl, gnb_contam = naive_bayes(rzg_oii,objtype)
    print('Completeness = {}, contamination = {}'.format(gnb_compl,gnb_contam))

    # Gaussian Mixture
    print('Working on Gaussian mixture Gauss classifier...')
    gmm_clf, gmm_compl, gmm_contam = gaussmix_bayes(rzg_oii,objtype)
    print('Completeness = {}, contamination = {}'.format(gmm_compl,gmm_contam))

    # K Neighbors
    print('Working on Kneighbors classifier...')
    knc_clf, knc_compl, knc_contam = kneighbor(rzg_oii,objtype)
    print('Completeness = {}, contamination = {}'.format(knc_compl,knc_contam))

    # Kernel SVM
    print('Working on Kernel SVM classifier...')
    svm_clf, svm_compl, svm_contam = kernel_svm(rzg_oii,objtype)
    print('Completeness = {}, contamination = {}'.format(svm_compl,svm_contam))

    # Make the plot!

   # sns.set(palette='Reds',style='ticks')
    pal = sns.color_palette("Blues_r")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
    plt.xlim(rzlim)
    plt.ylim(grlim)
    plt.tight_layout(pad=0.9, w_pad=0.5, h_pad=1.0)
    ax3.set_xlabel('r-z')
    ax4.set_xlabel('r-z')
    ax3.set_ylabel('g-r')
    ax1.set_ylabel('g-r')
   
    # Gaussian NB
    ax1.scatter(rzg_oii[:,0],rzg_oii[:,1],c=objtype)
    plot_boundary(gnb_clf,rzlim,grlim,ax1)
    sns.kdeplot(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,ax=ax1,color=pal)
    
    # Gaussian Mixed
    ax2.scatter(rzg_oii[:,0],rzg_oii[:,1],c=objtype)
    plot_boundary(gmm_clf,rzlim,grlim,ax2,useproba=True)
    #ax2.scatter(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,color='g',s=50)
    sns.kdeplot(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,ax=ax2,bw="silverman",color=pal)


    #KNeighbors
    ax3.scatter(rzg_oii[:,0],rzg_oii[:,1],c=objtype)
    plot_boundary(knc_clf,rzlim,grlim,ax3)
    #ax3.scatter(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,color='g',s=50)
    sns.kdeplot(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,ax=ax3, color=pal)

    
    #Kernel SVM
    ax4.scatter(rzg_oii[:,0],rzg_oii[:,1],c=objtype)
    plot_boundary(svm_clf,rzlim,grlim,ax4)
    #ax4.scatter(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,color='g',s=50)
    sns.kdeplot(rzg_stars[:,0],rzg_stars[:,1],c=objtype_stars,ax=ax4, color=pal)

     plt.savefig(os.getenv('HOME')+'/classifiers/deep2_classifiers.png')

     
    #Plotting the Completeness and the contamination
    classifier = [0,1,2,3]
    compl = [gnb_compl,gmm_compl,knc_compl,svm_compl]
    contam = [gnb_contam,gmm_contam,knc_contam,svm_contam]
    markers = ['x','o','v','d']

    fig, ax = plt.subplots(1,1)
    ax.plot(compl,'bo',markersize=15,label='Completeness')
    ax.set_xlabel('classifier')
    ax.set_ylabel(r'Completeness')
    ax.set_xlim(-0.25,3.25)
    ax.set_ylim(0,1)
    ax.legend(loc='upper left',frameon=True)
    myxticks = (['Niave Bayes Gaussian','Gaussian Mixed','K-Neighbors','Kernel SVM'])
    plt.xticks(classifier,myxticks,rotation=45)
    ax2 = ax.twinx()
    ax2.plot(contam,'rd',markersize=15,label='Contamination')
    ax2.set_ylabel(r'Contamination')
    ax2.set_xlim(-0.25,3.25)
    ax2.set_ylim(0,1)
    ax2.margins(0.2)
    ax2.legend(frameon=True)


    plt.savefig(os.getenv('HOME')+'/classifiers/deep2_contam.png') 
    plt.show()
 
if __name__ == '__main__':
    main()
    
