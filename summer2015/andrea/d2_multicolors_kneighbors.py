#!/usr/bin/env python

import os
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import svm
from astroML.utils import completeness_contamination


#def kernel_svm(colors, labels):
 #   from sklearn.svm import SVC
 #   clf = SVC(kernel='rbf')
 #   clf.fit(colors, labels)
 #   pred = clf.predict(colors)
 #   compl, contam = completeness_contamination(pred,labels)
 #   return clf, compl, contam

def kneighbor(colors,labels):
    from sklearn.neighbors import KNeighborsClassifier
    kvals = 3 #number of neighbors the point looks at - makes run time longer
              #as the number decreases
    clf = KNeighborsClassifier(n_neighbors=kvals)
    clf.fit(colors, labels)
    pred = clf.predict(colors)
    compl, contam = completeness_contamination(pred,labels)
    return clf, compl, contam

def get_deep2_multicolors(usesim=False):
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
        filename = '/home/desi2/data'+'/deep2egs-oii.fits.gz'
        oii = fits.getdata(filename,1)
        rmagcut = (oii['CFHTLS_R']<23.4)*1
        oii = oii[:][np.where(rmagcut==1)]
        print(oii['CFHTLS_R'].max())
        rz = oii['CFHTLS_R']-oii['CFHTLS_Z']
        gr = oii['CFHTLS_G']-oii['CFHTLS_R']
        zw1 = oii['CFHTLS_Z']-oii['W1']
        w1w2 = oii['W1']-oii['W2']
        zcut = (oii['z']>0.6)*1
        oiicut = (oii['oii_3727']>8E-17)*1
        rzg = np.vstack((rz,gr)).T
        rzgw1 = np.vstack((rz,gr,zw1)).T
        rzgw1w2 = np.vstack((rz,gr,zw1,w1w2)).T
        objtype = zcut & oiicut

        rzlim = [-0.5,2.0]
        grlim = [-0.5,1.5]
        w1lim = [-0.5,2.0]
        w2lim = [-0.5,2.0]

    return rzg, rzgw1, rzgw1w2, objtype, rzlim, grlim, w1lim, w2lim

def main():
    """Document me"""

    # Calls the above functions
    rzg, rzgw1, rzgw1w2, objtype, rzlim, grlim, w1lim, w2lim = get_deep2_multicolors(usesim=False)

    # K Neighbors
    print('Working on gr vs. rz...')
    knc_clf1, knc_compl1, knc_contam1 = kneighbor(rzg,objtype)
    print('Completeness = {}, contamination = {}'.format(knc_compl1,knc_contam1))

     # K Neighbors
    print('Working on gr vs. rz vs. w1...')
    knc_clf2, knc_compl2, knc_contam2 = kneighbor(rzgw1,objtype)
    print('Completeness = {}, contamination = {}'.format(knc_compl2,knc_contam2))

 # K Neighbors
    print('Working on gr vs. rz vs. w1 vs. w2...')
    knc_clf3, knc_compl3, knc_contam3 = kneighbor(rzgw1w2,objtype)
    print('Completeness = {}, contamination = {}'.format(knc_compl3,knc_contam3))

 # Kernel SVM
   # print('Working on gr vs. rz...')
   # svm_clf1, svm_compl1, svm_contam1 = kernel_svm(rzg,objtype)
   # print('Completeness = {}, contamination = {}'.format(svm_compl1,svm_contam1))

 # Kernel SVM
   # print('Working on gr vs. rz vs. w1...')
   # svm_clf2, svm_compl2, svm_contam2 = kernel_svm(rzgw1,objtype)
   # print('Completeness = {}, contamination = {}'.format(svm_compl2,svm_contam2))

  # Kernel SVM
   # print('Working on gr vs. rz vs. w1 vs. w2...')
   # svm_clf3, svm_compl3, svm_contam3 = kernel_svm(rzgw1w2,objtype)
   # print('Completeness = {}, contamination = {}'.format(svm_compl3,svm_contam3))
          
    #Plotting the Completeness and the contamination
    classifier = [0,1,2]
    compl = [knc_compl1,knc_compl2,knc_compl3]
    contam = [knc_contam1,knc_contam2,knc_contam3]
    markers = ['d','o','v']
    fig, ax = plt.subplots(1,1)
    for ii,mm in enumerate(markers):
        ax.plot(ii,compl[ii],'b'+mm,label='Completeness',markersize=15)
    ax.set_xlabel('classifier')
    ax.set_ylabel(r'Completeness')
    ax.set_xlim(-0.25,2.25)
    #ax.set_ylim(0.2,0.8)
    ax.set_ylim(0,1)
    ax.legend(loc='best', fancybox=True, framealpha=0.5)
    myxticks = (['gr-rz','gr-rz-w1','gr-rz-w1-w2'])
    plt.xticks(classifier,myxticks,rotation=45)
    ax2 = ax.twinx()
    for ii,mm in enumerate(markers):
        ax2.plot(ii,contam[ii],'r'+mm,label='Contamination',markersize=15)
    ax2.set_ylabel(r'Contamination')
    ax2.set_xlim(-0.25,2.25)
    #ax2.set_ylim(0.2,0.8)
    ax2.set_ylim(0,1)
    ax2.margins(0.2)


    plt.savefig(os.getenv('HOME')+'/classifiers/deep2_multicolors_kneighbors.png')
    plt.show()
 
if __name__ == '__main__':
    main()
    
