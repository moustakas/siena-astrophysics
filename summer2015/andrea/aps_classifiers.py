#!/usr/bin/env python

import os, sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import svm, tree
from astroML.utils import completeness_contamination

def plot_boundary(clf,xlim,ylim,ax,useproba=False):
    sns.set(palette='Reds',style='white')
    hh = 0.01  # step size in the mesh
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

def tree(colors, labels):
    from sklearn.tree import DecisionTreeClassifier as dtc
    from sklearn.cross_validation import StratifiedKFold
    clf = dtc(max_depth=12,criterion='entropy')
    nfold = 9
    skf = StratifiedKFold(labels, nfold)
    compl = np.zeros(nfold)
    contam = np.zeros(nfold)
    ii = 0
    #import pdb ; pdb.set_trace()
    for train, test in skf:
    #for ii, enumerate(zip(train, test)) in skf:
        clf.fit(colors[train,:], labels[train])
        pred = clf.predict(colors[test,:])
        compl1, contam1 = completeness_contamination(pred, labels[test])
        #print(compl1, contam1)
        compl[ii] = compl1
        contam[ii] = contam1
        ii += 1
    #import pdb ; pdb.set_trace()
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
    import pandas
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
        filename = '/home/desi2/aps'+'/deep2-field1-oii.fits.gz'
        #filestars = '/home/desi2/data'+'/deep2-stars.fits.gz'
        oii = fits.getdata(filename,1)
        #stars = fits.getdata(filestars,1)

        minmag = 23.4
        oii = oii[np.where((oii['CFHTLS_R']<minmag)*1)]
        #stars = stars[np.where((stars['CFHTLS_R']<minmag)*1)]

        ngal = len(oii)
        #nstar = len(stars)

        gr = oii['CFHTLS_G']-oii['CFHTLS_R']
        rz = oii['CFHTLS_R']-oii['CFHTLS_Z']
        oiiflux = oii['OII_3727']
        redshift = oii['ZHELIO']
        field = oii['FIELD']
        igood = ((redshift>1.0)*1)*((oiiflux>8E-17)*1)
        colorgood = ((gr>-0.3)*1)*((gr<1.8)*1)*((rz>-0.3)*1)*((rz<2.1)*1)
        istar = np.zeros(ngal)
        
        data = np.vstack((gr,rz,oiiflux,redshift,field,istar,igood,colorgood)).T
        pdata = pandas.DataFrame(data,columns=['g-r','r-z','oiiflux','redshift',
                                               'field','istar','igood','colorgood'])

        #gr = np.concatenate((oii['CFHTLS_G']-oii['CFHTLS_R'],stars['CFHTLS_G']-stars['CFHTLS_R']))
        #rz = np.concatenate((oii['CFHTLS_R']-oii['CFHTLS_Z'],stars['CFHTLS_R']-stars['CFHTLS_Z']))
        #oiiflux = np.concatenate((oii['OII_3727'],np.zeros(nstar)))
        #redshift = np.concatenate((oii['ZBEST'],np.zeros(nstar)))
        #field = np.concatenate((oii['FIELD'],stars['FIELD']))
        #istar = np.concatenate((np.zeros(ngal),np.ones(nstar)))
        #igood = ((redshift>0.8)*1)*((oiiflux>8E-17)*1)
        #colorgood = ((gr>-0.3)*1)*((gr<1.8)*1)*((rz>-0.3)*1)*((rz<2.1)*1)
        #
        #data = np.vstack((gr,rz,oiiflux,redshift,field,istar,igood,colorgood)).T
        #pdata = pandas.DataFrame(data,columns=['g-r','r-z','oiiflux','redshift',
        #                                       'field','istar','igood','colorgood'])

        rzlim = [-0.3,2.1]
        grlim = [-0.3,1.8]

    return pdata, rzlim, grlim

def main():
    """Document me"""

    # Comment
    pdata, rzlim, grlim = get_deep2(usesim=False)

    f1 = pdata.loc[(pdata['field']==1)&(pdata['colorgood']==1)]
    rzg = np.array([f1['r-z'].values,f1['g-r'].values]).T
    igood = f1['igood'].values

    # Decision Tree
    tree_clf, tree_compl, tree_contam = tree(rzg, igood)
    print(tree_compl, tree_contam)

    galgood = pdata.loc[(pdata['field']==1)&(pdata['istar']==0)&
                        (pdata['igood']==1)&(pdata['colorgood']==1)]
    galbad = pdata.loc[(pdata['field']==1)&(pdata['istar']==0)&
                       (pdata['igood']==0)&(pdata['colorgood']==1)]

    fig, ax = plt.subplots(1,1)
    ax.scatter(galbad['r-z'].values,galbad['g-r'].values,c='blue',marker='o',s=4)
    ax.scatter(galgood['r-z'].values,galgood['g-r'].values,c='magenta',marker='d',s=20)
    #sns.kdeplot(stars['r-z'].values,stars['g-r'].values,ax=ax2,gridsize=50)
    plot_boundary(tree_clf,rzlim,grlim,ax,useproba=True)
    ax.set_xlim(rzlim)
    ax.set_ylim(grlim)
    plt.show()

    #import pdb ; pdb.set_trace()

    sys.exit(1)

    ## Gaussian NB
    #print('Working on Gaussian NB classifier...')
    #gnb_clf, gnb_compl, gnb_contam = naive_bayes(rzg,igood)
    #print('Completeness = {}, contamination = {}'.format(gnb_compl,gnb_contam))
    #
    ## Gaussian Mixture
    #print('Working on Gaussian mixture Gauss classifier...')
    #gmm_clf, gmm_compl, gmm_contam = gaussmix_bayes(rzg,igood)
    #print('Completeness = {}, contamination = {}'.format(gmm_compl,gmm_contam))
    #
    ## K Neighbors
    #print('Working on Kneighbors classifier...')
    #knc_clf, knc_compl, knc_contam = kneighbor(rzg,igood)
    #print('Completeness = {}, contamination = {}'.format(knc_compl,knc_contam))
    #
    ## Kernel SVM
    #print('Working on Kernel SVM classifier...')
    #svm_clf, svm_compl, svm_contam = kernel_svm(rzg,igood)
    #print('Completeness = {}, contamination = {}'.format(svm_compl,svm_contam))

    # Make the plot!
    galgood = pdata.loc[(pdata['field']==1)&(pdata['istar']==0)&
                        (pdata['igood']==1)&(pdata['colorgood']==1)]
    galbad = pdata.loc[(pdata['field']==1)&(pdata['istar']==0)&
                       (pdata['igood']==0)&(pdata['colorgood']==1)]
    
    gal = pdata.loc[(pdata['istar']==0)]
    stars = pdata.loc[(pdata['istar']==1)&(pdata['colorgood']==1)]

    pal = sns.color_palette("Greys_r")
    
    sns.set_style(style='white')
    fig, ((ax1, ax2), (ax3, ax4))= plt.subplots(2,2,sharex=True,sharey=True)
    
    plt.xlim(rzlim)
    plt.ylim(grlim)
    ax1.tick_params(labelsize=20)
    ax2.tick_params(labelsize=20)
    ax3.tick_params(labelsize=20)
    ax4.tick_params(labelsize=20)

    ax1.spines['bottom'].set_linewidth(3)
    ax1.spines['left'].set_linewidth(3)
    ax1.spines['top'].set_linewidth(3)
    ax1.spines['right'].set_linewidth(3)
    ax2.spines['bottom'].set_linewidth(3)
    ax2.spines['left'].set_linewidth(3)
    ax2.spines['top'].set_linewidth(3)
    ax2.spines['right'].set_linewidth(3)
    ax3.spines['bottom'].set_linewidth(3)
    ax3.spines['left'].set_linewidth(3)
    ax3.spines['top'].set_linewidth(3)
    ax3.spines['right'].set_linewidth(3)
    ax4.spines['bottom'].set_linewidth(3)
    ax4.spines['left'].set_linewidth(3)
    ax4.spines['top'].set_linewidth(3)
    ax4.spines['right'].set_linewidth(3)
    
    plt.tight_layout(pad=0.9, w_pad=0.5, h_pad=1.0)
    ax3.set_xlabel('r-z',fontsize=20)
    ax4.set_xlabel('r-z',fontsize=20)
    ax3.set_ylabel('g-r',fontsize=20)
    ax1.set_ylabel('g-r',fontsize=20)

    
   
    # Gaussian NB
    ax1.scatter(galbad['r-z'].values,galbad['g-r'].values,c='blue',marker='o',s=4)
    ax1.scatter(galgood['r-z'].values,galgood['g-r'].values,c='magenta',marker='d',s=20)
    sns.kdeplot(stars['r-z'].values,stars['g-r'].values,ax=ax1,gridsize=50)
    plot_boundary(gnb_clf,rzlim,grlim,ax1)

    
    # Gaussian Mixed
    ax2.scatter(galbad['r-z'].values,galbad['g-r'].values,c='blue',marker='o',s=4)
    ax2.scatter(galgood['r-z'].values,galgood['g-r'].values,c='magenta',marker='d',s=20)
    sns.kdeplot(stars['r-z'].values,stars['g-r'].values,ax=ax2,gridsize=50)
    plot_boundary(gmm_clf,rzlim,grlim,ax2,useproba=True)

 
    #KNeighbors
    ax3.scatter(galbad['r-z'].values,galbad['g-r'].values,c='blue',marker='o',s=4)
    ax3.scatter(galgood['r-z'].values,galgood['g-r'].values,c='magenta',marker='d',s=40)
    sns.kdeplot(stars['r-z'].values,stars['g-r'].values,ax=ax3,gridsize=50)
    plot_boundary(knc_clf,rzlim,grlim,ax3)

    
    #Kernel SVM
    ax4.scatter(galbad['r-z'].values,galbad['g-r'].values,c='blue',marker='o',s=4)
    ax4.scatter(galgood['r-z'].values,galgood['g-r'].values,c='magenta',marker='d',s=20)
    sns.kdeplot(stars['r-z'].values,stars['g-r'].values,ax=ax4,gridsize=50)
    plot_boundary(svm_clf,rzlim,grlim,ax4)


    fig.subplots_adjust(bottom=0.15,left=0.12)
    plt.savefig(os.getenv('HOME')+'/classifiers/deep2_classifiers.png')
     
    #Plotting the Completeness and the contamination
    classifier = [0,1,2,3]
    compl = [gnb_compl,gmm_compl,knc_compl,svm_compl]
    contam = [gnb_contam,gmm_contam,knc_contam,svm_contam]
    markers = ['x','o','v','d']

    fig, ax = plt.subplots(1,1)
    fig.subplots_adjust(bottom=0.24,left=0.12)

    ax.spines['bottom'].set_linewidth(3)
    ax.spines['left'].set_linewidth(3)
    ax.spines['top'].set_linewidth(3)
    ax.spines['right'].set_linewidth(3)
    ax.plot(compl,'bo',markersize=20,label='Completeness')
    #ax.set_xlabel('Classifier',fontsize=20)
    ax.set_ylabel(r'Completeness',fontsize=20)
    ax.set_xlim(-0.25,3.25)
    ax.set_ylim(0,1)
    ax.legend(loc='upper left',frameon=True)
    myxticks = (['Niave Bayes','Gaussian Mixed','K-Neighbors','Kernel SVM'])
    plt.xticks(classifier,myxticks,rotation=20)
    ax2 = ax.twinx()
    ax2.spines['bottom'].set_linewidth(3)
    ax2.spines['left'].set_linewidth(3)
    ax2.spines['top'].set_linewidth(3)
    ax2.spines['right'].set_linewidth(3)
    ax.tick_params(labelsize=20)
    ax2.tick_params(labelsize=20)
    ax2.plot(contam,'rd',markersize=20,label='Contamination')
    ax2.set_ylabel(r'Contamination',fontsize=20)
    ax2.set_xlim(-0.25,3.25)
    ax2.set_ylim(0,1)
    ax2.margins(0.2)
    ax2.legend(frameon=True)


    plt.savefig(os.getenv('HOME')+'/classifiers/deep2_contam.png') 
    plt.show()
 
if __name__ == '__main__':
    main()
    
