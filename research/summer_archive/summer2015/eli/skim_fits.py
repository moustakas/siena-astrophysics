import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import matplotlib.pylab as plt
#hdulist=fits.open('specObj-dr12.fits' , ignore_missing_end=True)

infilename = sys.argv[1]

tag = "default"

if len(sys.argv)>2:
        tag = sys.argv[2]


hdulist=fits.open(infilename , ignore_missing_end=True)
hdulist.info()
h=hdulist[1]
print 'Reading in data...'
data=h.data

nentries = len(data['z'])

print "Read in %d entries." % (nentries)
'''
################################################################################
# BOSS TARGET CUT
################################################################################
print 'Selecting BOSS Targets'
index_bt1 = data['BOSS_TARGET1']>0
npass = len(index_bt1[index_bt1==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

################################################################################
#CMASS CUT
################################################################################
print 'Selecting CMASS targets'
bt=data['BOSS_TARGET1']
print bt[bt>0]
place=1
##################### CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!
digit1=2**place; digit2=2**(place+1)
cmassbin=(bt%digit2)/digit1
index_cmass = cmassbin==1
print index_cmass
npass = len(index_cmass[index_cmass==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

################################################################################
print 'Selecting ifib2<21.5 targets'
#bt=data['BOSS_TARGET1']
place=8
##################### CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!
digit1=2**place; digit2=2**(place+1)
ifib2bin=(bt%digit2)/digit1
index_ifib2 = ifib2bin==0
print index_ifib2
npass = len(index_ifib2[index_ifib2==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 
'''
################################################################################
#GALAXY CUT
################################################################################
#print 'Selecting Galaxies'
#index_gal=data['CLASS']=='GALAXY'
#npass = len(index_gal[index_gal==True])
#print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
#print 

#galnew=data['CLASS'][gal]
#notgal=len(data)-len(galnew)
#print " %d targets are not known to be galaxies " % notgal

################################################################################
# PRIMARY OBSERVATION CUT                
################################################################################
#print 'Primary Observation Cut'
#index_po=data['SPECPRIMARY']==1
#npass = len(index_po[index_po==True])
#print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
#print 

#ponew=data['SPECPRIMARY'][po]
#notpo=len(data)-len(ponew)
#print " %d targets were not primarily observed " % notpo

################################################################################
# WARNING VALUES CUT
################################################################################
#print 'Cutting out Warning Values (selecting ones with value==0)'
#index_warning=data['ZWARNING_NOQSO']==0
#npass = len(index_warning[index_warning==True])
#print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
#print 

#wvnew=data['ZWARNING_NOQSO'][wv]
#notwv=len(data)-len(wvnew)
#print " %d targets contained redshift errors" % notwv

################################################################################
# Z CUT #
################################################################################
zlo = 0.43
zhi = 0.70

print 'Selecting values with redshift z %f-%f' % (zlo,zhi)
index_zcut=data['Z']>zlo
index_zcut *= data['Z']<zhi
npass = len(index_zcut[index_zcut==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

################################################################################
# Z CUT #
################################################################################
ralo = 90.
rahi = 270.

print 'Selecting values with right ascension %f-%f' % (ralo,rahi)
index_racut=data['RA']>ralo
index_racut *= data['RA']<rahi
npass = len(index_racut[index_racut==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

################################################################################
print "Calculating the total cuts!"
################################################################################
index_total = index_zcut*\
              index_racut
              #index_cmass
              #index_po
              #index_warning
              #index_gal

npass = len(index_total[index_total==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

totcmass=data[index_total]

points = 200000
aa=np.arange(0,len(totcmass))
np.random.shuffle(aa)
finaldat=totcmass[aa[0:points]]


ra = finaldat['RA']
dec = finaldat['DEC']
z = finaldat['Z']
sys = finaldat['weight_systot']
cp = finaldat['weight_cp']
rf = finaldat['weight_noz']
weight = sys*(rf + cp - 1)
print len(weight)
print "Writing out data...\n"
mockdata=np.column_stack((ra,dec,z,weight))
print len(mockdata)
outname = "%s.dat" % (tag)
np.savetxt(outname,mockdata)
