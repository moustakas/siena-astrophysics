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

nentries = len(data['BOSS_TARGET1'])

print "Read in %d entries." % (nentries)

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
#bt=data['BOSS_TARGET1']
#place=1
##################### CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!
#digit1=2**place; digit2=2**(place+1)
#cmassbin=(bt%digit2)/digit1
#cmass=cmassbin==1
#cmassnew=cmassbin[cmass]
#notcmass=len(data)-len(cmassnew)

bt=data['BOSS_TARGET1']
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
#GALAXY CUT
################################################################################
print 'Selecting Galaxies'
index_gal=data['CLASS']=='GALAXY'
npass = len(index_gal[index_gal==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

#galnew=data['CLASS'][gal]
#notgal=len(data)-len(galnew)
#print " %d targets are not known to be galaxies " % notgal

################################################################################
# PRIMARY OBSERVATION CUT                
################################################################################
print 'Primary Observation Cut'
index_po=data['SPECPRIMARY']==1
npass = len(index_po[index_po==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

#ponew=data['SPECPRIMARY'][po]
#notpo=len(data)-len(ponew)
#print " %d targets were not primarily observed " % notpo

################################################################################
# WARNING VALUES CUT
################################################################################
print 'Cutting out Warning Values (selecting ones with value==0)'
index_warning=data['ZWARNING_NOQSO']==0
npass = len(index_warning[index_warning==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

#wvnew=data['ZWARNING_NOQSO'][wv]
#notwv=len(data)-len(wvnew)
#print " %d targets contained redshift errors" % notwv
'''
# CHUNK CUT
print 'Chunk Cut'
c1=data['CHUNK']!="boss1"
c1new=data['CHUNK'][c1]
c2=data['CHUNK']!="boss2"
c2new=data['CHUNK'][c2]
notc=(len(data)-len(c1new))+(len(data)-len(c2new))
print " %d targets were removed due to differing flag meanings " % notc

#IFIBER CUT
print 'IFiber Cut'
place=8
digit1=2**place; digit2=2**(place+1)
ifibbin=(bt%digit2)/digit1
ifib=ifibbin==0
ifibnew=ifibbin[ifib]
notifib=len(data)-len(ifibnew)
print " %d targets had a fiber magnitude over 21.5 " % notifib
'''
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
index_racut=data['PLUG_RA']>ralo
index_racut *= data['PLUG_RA']<rahi
npass = len(index_racut[index_racut==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

################################################################################
print "Calculating the total cuts!"
################################################################################
index_total = index_bt1*\
              index_cmass*\
              index_po*\
              index_warning*\
              index_gal*\
              index_zcut*\
              index_racut

npass = len(index_total[index_total==True])
print "%d (%f) targets pass cuts" % (npass,npass/float(nentries))
print 

totcmass=data[index_total]

ra = totcmass['PLUG_RA']
dec = totcmass['PLUG_DEC']
z = totcmass['Z']

print "Writing out data...\n"
mockdata=np.column_stack((ra,dec,z))
print len(mockdata)
outname = "%s.dat" % (tag)
np.savetxt(outname,mockdata)


print 'Making Plot'
##### Scatter Plot ######
import matplotlib.pylab as plt
plt.plot(ra,dec,'o',markersize=0.2)
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Declination (Degrees)')
plt.title('CMASS DR10 Data')
#plt.xlim(90,280)
#plt.ylim(-10,80)
plt.show()



