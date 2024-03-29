import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import matplotlib.pylab as plt
#hdulist=fits.open('specObj-dr12.fits' , ignore_missing_end=True)

infilename = sys.argv[1]
hdulist=fits.open(infilename , ignore_missing_end=True)
hdulist.info()
h=hdulist[1]
print 'Reading in Data'
data=h.data

#Test#
#data=data[0:10000]
'''
print 'Data is Read'
##### CMASS Galaxies ######
'''
# BOSS TARGET CUT
print 'Cutting BOSS Targets'
indbt1=data['BOSS_TARGET1']>0
#bt1new=data['BOSS_TARGET1'][indbt1]
notboss=len(data)-len(indbt1[indbt1==True])
print " %d targets are not BOSS Targets " % notboss

#CMASS CUT
print 'Cutting CMASS Targets'
bt=data['BOSS_TARGET1']
place=1
##################### CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!
digit1=2**place; digit2=2**(place+1)
cmassbin=(bt%digit2)/digit1
cmass=cmassbin==1
cmassnew=cmassbin[cmass]
notcmass=len(data)-len(cmassnew)
print " %d targets are not CMASS Targets " % notcmass

#GALAXY CUT
print 'Cutting Galaxies'
gal=data['CLASS']=='GALAXY'
galnew=data['CLASS'][gal]
notgal=len(data)-len(galnew)
print " %d targets are not known to be galaxies " % notgal

# PRIMARY OBSERVATION CUT                
print 'Primary Observation Cut'
po=data['SPECPRIMARY']==1
ponew=data['SPECPRIMARY'][po]
notpo=len(data)-len(ponew)
print " %d targets were not primarily observed " % notpo

# WARNING VALUES CUT
print 'Cutting Warning Values'
wv=data['ZWARNING_NOQSO']==0
wvnew=data['ZWARNING_NOQSO'][wv]
notwv=len(data)-len(wvnew)
print " %d targets contained redshift errors" % notwv
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
# Z CUT #
print 'Z Cut'
zcut=data['Z']>0.43
zcutnew=data['Z'][zcut]

zcut1=data['Z']<0.70
zcutnew1=data['Z'][zcut1]
notzcut=(len(data)-len(zcutnew))+(len(data)-len(zcutnew1))
print " %d targets failed the Z cut" % notzcut
'''
'''
tot=indbt1*cmass*po*wv*gal*zcut1*zcut
totcmass=data[tot]

'''
# TOTAL CUT
print 'Finalizing CMASS Cut'
tot1=indbt1
print " %d BOSS targets " % len(data[tot1])
tot2=cmass*tot1
print " %d of those are CMASS targets" % len(data[tot2])
tot3=tot2*gal
print " %d of those are also galaxies" % len(data[tot3])
tot4=po*tot3
print " %d of those were primarily observed" % len(data[tot4])
tot5=wv*tot4
print " %d of those did not contain redshift errors " % len(data[tot5])
tot6=c1*c2*tot5
print " %d of those had standard flag meanings" % len(data[tot6])
tot7=ifib*tot6
print " %d of those had a fiber magnitude under 21.5" % len(data[tot7])
tot8=tot7*zcut1*zcut
print " %d of those passed the Z cut" % len(data[tot8])

totcmass=data[tot8]
'''
#print "Remaining CMASS Galaxies: %d" % len(totcmass)
plt.plot(totcmass['PLUG_RA'],totcmass['PLUG_DEC'],'o',markersize=0.2)
plt.xlim(90,280)
plt.ylim(-10,80)
plt.title('BOSS,CMASS,Primary Obs,ZWarning,Galaxy, and Z cut')
plt.show()
##### LOWZ GALAXIES #####

'''
# BOSS TARGET CUT
print 'Cutting BOSS Targets'
indbt1=data['BOSS_TARGET1']>0
bt1new=data['BOSS_TARGET1'][indbt1]
notboss=len(data)-len(bt1new)
print " %d targets are not BOSS Targets " % notboss

#LOWZ CUT
print 'Cutting LOWZ Targets'
bt=data['BOSS_TARGET1']
place=0
digit1=2**place; digit2=2**(place+1)
lowzbin=(bt%digit2)/digit1
lowz=lowzbin==1
lowznew=lowzbin[lowz]
notlowz=len(data)-len(lowznew)
print " %d targets are not LOWZ Targets " % notlowz

#GALAXY CUT
print 'Cutting Galaxies'
gal=data['CLASS']=='GALAXY'
galnew=data['CLASS'][gal]
notgal=len(data)-len(galnew)
print " %d targets are not known to be galaxies " % notgal

# PRIMARY OBSERVATION CUT
print 'Primary Observation Cut'
po=data['SPECPRIMARY']==1
ponew=data['SPECPRIMARY'][po]
notpo=len(data)-len(ponew)
print " %d targets were not primarily observed " % notpo

# WARNING VALUES CUT
print 'Cutting Warning Values'
wv=data['ZWARNING_NOQSO']==0
wvnew=data['ZWARNING_NOQSO'][wv]
notwv=len(data)-len(wvnew)
print " %d targets contained redshift errors" % notwv

#IFIBER CUT
print 'IFiber Cut'
place=8
digit1=2**place; digit2=2**(place+1)
ifibbin=(bt%digit2)/digit1
ifib=ifibbin==0
ifibnew=ifibbin[ifib]
notifib=len(data)-len(ifibnew)
print " %d targets had a fiber magnitude over 21.5 " % notifib

#TILE ID CUT
print 'Tile ID Cut'
tile=data['TILE']>=10324
tilenew=data['TILE'][tile]
nottile=len(data)-len(tilenew)
print " %d targets were incorrectly targeted" % nottile

# TOTAL CUT
print 'Finalizing CMASS Cut'
tot=indbt1
print " %d BOSS targets " % len(data[tot])
tot*=lowz
print " %d of those are LOWZ targets" % len(data[tot])
tot*=gal
print " %d of those are also galaxies" % len(data[tot])
tot*=po
print " %d of those were primarily observed" % len(data[tot])
tot*=wv
print " %d of those did not contain redshift errors " % len(data[tot])
tot*=tile
print " %d of those were correctly targeted" % len(data[tot])
tot*=ifib
print " %d of those had a fiber magnitude under 21.5" % len(data[tot])
totlowz=data[tot]
print "Remaining LOWZ Galaxies: %d" % len(totlowz)



'''
'''
#####Histogram######
import matplotlib.pylab as plt
print 'Creating CMASS Histogram'

plt.figure()
plt.hist(data['Z'],bins=100,range=(0.0,1.0),label='Total CMASS Data',facecolor='orange') 
print 'a'
plt.hist(data[tot1]['Z'],bins=100,range=(0.0,1.0),label='BOSS Targets',facecolor='yellow')
print 'b'
plt.hist(data[tot2]['Z'],bins=100,range=(0.0,1.0),label='CMASS Targets',facecolor='magenta')
print 'c'
plt.hist(data[tot3]['Z'],bins=100,range=(0.0,1.0),label='Galaxies',facecolor='cyan')
print 'd'
plt.hist(data[tot4]['Z'],bins=100,range=(0.0,1.0),label='Primarily Observed', facecolor='red')
print 'e'
plt.hist(data[tot5]['Z'],bins=100,range=(0.0,1.0),label='No Redshit Error',facecolor='green')
print 'f'
plt.hist(data[tot6]['Z'],bins=100,range=(0.0,1.0),label='Standard Flag Meanings',facecolor='blue')
print 'g'
plt.hist(data[tot7]['Z'],bins=100,range=(0.0,1.0),label='Final CMASS Targets',facecolor='black')
print 'h'
plt.ylim(0,75000)
plt.xlabel('RedShift Value')
plt.ylabel('Frequency')
plt.title('Z Values of SDSS DR12 Data and Corresponding CMASS Target Cuts')
plt.legend(prop={'size':6.5})
plt.show()


print 'Making Plot'
##### Scatter Plot ######
import matplotlib.pylab as plt
plt.plot(totcmass['PLUG_RA'],totcmass['PLUG_DEC'],'o',markersize=0.2)
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Declination (Degrees)')
plt.title('CMASS DR10 Data')
plt.xlim(90,280)
plt.ylim(-10,80)
plt.show()
'''




