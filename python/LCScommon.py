#!/usr/bin/env python
from pylab import *
import os
from scipy.stats.stats import spearmanr
from scipy.stats import ks_2samp


pscale24=2.45#arcsec per pixel

pscalesdss=1.#arcsec per pixel
sdsspixelscale=0.396127#conversion for isophotal radii from pixels to arcseconds
mipspixelscale=pscale24
mipsconv_MJysr_to_uJy=141.086

mingalaxysize=2.*mipspixelscale


clusternames=['MKW11', 'MKW8', 'AWM4', 'A2063', 'A2052', 'NGC6107', 'Coma', 'A1367', 'Hercules']

clusterRA={'MKW11':202.3800, 'MKW8':220.1796, 'AWM4':241.2375, 'A2063':230.7578, 'A2052':229.1896, 'NGC6107':244.333750, 'Coma':194.9531, 'A1367':176.1231, 'Hercules':241.3125}
clusterDec={'MKW11':11.78861, 'MKW8':3.4530, 'AWM4':23.9206, 'A2063':8.6394, 'A2052':7.0003, 'NGC6107':34.901389, 'Coma':27.9807, 'A1367':19.8391, 'Hercules':17.7485}

clustervel={'MKW11':6854., 'MKW8':8100., 'AWM4':9526., 'A2063':10481., 'A2052':10647., 'NGC6107':9197., 'Coma':6900., 'A1367':8400., 'Hercules':11100.}

clustersigma={'MKW11':361, 'MKW8':325., 'AWM4':500., 'A2063':660., 'A2052':562., 'NGC6107':500., 'Coma':1000., 'A1367':745., 'Hercules':689.}

clusterf80MJysr={'MKW11':4., 'MKW8':3.75, 'AWM4':3.5, 'A2063':4., 'A2052':4., 'NGC6107':3.25, 'Coma':2.25, 'A1367':3.5, 'Hercules':3.25}

clusterz={'MKW11':.022849,'MKW8':.027,'AWM4':.031755,'A2063':.034937,'A2052':.035491,'NGC6107':.030658,'Coma':.023,'A1367':.028,'Hercules':.037}

john_prefix={'MKW11':'mkw11','MKW8':'mkw8','AWM4':'awm4','A2063':'abell2063','A2052':'abell2052','NGC6107':'ngc6107','Coma':'coma','A1367':'abell1367','Hercules':'hercules'}

xraycontourlevels={'MKW11':[.85,1.69,2.54],'MKW8':[.49,.99,1.48,1.98],'AWM4':[.8,1.6,2.4],'NGC6107':[1.43,2.85,4.27],'A2052':[.9,1.8,2.7,3.6],'A2063':[.9,1.8,2.7,3.6],'Hercules':[.9,1.92,2.9,3.8],'A1367':[.6,1.17,1.76,2.35],'Coma':[.88,1.76,2.63,3.51]}#used contour option in ds9 to derive these

  

spiral_nozoo={'MKW11':[70685, 143485, 143530],'MKW8':[18127],'AWM4':[68283,  68288,  68338,  68341, 166624],'NGC6107':[43707,  43712,  43787,  43857,  69538],'A2052':[79646,145994,166042],'A2063':[166124],'Hercules':[99840, 146607],'A1367':[140124,140145, 140176, 140194],'Coma':[104125, 104232,142527,142561,142563,142568,142572, 142627,142642,142653,142656,142663,142666,142668,142669,142676,142682,142687,142689,142693,142695,142706,142710,142713,142715,142716,142723,142727,142737,142740,142745,142750,142755,142758,142765,142767,142769,142774,142779,142781,142793,142795,142801,142804,142806,142808,142809,142810,142815,142819,142825,142837,142847,142855,142873,142914,162740]}#used contour option in ds9 to derive these

elliptical_nozoo={'MKW11':[],'MKW8':[],'AWM4':[],'NGC6107':[146832, 146876, 146878,146880, 166860],'A2052':[79600,  79610, 79705,79710,146012, 146037, 146041],'A2063':[72627,72751,146137],'Hercules':[146638,146664],'A1367':[113076,113458,  140164],'Coma':[103978,104022,104061,104115,142531,142552,142584,142585,142604,142605,142609,142611,142614,142615,142616,142622,142623,142628,142636,142637,142638,142647,142648,142649,142651,142658,142660,142661,142675,142677,142678,142681,142684,142690,142699,142705,142717,142721,142725,142729,142741,142743,142761,142787,142803,142813,142832,142852,142866,162659]}#used contour option in ds9 to derive these

irreg_nozoo={'MKW11':[],'MKW8':[],'AWM4':[],'NGC6107':[],'A2052':[79550, 79593, 79665],'A2063':[],'Hercules':[146635,146673, 146680, 166679],'A1367':[140170, 140177, 140183, 140184, 140186, 160496],'Coma':[142559,142560,142578,142590,142593,142613,142620,142631,142645,142652,142667,142673,142679,142697,142718,142733,142753,142762,142771,142786,142821,142823,142826,142831,142834,142849,162689]}#used contour option in ds9 to derive these

# galaxies to cut from sample, at least initially
# galaxies that have contamination be nearby neighbor.  see notes.
visual_cut={'MKW11':[70639,70694,143485,171004],'MKW8':[18111, 18171],'AWM4':[82134, 82188, 82209, 146626, 166655, 166699],'NGC6107':[43782, 43814, 69617, 69618],'A2052':[ 79388, 166086],'A2063':[72631, 72710, 72745, 72782, 146106, 146107, 146124, 146128, 146130,146135],'Hercules':[99056, 99644, 99822, 99859, 99872, 146607, 146659, 166638],'A1367':[113058, 113404,  140197],'Coma':[103612, 103628, 103648, 103784, 103831, 103833, 103844, 103924,  103933, 104001, 104004, 104035, 104126, 142655, 142840, 162793, 162831]}


#Group names
groupnames=['NRGb041','NRGb151','NRGb157','NRGb168','NRGb206','NRGb247','NRGb282','NRGb301','MKW8','NCG5846','NRGs076','NRGs272','NRGs385']
altgroupnames=['WBL226','MKW10','HCG59','WBL368','WBL404','MKW11test','Zw1400','WBL509','MKW8','NGC5846','WBL251','WBL477','NGC6107']
#location of Final images


# central biweight location as calculated from findbestbiweight code
clusterbiweightcenter={'MKW11':6906,'MKW8':8098,'AWM4':9650,'A2063':10422,'A2052':10354.5,'NGC6107':9429,'Coma':6999,'A1367':6481,'Hercules':10957.5}

 #sbi values output from +/- 4000km/s and 1 degree velocity cut from findbestbiweight code
clusterbiweightscale={'MKW11':392.37,'MKW8':491.32,'AWM4':476.67,'A2063':727.06,'A2052':626.32,'NGC6107':616.86,'Coma':937.03,'A1367':794.61,'Hercules':772.74}

# X-ray luminosity in 10^43 ergs/s
# from Bohringer et al 2000, and Mahdavi et Geller
clusterLx={'MKW11':0.033,'MKW8':0.096,'AWM4':0.550,'A2063':1.940,'A2052':2.580,'NGC6107':0.33,'Coma':7.010,'A1367':1.510,'Hercules':0.980}

# Tx, errdown, errup
clusterTx={'MKW11':[0],'MKW8':[3.,.12,.12],'AWM4':[0],'A2063':[3.77,.06,.06],'A2052':[3.35,.02,.02],'NGC6107':[0],'Coma':[9.15,.17,.17],'A1367':[3.58,.06,.06],'Hercules':[0]}
# X-ray temp in keV; from Mittal et al 2011

clusterLx2={'MKW11':0.033,'MKW8':(.692,.058),'AWM4':0.550,'A2063':(2.06,.027),'A2052':(2.18,.022),'NGC6107':0.083,'Coma':(11.1,.156),'A1367':(1.13,.009),'Hercules':0.980}
# list of L500 (1.e37 W), M500(1.e14 Msun) and R500 (Mpc) from Piffaretti+ 2011
clusterXray={'MKW11':[0.065077,	0.3805,	0.5078],'MKW8':[0.192567,0.7352,0.6316],'AWM4':[0.284521,0.9289,0.6815],'A2063':[1.138819,2.1598,0.9020],'A2052':[1.442058,2.4945,0.9465],'NGC6107':[0.168099,0.6744,0.6127],'Coma':[3.455556,4.2846,1.1378],'A1367':[1.104603,2.1398,0.9032],'Hercules':[0.508824,1.3202,0.7652]}


# these correpond to area w/uniform 24um coverage
# center x,y,dx,dy,rotation E of N, all in degrees
cluster24Box={'MKW11':array([202.36239,11.752736,1.3138054,3.046197,27.0001],'f'), 'MKW8':array([220.18764,3.4955922,1.3188409,3.040413,13.5],'f'), 'AWM4':array([241.21434,23.872723,1.3441978,3.0241238,10],'f'), 'A2063':array([230.77172,8.6817732,1.3126447,3.0415136,13.5001],'f'), 'A2052':array([229.19761,7.0403283,1.3194664,3.0412907,13.25],'f'), 'NGC6107':array([244.30039,34.934184,1.3199655,3.0435265,322],'f'), 'Coma':array([194.86318,27.865896,1.5391027,1.976467,29.5002],'f'), 'A1367':array([176.1019,19.799614,.51080152,.90025557,31.5],'f'), 'Hercules':array([241.3065,17.771646,.51029561,.93431905,19.5001],'f')}


#solar magnitude in SDSS filters
SolarMag={'u':6.39,'g':5.07,'r':4.62,'i':4.52,'z':4.48}

#cosmology
H0=70
OmegaL=0.7
OmegaM=0.3
h=H0/100.

#bell stellar mass coefficients for sdss filters
bellug={'g':[-.221,0.485],'r':[-.099,0.345],'i':[-.053,0.268],'z':[-.105,0.226]}
bellur={'g':[-.390,0.417],'r':[-.223,0.229],'i':[-.151,0.233],'z':[-.178,0.192]}
bellui={'g':[-.375,0.359],'r':[-.212,0.257],'i':[-.144,0.201],'z':[-.171,0.165]}
belluz={'g':[-.400,0.332],'r':[-.232,0.239],'i':[-.161,0.187],'z':[-.179,0.151]}
bellgr={'g':[-.499,1.519],'r':[-.306,1.097],'i':[-.222,0.864],'z':[-.223,0.689]}
bellgi={'g':[-.379,0.914],'r':[-.220,0.661],'i':[-.152,0.518],'z':[-.175,0.421]}
bellgz={'g':[-.367,0.698],'r':[-.215,0.508],'i':[-.153,0.402],'z':[-.171,0.322]}
bellri={'g':[-.106,1.982],'r':[-.022,1.431],'i':[0.006,1.114],'z':[-.952,0.923]}
bellrz={'g':[-.124,1.067],'r':[-.041,0.780],'i':[-.018,0.623],'z':[-.041,0.463]}

snr24cut=2.
deltaCutout=100.#width of cutouts in arcsec
ramin=170.#cuts for culling the ac
ramax=250.#cuts for culling the ac
decmin=0.
decmax=38.#cuts for culling the ac
zmin=0.01366#min z cut, z(coma)-3 sigma
zmax=0.04333#max z cut, z(A2052, which is  10900 km/s)+ 4*sigma
vmin=zmin*3.e5
vmax=zmax*3.e5
#cutoutpath='/home/rfinn/research/LocalClusters/cutouts/'
cutoutpath='/home/rfinn/research/LocalClusters/cutouts/'

Lsol=3.826e33#normalize by solar luminosity
bellconv=9.8e-11#converts Lir (in L_sun) to SFR/yr
bellconv=4.5e-44#Kenn 98 conversion fro erg/s to SFR/yr
catalog_radial_cut = 3. # mastertable radial cut in degrees
mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/rfinn/'

mipsflux2umJyconv=141.086


def multiplotaxes(i):
    ax=gca()
    noylabel=[2,3,5,6,8,9]
    if i < 7:
        ax.set_xticklabels(([]))
    if i in noylabel:
        ax.set_yticklabels(([]))
def multiplotlabels(xl,yl):
    ax=gca()
    text(-.5,-.25,xl,fontsize=22,horizontalalignment='center',transform=ax.transAxes)
    text(-2.4,1.5,yl,fontsize=22,verticalalignment='center',rotation=90,transform=ax.transAxes,family='serif')


def spearman(x,y):
    rho,pvalue=spearmanr(x,y)
    print 'Spearman Rank Test:'
    print 'rho = %6.2f'%(rho)
    print 'p-vale = %6.5f (prob that samples are uncorrelated)'%(pvalue) 
    return rho,pvalue

def ks(x,y):
    D,pvalue=ks_2samp(x,y)
    print 'KS Test:'
    print 'D = %6.2f'%(D)
    print 'p-vale = %6.5f (prob that samples are from same distribution)'%(pvalue) 
    return D,pvalue
    
def findnearest(x1,y1,x2,y2,delta):#use where command
	matchflag=1
	nmatch=0
	d=sqrt((x1-x2)**2 + (y1-y2)**2)#x2 and y2 are arrays
	index=arange(len(d))
	t=index[d<delta]
	matches=t
	if len(matches) > 0:
		nmatch=len(matches)
		if nmatch > 1:
			imatch=index[(d == min(d[t]))]
		else:
			imatch=matches[0]			
	else:
		imatch = 0
		matchflag = 0

	return imatch, matchflag,nmatch

def drawbox(data,style):#feed in center x,y,dx,dy,rotation E of N
    #xcoords of unrotated box, going around CCW
    xl=array([data[0]-0.5*data[2],data[0]+0.5*data[2],data[0]+0.5*data[2],data[0]-0.5*data[2],data[0]-0.5*data[2]],'d')
    yl=array([data[1]-0.5*data[3],data[1]-0.5*data[3],data[1]+0.5*data[3],data[1]+0.5*data[3],data[1]-0.5*data[3] ],'d')

    xl=array([-0.5*data[2],+0.5*data[2],+0.5*data[2],-0.5*data[2],-0.5*data[2]],'d')
    yl=array([-0.5*data[3],-0.5*data[3],+0.5*data[3],+0.5*data[3],-0.5*data[3] ],'d')

    ang=data[4]*pi/180.*-1.#convert rotation to radians
    #rotate coordinates
    xp=cos(ang)*xl-sin(ang)*yl
    yp=sin(ang)*xl+cos(ang)*yl

    #put back on absolute scale
    xp=data[0]+xp
    yp=data[1]+yp
    #draw rotated box
    plot(xp,yp,style)

