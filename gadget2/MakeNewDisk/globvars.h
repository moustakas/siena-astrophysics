
#define RMASSBINS   2048
#define ZMASSBINS   32
#define PHIMASSBINS 64


#define  GAMMA (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define  PI  3.1415926

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10


#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  HYDROGEN_MASSFRAC  0.76

#define  HUBBLE      1.0     /* Hubble constant in 100km/sec/Mpc */  



extern int NumPart;  /* this is the particle number used for the tree */
extern struct part_data
{
  float Pos[3];
  float Mass;
}
*P;
extern double NumGasScaleLengths;

extern double U4;



extern double ErrTolTheta;

extern double MinGasTemp;

extern double MaxGasDiskHeight;
 
/*** integration parameters ***/


#define  RSIZE 512
#define  ZSIZE 512


#define NSHEETS 100  /* 100 */


#define TABSIZE 512

extern double  rhocentral[TABSIZE];
extern double  sigmatab[TABSIZE];

extern double FactorEVP;
extern double PhysDensThresh;
extern double EgySpecSN;
extern double MaxSfrTimescale;
extern double FactorSN;
extern double TempSupernova;
extern double TempClouds;
extern double EgySpecCold;
extern double EgySpecSN;
extern double FeedbackEnergy;
extern double FactorForSofterEQS;
extern double RadialDispersionFactor;



/***********  INPUT PARAMETERS *********/

extern double  CC;      /* halo concentration */
extern double  V200;    /* circular velocity v_200 */
extern double  LAMBDA;    /* spin parameter  */
extern double  MD;        /* disk mass fraction */
extern double  JD;        /* disk spin fraction */
extern double  MB;        /* bulge mass fraction */
extern double  MBH;       /* central black hole mass fraction */
extern double  GasFraction;  
extern double  DiskHeight; 
extern double  BulgeSize;
extern double  HI_GasMassFraction;    /* in terms of the total gas mass */
extern double  HI_GasDiskScaleLength;  /* in terms of scale length of the disk */ 

extern int     N_HALO;    /* desired number of particles in halo */
extern int     N_DISK;    /* desired number of collsionless particles in disk */
extern int     N_GAS;     /* number of gas particles in stellar disk */ 
extern int     N_BULGE;   /* number of gas particles in stellar disk */ 
extern int     N_BLACKHOLE;

extern double  REDSHIFT;  /*redshift to scale galaxy to, added by Brant 10/1/04 */
extern double  Omega_m0;  /*Omega_m(z=0), added by Brant 10/1/04 */
extern double  Omega_L0;  /*Omega_L(z=0), added by Brant 10/1/04 */
extern double  Hz;        /*Hubble param at z=REDSHIFT, added by Brant 10/1/04 */
/*********************************************/







extern double  M200;     /* virial mass */

extern double  RH;       /* scale radius in Hernquist profile */

extern double  RS;       /* scale radius for halo */
extern double  R200;     /* virial radius */
extern double  H;        /* disk scale length */
extern double  Z0;       /* disk thickness */
extern double  A;        /* bulge scale radius */


extern double  M_HALO;      /* total dark mass */
extern double  M_DISK;      /* mass of stellar disk (collisionless part) */
extern double  M_GAS;       /* gas mass in disk */
extern double  M_BULGE;     /* mass of bulge */
extern double  M_BLACKHOLE; /* mass of black-hole */

extern double  halo_spinfactor;  /* computed streamin of dark matter */




extern double HubbleParam;

extern double G;            /* gravitational constant */
extern double H0;           /* Hubble constant */
extern double UnitTime_in_s;
extern double UnitMass_in_g;
extern double UnitLength_in_cm;
extern double UnitVelocity_in_cm_per_s;
extern double UnitTime_in_Megayears;


extern double UnitPressure_in_cgs,
    UnitDensity_in_cgs,
    UnitCoolingRate_in_cgs,
    UnitEnergy_in_cgs;


extern char OutputDir[500], OutputFile[500];






/* particle data */

extern double    *vmax2_halo,*vmax2_disk,*vmax2_bulge,*vmax2_gas;

extern double    *xp_halo,*yp_halo,*zp_halo,*mp_halo;
extern double    *xp_disk,*yp_disk,*zp_disk,*mp_disk;
extern double    *xp_bulge,*yp_bulge,*zp_bulge,*mp_bulge;
extern double    *xp_gas,*yp_gas,*zp_gas,*mp_gas,*u_gas;

extern double    *vxp_halo,*vyp_halo,*vzp_halo;
extern double    *vxp_disk,*vyp_disk,*vzp_disk;
extern double    *vxp_bulge,*vyp_bulge,*vzp_bulge;
extern double    *vxp_gas,*vyp_gas,*vzp_gas;



extern double  LL;       /* LL = extension of fields in R and z. */
extern double  Baselen;



extern double **RhoGas, **CumulMassGas;
extern double *Zrho, *Zcumul;


extern double **Dphi_z,**Dphi_R,**Dphi_z_dR;  /* derivatives of total potential */

extern double *epi_gamma2,*epi_kappa2;  /* epicycle gamma^2  */ 



/* halo velocity fields */

extern double **VelDispRz_halo;
extern double **VelDispPhi_halo;
extern double **VelVc2_halo;
extern double **VelStreamPhi_halo;
extern double **VelDispRz_dR_halo;


/* bulge velocity fields */

extern double **VelDispRz_bulge;
extern double **VelDispPhi_bulge;
extern double **VelVc2_bulge;
extern double **VelStreamPhi_bulge;
extern double **VelDispRz_dR_bulge;


/* disk velocity fields */

extern double **VelDispRz_disk;
extern double **VelDispPhi_disk;
extern double **VelVc2_disk;
extern double **VelStreamPhi_disk;
extern double **VelDispRz_dR_disk;


/* gas velocity fields */

extern double **VelStreamGas;


/* auxiliary field */

extern double *xl,*yl,*D2yl;
extern double *list_z,*list_R,*list_RplusdR;






