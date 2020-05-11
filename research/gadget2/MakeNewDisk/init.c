
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"





void init_units(void)
{

  UnitMass_in_g = 1e10 * SOLAR_MASS;	/* 10^10 solar masses */
  UnitLength_in_cm = CM_PER_MPC / 1000;	/* 1 kpc */
  UnitVelocity_in_cm_per_s = 1e5;	/* 1 km/sec */


  HubbleParam = 0.7;


  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  H0 = HUBBLE * 100 * 1e5 / CM_PER_MPC / UnitVelocity_in_cm_per_s * UnitLength_in_cm;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);

  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;

  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);


  printf("%g %g %g \n", G, H0, UnitTime_in_Megayears);
}







void init(void)
{
  int i;


  if(N_HALO > 0)
    {
      xp_halo = vector(1, N_HALO);
      yp_halo = vector(1, N_HALO);
      zp_halo = vector(1, N_HALO);
      mp_halo = vector(1, N_HALO);
      vxp_halo = vector(1, N_HALO);
      vyp_halo = vector(1, N_HALO);
      vzp_halo = vector(1, N_HALO);
      vmax2_halo = vector(1, N_HALO);
    }

  if(N_DISK > 0)
    {
      xp_disk = vector(1, N_DISK);
      yp_disk = vector(1, N_DISK);
      zp_disk = vector(1, N_DISK);
      mp_disk = vector(1, N_DISK);
      vxp_disk = vector(1, N_DISK);
      vyp_disk = vector(1, N_DISK);
      vzp_disk = vector(1, N_DISK);
      vmax2_disk = vector(1, N_DISK);
    }

  if(N_BULGE > 0)
    {
      xp_bulge = vector(1, N_BULGE);
      yp_bulge = vector(1, N_BULGE);
      zp_bulge = vector(1, N_BULGE);
      mp_bulge = vector(1, N_BULGE);
      vxp_bulge = vector(1, N_BULGE);
      vyp_bulge = vector(1, N_BULGE);
      vzp_bulge = vector(1, N_BULGE);
      vmax2_bulge = vector(1, N_BULGE);
    }

  if(N_GAS > 0)
    {
      xp_gas = vector(1, N_GAS);
      yp_gas = vector(1, N_GAS);
      zp_gas = vector(1, N_GAS);
      mp_gas = vector(1, N_GAS);
      vxp_gas = vector(1, N_GAS);
      vyp_gas = vector(1, N_GAS);
      vzp_gas = vector(1, N_GAS);
      u_gas = vector(1, N_GAS);
      vmax2_gas = vector(1, N_GAS);
    }



  Dphi_z = matrix(0, RSIZE, 0, ZSIZE);
  Dphi_R = matrix(0, RSIZE, 0, ZSIZE);
  Dphi_z_dR = matrix(0, RSIZE, 0, ZSIZE);

  VelDispRz_halo = matrix(0, RSIZE, 0, ZSIZE);
  VelDispRz_dR_halo = matrix(0, RSIZE, 0, ZSIZE);
  VelDispPhi_halo = matrix(0, RSIZE, 0, ZSIZE);
  VelVc2_halo = matrix(0, RSIZE, 0, ZSIZE);
  VelStreamPhi_halo = matrix(0, RSIZE, 0, ZSIZE);

  VelDispRz_disk = matrix(0, RSIZE, 0, ZSIZE);
  VelDispRz_dR_disk = matrix(0, RSIZE, 0, ZSIZE);
  VelDispPhi_disk = matrix(0, RSIZE, 0, ZSIZE);
  VelVc2_disk = matrix(0, RSIZE, 0, ZSIZE);
  VelStreamPhi_disk = matrix(0, RSIZE, 0, ZSIZE);


  VelDispRz_dR_bulge = matrix(0, RSIZE, 0, ZSIZE);
  VelDispRz_bulge = matrix(0, RSIZE, 0, ZSIZE);
  VelDispPhi_bulge = matrix(0, RSIZE, 0, ZSIZE);
  VelVc2_bulge = matrix(0, RSIZE, 0, ZSIZE);
  VelStreamPhi_bulge = matrix(0, RSIZE, 0, ZSIZE);

  VelStreamGas = matrix(0, RSIZE, 0, ZSIZE);

  RhoGas = matrix(0, RSIZE, 0, ZSIZE);
  CumulMassGas = matrix(0, RSIZE, 0, ZSIZE);

  Zrho = vector(0, ZSIZE);
  Zcumul = vector(0, ZSIZE);


  if(N_GAS > 0 && N_DISK > 0)
    NumPart = 2 * RMASSBINS * ZMASSBINS * PHIMASSBINS;
  else
    NumPart = RMASSBINS * ZMASSBINS * PHIMASSBINS;

  P = malloc(sizeof(struct part_data) * NumPart);

  force_treeallocate(1.5 * NumPart, NumPart);



  list_R = vector(0, RSIZE);
  list_RplusdR = vector(0, RSIZE);
  list_z = vector(0, ZSIZE);

  Baselen = 0.001 * H;
  LL = 10.0*R200;
  for(i = 1, list_R[0] = 0; i <= RSIZE; i++)
    {
      list_R[i] = exp(log(Baselen) + (i - 1) * (log(LL) - log(Baselen)) / (RSIZE - 1));
    }

  for(i = 1, list_RplusdR[0] = 0; i <= RSIZE; i++)
    {
      list_RplusdR[i] = exp(log(Baselen) + (i - 1 + 0.01) * (log(LL) - log(Baselen)) / (RSIZE - 1));
    }


  for(i = 1, list_z[0] = 0; i <= ZSIZE; i++)
    {
      list_z[i] = exp(log(Baselen) + (i - 1) * (log(LL) - log(Baselen)) / (ZSIZE - 1));
    }



  epi_gamma2 = dvector(0, RSIZE);	/* epicycle gamma^2  */
  epi_kappa2 = dvector(0, RSIZE);


  xl = vector(1, ZSIZE + 1);
  yl = vector(1, ZSIZE + 1);
  D2yl = vector(1, ZSIZE + 1);
}
