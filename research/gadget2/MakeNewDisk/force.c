#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"



void dump_gas_density(void)
{
  FILE *fd;
  int si;
  char gsdfile[50]="";

  strncpy(gsdfile, OutputFile, strlen(OutputFile)-3);
  strcat(gsdfile, ".gsd");
  fd = fopen(gsdfile, "w");
/*fd = fopen("gas_density.txt", "w");*/

  si = RSIZE + 1;
  fwrite(&si, sizeof(int), 1, fd);
  si = ZSIZE + 1;
  fwrite(&si, sizeof(int), 1, fd);

  fwrite(list_R, sizeof(double), RSIZE + 1, fd);
  fwrite(list_z, sizeof(double), ZSIZE + 1, fd);

  fwrite(&RhoGas[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
  fwrite(&CumulMassGas[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);

  fclose(fd);
}


double surface_density_gasdisk(double r)
{

  if(r / H > NumGasScaleLengths)
    return 0;
  else
    return (GasFraction * M_DISK / (2 * M_PI * H * H)) * exp(-r / H)
      / (1 - (1 + NumGasScaleLengths) * exp(-NumGasScaleLengths));
}

void init_central_densities(void)
{
  int tabbin, rbin;
  double sigma;

  if(N_GAS == 0)
    return;

  tabbin = 0;

  for(rbin = RSIZE; rbin >= 0; rbin--)
    {
      RhoGas[rbin][0] = 0;

      sigma = surface_density_gasdisk(list_R[rbin]);

      if(sigma > 0)
	{
	  if(sigma < sigmatab[0] || sigma >= sigmatab[TABSIZE - 1])
	    {
	      if(sigma < sigmatab[0])
		sigma = 0;
	      else
		sigma = 0.999 * sigmatab[TABSIZE - 1];
	    }

	  if(sigma >= sigmatab[0] && sigma < sigmatab[TABSIZE - 1])
	    {
	      while(tabbin < TABSIZE - 1)
		{
		  if(sigma > sigmatab[tabbin + 1])
		    tabbin++;
		  else
		    break;
		}

	      RhoGas[rbin][0] = rhocentral[tabbin];


	      printf("r=%g  RhoGas[rbin][0] = %g  sigma=%g   rbin=%d\n", list_R[rbin], RhoGas[rbin][0], sigma, rbin);
	    }
	}
    }
}


void determine_cumulative_gasdensity(void)
{
  int rbin, zbin;

  if(N_GAS == 0)
    return;

  for(rbin = 0; rbin <= RSIZE; rbin++)
    {
      for(zbin = 0; zbin <= ZSIZE; zbin++)
	CumulMassGas[rbin][zbin] = 0;

      for(zbin = 1; zbin <= ZSIZE; zbin++)
	CumulMassGas[rbin][zbin] = CumulMassGas[rbin][zbin - 1] +
	  0.5 * (RhoGas[rbin][zbin] + RhoGas[rbin][zbin - 1]) * (list_z[zbin] - list_z[zbin - 1]);

      for(zbin = 0; zbin <= ZSIZE; zbin++)
	CumulMassGas[rbin][zbin] *= 2;	/* to account for other side */

      if(surface_density_gasdisk(list_R[rbin]) > 0)
	printf("target=%g  found=%g\n", surface_density_gasdisk(list_R[rbin]), CumulMassGas[rbin][ZSIZE]);
    }
}



void integrate_and_adjust(void)
{
  double rho0, rho, rho2, q, dz, gam, P1;
  double sigmacheck, rhoold;
  double P, P2, drho, target, found;
  double meanweight, u4, z;
  int rbin, zbin, rep;

  if(N_GAS == 0)
    return;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= UnitMass_in_g / UnitEnergy_in_cgs;

  printf("integrate and adjust...\n");

  for(rep = 0; rep < 5; rep++)
    {
      for(rbin = 0; rbin <= RSIZE; rbin++)
	{
	  rho0 = RhoGas[rbin][0];
	  for(zbin = 1; zbin <= ZSIZE; zbin++)
	    RhoGas[rbin][zbin] = 0;

	  if(rho0 > 0)
	    {
	      zbin = 0;
	      z = 0;
	      rho = rho0;
	      q = 0;
	      drho = 0;

	      sigmacheck = 0;

	      while(zbin < ZSIZE && rho > 0.00001 * rho0)
		{
		  dz = list_z[zbin + 1] - list_z[zbin];

		  rhoold = rho;

		  sigmacheck += 0.5 * dz * rho;

		  rho += 0.5 * dz * drho;

		  if(rho > PhysDensThresh)
		    {
		      P = P1 = eqn_of_state(rho);

		      rho2 = 1.1 * rho;
		      P2 = eqn_of_state(rho2);

		      gam = log(P2 / P1) / log(rho2 / rho);
		    }
		  else
		    {
		      P = GAMMA_MINUS1 * rho * u4;
		      gam = 1.0;
		    }

		  drho = -0.5 * (Dphi_z[rbin][zbin] + Dphi_z[rbin][zbin + 1]) * rho * rho / P / gam;

		  rho = rhoold + drho * dz;

		  sigmacheck += 0.5 * rho * dz;

		  zbin++;

		  RhoGas[rbin][zbin] = rho;
		}

	      sigmacheck *= 2;

	    }
	}

      determine_cumulative_gasdensity();

      for(rbin = 0; rbin <= RSIZE; rbin++)
	{
	  if(RhoGas[rbin][0] != 0)
	    {
	      target = surface_density_gasdisk(list_R[rbin]);
	      found = CumulMassGas[rbin][ZSIZE];

	      RhoGas[rbin][0] = RhoGas[rbin][0] * target / found;
	    }
	}


      printf("\n\n");

    }
}





void integrate_gasdensity(void)
{
  double rho0, rho, rho2, q, dz, gam, P1;
  double sigmacheck, rhoold;
  double P, P2, drho;
  double meanweight, u4, z;
  int rbin, zbin;

  if(N_GAS == 0)
    return;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= UnitMass_in_g / UnitEnergy_in_cgs;


  for(rbin = 0; rbin <= RSIZE; rbin++)
    {
      rho0 = RhoGas[rbin][0];
      for(zbin = 1; zbin <= ZSIZE; zbin++)
	RhoGas[rbin][zbin] = 0;

      if(rho0 > 0)
	{
	  zbin = 0;
	  z = 0;
	  rho = rho0;
	  q = 0;
	  drho = 0;

	  sigmacheck = 0;

	  while(zbin < ZSIZE && rho > 0.00001 * rho0)
	    {
	      dz = list_z[zbin + 1] - list_z[zbin];

	      rhoold = rho;

	      sigmacheck += 0.5 * dz * rho;

	      rho += 0.5 * dz * drho;

	      if(rho > PhysDensThresh)
		{
		  P = P1 = eqn_of_state(rho);

		  rho2 = 1.1 * rho;
		  P2 = eqn_of_state(rho2);

		  gam = log(P2 / P1) / log(rho2 / rho);
		}
	      else
		{
		  P = GAMMA_MINUS1 * rho * u4;
		  gam = 1.0;
		}

	      drho = -4 * PI * G * sigmacheck * rho * rho / P / gam;

	      rho = rhoold + drho * dz;

	      sigmacheck += 0.5 * rho * dz;

	      zbin++;

	      RhoGas[rbin][zbin] = rho;
	    }

	  sigmacheck *= 2;

	  printf("rbin=%d sigmacheck=%g\n", rbin, sigmacheck);
	}
    }
}




void set_dummy_particles(void)
{
  int i, j, k, n, rbin, zbin;
  double R, Rold, f, f_, q, z, phi, u;

  /*
     FILE *fd;
   */


  printf("set dummy particles...\n");

  /* first, let's set the gas particles */

  n = 0;

  if(N_GAS > 0)
    {
      for(i = 0; i < RMASSBINS; i++)
	{
	  q = (i + 0.5) / RMASSBINS;

	  R = 1.0;
	  do
	    {
	      f = (-1 + (1 + R) * exp(-R)) / ((1 - (1 + NumGasScaleLengths) * exp(-NumGasScaleLengths))) + q;
	      f_ = -R * exp(-R) / ((1 - (1 + NumGasScaleLengths) * exp(-NumGasScaleLengths)));

	      Rold = R;
	      R = R - f / f_;
	    }
	  while(fabs(R - Rold) / R > 1e-7);

	  R *= H;

	  for(rbin = 0; rbin < RSIZE - 1; rbin++)
	    {
	      if((R >= list_R[rbin] && R < list_R[rbin + 1]))	/* || RhoGas[rbin+1][0]==0) */
		break;
	    }

	  u = (R - list_R[rbin]) / (list_R[rbin + 1] - list_R[rbin]);

	  do
	    {
	      for(zbin = 0; zbin <= ZSIZE; zbin++)
		Zrho[zbin] = (1 - u) * RhoGas[rbin][zbin] + u * RhoGas[rbin + 1][zbin];

	      for(zbin = 1, Zcumul[0] = 0; zbin <= ZSIZE; zbin++)
		Zcumul[zbin] =
		  Zcumul[zbin - 1] + 0.5 * (Zrho[zbin] + Zrho[zbin - 1]) * (list_z[zbin] - list_z[zbin - 1]);

	      /*
	         if(Zcumul[ZSIZE]==0)
	         {
	         printf("Problem!  rbin=%d  i=%d/%d\n", rbin, i, RMASSBINS);
	         exit(0);
	         }
	       */
	      if(Zcumul[ZSIZE] == 0)
		rbin--;
	    }
	  while(Zcumul[ZSIZE] == 0);


	  for(zbin = 0; zbin <= ZSIZE; zbin++)
	    Zcumul[zbin] = Zcumul[zbin] / Zcumul[ZSIZE];

	  for(j = 0; j < ZMASSBINS; j++)
	    {
	      if(j >= ZMASSBINS / 2)
		q = (j - ZMASSBINS / 2 + 0.5) / (ZMASSBINS / 2);
	      else
		q = (j + 0.5) / (ZMASSBINS / 2);

	      for(zbin = 0; zbin < ZSIZE - 1; zbin++)
		{
		  if(q >= Zcumul[zbin] && q < Zcumul[zbin + 1])
		    break;
		}

	      u = (q - Zcumul[zbin]) / (Zcumul[zbin + 1] - Zcumul[zbin]);

	      z = (1 - u) * list_z[zbin] + u * list_z[zbin + 1];

	      if(j >= ZMASSBINS / 2)
		z = -z;


	      for(k = 0; k < PHIMASSBINS; k++)
		{
		  if(n == 1036352)
		    printf("z=%g zbin=%d Zcumul[zbin+1]=%g Zcumul[zbin]=%g  \n", z, zbin, Zcumul[zbin + 1],
			   Zcumul[zbin]);


		  phi = (k + 0.5) / PHIMASSBINS * 2 * M_PI;

		  P[n].Pos[0] = R * cos(phi);
		  P[n].Pos[1] = R * sin(phi);
		  P[n].Pos[2] = z;

		  P[n].Mass = GasFraction * M_DISK / (RMASSBINS * ZMASSBINS * PHIMASSBINS);
		  n++;
		}
	    }
	}

    }

  /* now let's set the star particles (note: n continues to count higher) */

  if(N_DISK > 0)
    {
      for(i = 0; i < RMASSBINS; i++)
	{
	  q = (i + 0.5) / RMASSBINS;

	  R = 1.0;
	  do
	    {
	      f = (1 + R) * exp(-R) + q - 1;
	      f_ = -R * exp(-R);

	      Rold = R;
	      R = R - f / f_;
	    }
	  while(fabs(R - Rold) / R > 1e-6);

	  R *= H;


	  for(j = 0; j < ZMASSBINS; j++)
	    {
	      q = (j + 0.5) / ZMASSBINS;

	      z = Z0 / 2 * log(q / (1 - q));

	      for(k = 0; k < PHIMASSBINS; k++)
		{
		  phi = (k + 0.5) / PHIMASSBINS * 2 * M_PI;

		  P[n].Pos[0] = R * cos(phi);
		  P[n].Pos[1] = R * sin(phi);
		  P[n].Pos[2] = z;


		  P[n].Mass = (1 - GasFraction) * M_DISK / (RMASSBINS * ZMASSBINS * PHIMASSBINS);
		  n++;
		}
	    }
	}

    }


  force_treebuild();

  ErrTolTheta = 0.35;


  /*
     fd = fopen("particles.dat", "w");
     fwrite(&NumPart, sizeof(int), 1, fd);
     for(n = 0; n < NumPart; n++)
     fwrite(&P[n].Pos[0], sizeof(float), 3, fd);
     fclose(fd);
   */
}



void compute_vertical_force_field(void)
{
  int i, j;

  printf("Start computing vertical force field.\n");

  for(i = 0; i <= RSIZE; i++)
    {
      printf("vert force field %d(%d)\n", i, RSIZE);

      for(j = 0; j <= ZSIZE; j++)
	{
	  if(j == 0)
	    Dphi_z[i][j] = 0;
	  else
	    Dphi_z[i][j] = comp_Dphi_z(list_R[i], list_z[j]);
	}

    }
}



void compute_radial_force_field(void)
{
  int i, j;
  double k2, dphi_R_dr;

  printf("Start computing force field.\n");

  for(i = 0; i <= RSIZE; i++)
    {
      printf("force field %d(%d)\n", i, RSIZE);

      for(j = 0; j <= ZSIZE; j++)
	{
	  if(j == 0)
	    {
	      Dphi_z[i][j] = 0;
	      Dphi_z_dR[i][j] = 0;
	    }
	  else
	    {
	      Dphi_z[i][j] = comp_Dphi_z(list_R[i], list_z[j]);
	      Dphi_z_dR[i][j] = comp_Dphi_z(list_RplusdR[i], list_z[j]);
	    }

	  if(i == 0)
	    Dphi_R[i][j] = 0;
	  else
	    Dphi_R[i][j] = comp_Dphi_R(list_R[i], list_z[j]);
	}
    }

  for(i = 1, epi_gamma2[0] = 1; i <= RSIZE; i++)
    {
      dphi_R_dr = comp_Dphi_R(list_RplusdR[i], 0);

      k2 = 3 / list_R[i] * Dphi_R[i][0] + (dphi_R_dr - Dphi_R[i][0]) / (list_RplusdR[i] - list_R[i]);

      epi_gamma2[i] = 4 / list_R[i] * Dphi_R[i][0] / k2;
      epi_kappa2[i] = k2;
    }

  epi_kappa2[0] = epi_kappa2[1];

  printf("Force field finished.\n");
}






double comp_Dphi_z(double R, double z)
{
  return comp_Dphi_Z_disk_tree(R, z) + comp_Dphi_z_halo(R, z) + comp_Dphi_z_bulge(R, z);
}

double comp_Dphi_R(double R, double z)
{
  return comp_Dphi_R_disk_tree(R, z) + comp_Dphi_R_halo(R, z) + comp_Dphi_R_bulge(R, z);
/*
  double Dp1, Dp2, Dp3;

  Dp1= comp_Dphi_R_disk_tree(R, z);
  Dp2= comp_Dphi_R_halo(R, z);
  Dp3= comp_Dphi_R_bulge(R, z);
printf("R= %g  z=%g      Dp1= %g  Dp2= %g  Dp3= %g\n",R,z,Dp1,Dp2,Dp3); fflush(stdout);
  return Dp1 + Dp2 + Dp3;
*/
}
