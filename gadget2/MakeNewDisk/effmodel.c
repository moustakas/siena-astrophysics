#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"
#include "cooling.h"



void set_sfr_parameters(void)
{
  double meanweight, feedbackenergyinergs;


  PhysDensThresh = 0;


  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * TempClouds;
  EgySpecCold *= UnitMass_in_g / UnitEnergy_in_cgs;
  printf("EgySpecCold= %g\n", EgySpecCold);

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  EgySpecSN = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * TempSupernova;
  EgySpecSN *= UnitMass_in_g / UnitEnergy_in_cgs;
  printf("EgySpecSN= %g\n", EgySpecSN);

  FeedbackEnergy = FactorSN / (1 - FactorSN) * EgySpecSN;

  feedbackenergyinergs = FeedbackEnergy / UnitMass_in_g * (UnitEnergy_in_cgs * SOLAR_MASS);


  printf("Feedback energy per formed solar mass in stars= %g  ergs\n", feedbackenergyinergs);
}

void init_clouds(void)
{
  double A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;


  set_sfr_parameters();
  InitCool();
  SetZeroIonization();


  A0 = FactorEVP;

  egyhot = EgySpecSN / A0;
  printf("egyhot= %g\n", egyhot);

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  printf("meanweight= %g\n",meanweight);

  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  printf("10^4 K\n");
  printf("u4= %g (cgs = ergs/g) \n", egyhot);
  u4 *= UnitMass_in_g / UnitEnergy_in_cgs;
  printf("u4= %g (gu)\n", egyhot);

  printf("u4=%g\n", u4);
  U4 = u4;


  dens = 1.0e6 * 3 * 0.1 * 0.1 / (8 * PI * G);
  printf("dens= %g\n", dens);

  ne = 1.0;
  printf("ne= %g\n",ne);
  SetZeroIonization();
  printf("ne= %g\n",ne);
  tcool = GetCoolingTime(egyhot, dens, &ne);
  printf("ne= %g\n",ne);
  printf("tcool= %g\n", tcool);

  coolrate = egyhot / tcool / dens;
  printf("coolrate= %g\n", coolrate);

  x = (egyhot - u4) / (egyhot - EgySpecCold);
  printf("x= %g\n",x);

  PhysDensThresh =
    x / pow(1 - x, 2) * (FactorSN * EgySpecSN - (1 - FactorSN) * EgySpecCold) / (MaxSfrTimescale * coolrate);


  printf("\nA0= %g  \n", A0);
  printf("Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n", PhysDensThresh,
	 PhysDensThresh / (PROTONMASS / HYDROGEN_MASSFRAC / UnitDensity_in_cgs));
  printf("EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n", x);
  printf("tcool=%g dens=%g egyhot=%g\n", tcool, dens, egyhot);


  if (M_GAS > 0)
	dump_eqn_of_state();
  integrate_surfacedensity();
}


void dump_eqn_of_state(void)
{
  double rho;
  FILE *fd;
  char eqsfile[50]="";

  strncpy(eqsfile, OutputFile, strlen(OutputFile)-4);
  strcat(eqsfile, ".eqs");
  if((fd = fopen(eqsfile, "w")))
    {
	for(rho = PhysDensThresh / 100; rho < PhysDensThresh * 1000000; rho *= 1.1)
	  {
		fprintf(fd, "%g %g\n", rho, eqn_of_state(rho));
	  }
	fclose(fd);
    }
  else
    {
	fprintf(stderr, "Can't open file %s\n",eqsfile);
    }
}



double eqn_of_state(double rho)
{
  double tsfr, factorEVP, egyhot, ne, y, x, egyeff, P, P4, Pmulti, tcool;
  double u4, meanweight;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= UnitMass_in_g / UnitEnergy_in_cgs;

  P4 = GAMMA_MINUS1 * rho * u4;


  if(rho >= PhysDensThresh)
    {

      tsfr = sqrt(PhysDensThresh / rho) * MaxSfrTimescale;

      factorEVP = pow(rho / PhysDensThresh, -0.8) * FactorEVP;

      egyhot = EgySpecSN / (1 + factorEVP) + EgySpecCold;

      ne = 1.0;
      tcool = GetCoolingTime(egyhot, rho, &ne);

      y = tsfr / tcool * egyhot / (FactorSN * EgySpecSN - (1 - FactorSN) * EgySpecCold);
      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + EgySpecCold * x;

      Pmulti = GAMMA_MINUS1 * rho * egyeff;

      P = FactorForSofterEQS * Pmulti + (1 - FactorForSofterEQS) * P4;
    }
  else
    {
      P = P4;
    }

  return P;
}





void integrate_surfacedensity(void)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, P1;
  double sigmacheck, rhoold;
  double x = 0, P, P2, drho, dq;
  double meanweight, u4, z, tsfr;
  double rhomin, rhomax;
  int bin, rep;
  FILE *fd;
  char sdfile[50]="";

  if(N_GAS == 0)
    return;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= UnitMass_in_g / UnitEnergy_in_cgs;


  NumGasScaleLengths = -log((GAMMA - 1) * u4 * H / G / MaxGasDiskHeight / (GasFraction * M_DISK));

  if(NumGasScaleLengths < 8.0)
    NumGasScaleLengths = 8.0;

  printf("NumGasScaleLengths = %g\n", NumGasScaleLengths);

  printf("integrating surface density...\n");

  strncpy(sdfile, OutputFile, strlen(OutputFile)-4);
  strcat(sdfile, ".sd");
  fd = fopen(sdfile, "w");
/*fd = fopen("surface_density.txt", "w");*/


  rhomin = PhysDensThresh / 100;
  rhomax = 10000 * PhysDensThresh;


  for(rep = 0; rep < TABSIZE; rep++)
    rhocentral[rep] = exp(log(rhomin) + rep * (log(rhomax) - log(rhomin)) / (TABSIZE - 1));

  for(rep = 0; rep < TABSIZE; rep++)
    {
      rho0 = rhocentral[rep];

      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.0025;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.00001 * rho0)
	{
	  if(rho > PhysDensThresh)
	    {
	      tsfr = sqrt(PhysDensThresh / rho) * MaxSfrTimescale;

	      P = P1 = eqn_of_state(rho);

	      rho2 = 1.1 * rho;
	      P2 = eqn_of_state(rho2);

	      gam = log(P2 / P1) / log(rho2 / rho);
	    }
	  else
	    {
	      tsfr = 0;

	      P = GAMMA_MINUS1 * rho * u4;
	      gam = 1.0;

	      sigma_u4 += rho * dz;
	    }

	  drho = q;
	  dq = -(gam - 2) / rho * q * q - 4 * PI * G / (gam * P) * rho * rho * rho;

	  sigma += rho * dz;
	  if(tsfr > 0)
	    {
	      sigmasfr += (1 - FactorSN) * rho * x / tsfr * dz;
	    }

	  rho += drho * dz;
	  q += dq * dz;
	}


      sigma *= 2;		/* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;


      bin = 0;
      z = 0;
      rho = rho0;
      q = 0;
      drho = 0;

      sigmacheck = 0;

      while(bin < ZSIZE && rho > 0.00001 * rho0)
	{
	  dz = list_z[bin + 1] - list_z[bin];

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

	  bin++;
	}

      sigmacheck *= 2;

      fprintf(fd, "%g %g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4, sigmacheck);

      sigmatab[rep] = sigma;
    }

  fclose(fd);

  printf("done.\n");
}










/* returns the maximum of two double
 */
double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}

/* returns the minimum of two double
 */
double dmin(double x, double y)
{
  if(x < y)
    return x;
  else
    return y;
}



/* returns the maximum of two integers
 */
int imax(int x, int y)
{
  if(x > y)
    return x;
  else
    return y;
}

/* returns the minimum of two integers
 */
int imin(int x, int y)
{
  if(x < y)
    return x;
  else
    return y;
}
