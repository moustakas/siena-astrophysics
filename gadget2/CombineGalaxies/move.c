#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "globvars.h"

#define WORKSIZE 100000
gsl_integration_workspace *Workspace;





void init_force(struct galaxy_data *g1, struct galaxy_data *g2);
double energy_kernel(double r, void *param);


int compare_kernel(const void *a, const void *b)
{
  if(*((float *) a) < *((float *) b))
    return -1;
  if(*((float *) a) > *((float *) b))
    return +1;
  return 0;
}


void determine_scalefactor(struct galaxy_data *g)
{
#define BINS 10
  float *r;
  int i, j, iter;
  double ri[BINS], frac[BINS];
  double ff, ff_, a, da;

  r = malloc(g->Ndm * sizeof(float));

  for(i = 0; i < g->Ndm; i++)
    {
      j = i + 1 + g->Ngas;
      r[i] = sqrt(g->pos[j][1] * g->pos[j][1] + g->pos[j][2] * g->pos[j][2] + g->pos[j][3] * g->pos[j][3]);
    }

  qsort(r, g->Ndm, sizeof(float), compare_kernel);

  for(i = 0; i < BINS; i++)
    {
      ri[i] = r[(int) (((double) g->Ndm) / (BINS + 1) * (i + 1))];
      frac[i] = ((double) (i + 1)) / (BINS + 1);
    }

  free(r);

  /*now newton iteration */

  a = ri[(int) 0.1 * g->Ndm];
  iter = 0;
  do
    {
      for(i = 0, ff = ff_ = 0; i < BINS; i++)
	{
	  ff += 4 * ri[i] * ri[i] * (frac[i] - ri[i] * ri[i] / pow(a + ri[i], 2)) / pow(a + ri[i], 3);
	  ff_ += 8 * pow(ri[i], 4) / pow(a + ri[i], 6) -
	    12 * ri[i] * ri[i] * (frac[i] - ri[i] * ri[i] / pow(a + ri[i], 2)) / pow(a + ri[i], 4);
	}

      da = -ff / ff_;

      if(abs(da) > 0.05 * a)
	da = abs(da) / da * 0.05 * a;

      a += da;
      iter++;
    }
  while(fabs(da / a) > 1.0e-3 && iter < 100);

  if(iter >= 100)
    {
      printf("iteration failed\n");
      exit(0);
    }

  g->Scalefac = a;
}

void get_force(struct galaxy_data *g1, struct galaxy_data *g2, double *pos, double *acc)
{
  gsl_function F;
  double M1, M2, Mreduced, Mb1, Mb2, Mdm1, Mdm2, dist, frc[3];
  double result, abserr, e1, e2, r1, r2;
  double params[3];
  int i;

  M1 = g1->Mtot;
  M2 = g2->Mtot;
  Mreduced = M1 * M2 / (M1 + M2);

  dist = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

#ifdef ORBITCORRECTION
  Mb1 = g1->Mtot - g1->Mdm;
  Mb2 = g2->Mtot - g2->Mdm;
  Mdm1 = g1->Mdm;
  Mdm2 = g2->Mdm;

  for(i = 0; i < 3; i++)
    acc[i] = -G * Mb1 * Mb2 / Mreduced * pos[i] / (dist * dist * dist);

  for(i = 0; i < 3; i++)
    acc[i] += -G * Mb1 * Mdm2 / Mreduced * pos[i] / (dist * (dist + g2->Scalefac) * (dist + g2->Scalefac));

  for(i = 0; i < 3; i++)
    acc[i] += -G * Mb2 * Mdm1 / Mreduced * pos[i] / (dist * (dist + g1->Scalefac) * (dist + g1->Scalefac));

  F.function = &energy_kernel;
  F.params = params;

  params[0] = dist;
  params[1] = g1->Scalefac;
  params[2] = g2->Scalefac;

  gsl_integration_qagiu(&F, 0, 0.000001, 1.0e-8, WORKSIZE, Workspace, &result, &abserr);
  r1 = params[0];
  e1 = result / r1;

  params[0] = 1.0001 * dist;
  gsl_integration_qagiu(&F, 0, 0.000001, 1.0e-8, WORKSIZE, Workspace, &result, &abserr);
  r2 = params[0];
  e2 = result / r2;

  for(i = 0; i < 3; i++)
    frc[i] = 0.5 * G * Mdm1 * Mdm2 * params[1] * (e2 - e1) / (r2 - r1) * pos[i] / dist;

  params[0] = dist;
  params[1] = g2->Scalefac;
  params[2] = g1->Scalefac;

  gsl_integration_qagiu(&F, 0, 0.0001, 1.0e-8, WORKSIZE, Workspace, &result, &abserr);
  r1 = params[0];
  e1 = result / r1;

  params[0] = 1.0001 * dist;
  gsl_integration_qagiu(&F, 0, 0.0001, 1.0e-8, WORKSIZE, Workspace, &result, &abserr);
  r2 = params[0];
  e2 = result / r2;

  for(i = 0; i < 3; i++)
    frc[i] += 0.5 * G * Mdm1 * Mdm2 * params[1] * (e2 - e1) / (r2 - r1) * pos[i] / dist;

  for(i = 0; i < 3; i++)
    acc[i] += frc[i] / Mreduced;

#else
  for(i = 0; i < 3; i++)
    acc[i] = -G * M1 * M2 / Mreduced * pos[i] / (dist * dist * dist);

#endif
}

void correct_orbit(struct galaxy_data *g1, struct galaxy_data *g2, double *pos, double *vel, double rtarget)
{
  double dist, dt;
  double acc[3];
  int i;

  Workspace = gsl_integration_workspace_alloc(WORKSIZE);

  do
    {
      dist = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

      dt = 0.0025 / sqrt(G * (g1->Mtot + g2->Mtot) / pow(dist, 3));

      get_force(g1, g2, pos, acc);

      for(i = 0; i < 3; i++)
	vel[i] += 0.5 * acc[i] * dt;

      for(i = 0; i < 3; i++)
	pos[i] += vel[i] * dt;

      dist = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

      get_force(g1, g2, pos, acc);

      for(i = 0; i < 3; i++)
	vel[i] += 0.5 * acc[i] * dt;

      /*  printf("%g\n", dist); */
    }
  while(dist > rtarget);
}

double energy_kernel(double r, void *param)
{
  double r0, a1, a2, *p;

  p = (double *) param;
  r0 = p[0];
  a1 = p[1];
  a2 = p[2];

  if(r > r0)
    return pow(r + a1, -3) * (2 * r0 + a2 * log((a2 + r - r0) / (a2 + r + r0)));
  else
    return pow(r + a1, -3) * (2 * r + a2 * log((a2 - r + r0) / (a2 + r + r0)));
}



void move_galaxies(struct galaxy_data *g1, struct galaxy_data *g2, double rmin, double rstart, double ecc)
{
  int i;
  double p, phistart, rdotstart, phidotstart;
  double M, Rstart;
  double xyz[3], pos[3], vel[3];
#ifdef ECC_NONZERO
  double rmax=0;
#endif

  void translate(struct galaxy_data *g, double vec[3]);
  void vel_translate(struct galaxy_data *g, double vec[3]);

  Rstart = rstart;

#ifdef ORBITCORRECTION
  determine_scalefactor(g1);
  printf("Galaxy 1 has Scalefactor = %g\n", g1->Scalefac);
  determine_scalefactor(g2);
  printf("Galaxy 2 has Scalefactor = %g\n", g2->Scalefac);

  Rstart = (g1->Scalefac + g2->Scalefac) * sqrt(0.99) / (1 - sqrt(0.99));
#endif

  p = 2 * rmin;

  phistart = -acos(p / Rstart - 1);	/* choose neg. sign here */

#if ECC_NONZERO
  printf("rstart= %g\n",rstart);
  printf("ecc= %g\n",ecc);
  p = (1+ecc) * rmin;
  if(ecc>0)
    {
      phistart = -acos((p / Rstart - 1)/ecc);	/* choose neg. sign here */
    }
#endif

  M = g1->Mtot + g2->Mtot;

  phidotstart = sqrt(G * M * p) / (Rstart * Rstart);

  if(rmin > 0)
    {
#ifndef ECC_NONZERO
      rdotstart = Rstart * Rstart / p * sin(phistart) * phidotstart;
#else
      rdotstart = ecc * Rstart * Rstart / p * sin(phistart) * phidotstart;
#endif
    }
  else
    {
      rdotstart = -sqrt(2 * G * M / Rstart);
    }

  printf("rdotstart = %g\n",rdotstart);
  printf("phidotstart = %g\n",phidotstart);
  printf("r*phidotstart = %g\n",rstart*phidotstart);

#ifdef ECC_NONZERO
  if(ecc<1) rmax= p/(1-ecc);
  printf("Rstart=%g\nrmax=%g\t",Rstart,rmax);
  if(Rstart>rmax)
    {
      printf("  ..... this might be a problem, maybe we should stop ....\n\n");
      /* exit(0); */
    }
  else
      printf("\n");
#endif

  pos[0] = Rstart * cos(phistart);
  pos[1] = Rstart * sin(phistart);
  pos[2] = 0;

  vel[0] = rdotstart * cos(phistart) - Rstart * phidotstart * sin(phistart);
  vel[1] = rdotstart * sin(phistart) + Rstart * phidotstart * cos(phistart);
  vel[2] = 0;

  printf("Relative V= %g\n",sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]));

#ifdef ORBITCORRECTION
  printf("correcting orbit for dark halo potential...\n");
  correct_orbit(g1, g2, pos, vel, rstart);
#endif

  printf("Relative V (post OC)= %g\n",sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]));

  for(i = 0; i <= 2; i++)
    xyz[i] = pos[i] * g2->Mtot / M;
  printf("galaxy 1 moved to: ( %g | %g | %g )\n", xyz[0], xyz[1], xyz[2]);
  translate(g1, &xyz[0]);

  for(i = 0; i <= 2; i++)
    xyz[i] = pos[i] * (-g1->Mtot / M);
  translate(g2, &xyz[0]);
  printf("galaxy 2 moved to: ( %g | %g | %g )\n", xyz[0], xyz[1], xyz[2]);

  for(i = 0; i <= 2; i++)
    xyz[i] = vel[i] * g2->Mtot / M;
  vel_translate(g1, &xyz[0]);
  printf("velocity galaxy 1: ( %g | %g | %g )\n", xyz[0], xyz[1], xyz[2]);

  for(i = 0; i <= 2; i++)
    xyz[i] = vel[i] * (-g1->Mtot / M);
  vel_translate(g2, &xyz[0]);
  printf("velocity galaxy 2: ( %g | %g | %g )\n", xyz[0], xyz[1], xyz[2]);

  printf("\n");
}




void translate(struct galaxy_data *g, double vec[3])
{
  int i, j;


  for(i = 1; i <= g->Ntot; i++)
    for(j = 0; j <= 2; j++)
      g->pos[i][j + 1] += vec[j];
}


void vel_translate(struct galaxy_data *g, double vec[3])
{
  int i, j;


  for(i = 1; i <= g->Ntot; i++)
    for(j = 0; j <= 2; j++)
      g->vel[i][j + 1] += vec[j];
}
