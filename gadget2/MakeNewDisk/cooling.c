#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "prototypes.h"
#include "globvars.h"
#include "cooling.h"


#define NCOOLTAB  2000

#define SMALLNUM 1.0e-60
#define COOLLIM  0.1
#define HEATLIM	 20.0

#define MAXITER 200


static double XH = HYDROGEN_MASSFRAC;	/* hydrogen abundance by mass */
static double yhelium;

static double mhboltz;		/* hydrogen mass over Boltzmann constant */
static double ethmin;		/* minimum internal energy for neutral gas */

static double Tmin = 0.0;	/* in log10 */
static double Tmax = 9.0;
static double deltaT;

static double *BetaH0, *BetaHep, *Betaff;
static double *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
static double *GammaeH0, *GammaeHe0, *GammaeHep;

static double J_UV = 0, gJH0 = 0, gJHep = 0, gJHe0 = 0, epsH0 = 0, epsHep = 0, epsHe0 = 0;

static double ne, necgs, nHcgs;
static double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
static double gJH0ne, gJHe0ne, gJHepne;
static double nH0, nHp, nHep, nHe0, nHepp;



static double DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;

/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */
double DoCooling(double u_old, double rho, double dt, double *ne_guess)
{
  double u, du;
  double u_lower, u_upper;
  double ratefact;
  double LambdaNet;
  int iter = 0;

  DoCool_u_old_input = u_old;
  DoCool_rho_input = rho;
  DoCool_dt_input = dt;
  DoCool_ne_guess_input = *ne_guess;


  rho *= UnitDensity_in_cgs * HubbleParam * HubbleParam;	/* convert to physical cgs units */
  u_old *= UnitPressure_in_cgs / UnitDensity_in_cgs;
  dt *= UnitTime_in_s / HubbleParam;

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  u = u_old;
  u_lower = u;
  u_upper = u;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess);

  /* bracketing */

  if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess) * dt < 0)
	{
	  u_upper *= 1.1;
	  u_lower *= 1.1;
	}

    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess) * dt > 0)
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
	}
    }

  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = CoolingRateFromU(u, rho, ne_guess);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
	{
	  u_upper = u;
	}
      else
	{
	  u_lower = u;
	}

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("u= %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      exit(10);
    }

  u *= UnitDensity_in_cgs / UnitPressure_in_cgs;	/* to internal units */

  return u;
}



/* returns cooling time. 
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime(double u_old, double rho, double *ne_guess)
{
  double u;
  double ratefact;
  double LambdaNet, coolingtime;

  DoCool_u_old_input = u_old;
  DoCool_rho_input = rho;
  DoCool_ne_guess_input = *ne_guess;

  rho *= UnitDensity_in_cgs * HubbleParam * HubbleParam;	/* convert to physical cgs units */
  u_old *= UnitPressure_in_cgs / UnitDensity_in_cgs;


  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  u = u_old;


  LambdaNet = CoolingRateFromU(u, rho, ne_guess);

  /* bracketing */

  if(LambdaNet >= 0)		/* ups, we have actually heating due to UV background */
    return 0;

  coolingtime = u_old / (-ratefact * LambdaNet);

  coolingtime *= HubbleParam / UnitTime_in_s;

  return coolingtime;
}


/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess)
{
  double m, dm;
  double m_lower, m_upper;
  double ratefact;
  double LambdaNet;
  int iter = 0;

  DoCool_u_old_input = u;
  DoCool_rho_input = rho;
  DoCool_dt_input = dt;
  DoCool_ne_guess_input = *ne_guess;

  if(fac <= 0)			/* the hot phase is actually colder than the cold reservoir! */
    {
      return 0.01 * m_old;
    }

  rho *= UnitDensity_in_cgs * HubbleParam * HubbleParam;	/* convert to physical cgs units */
  u *= UnitPressure_in_cgs / UnitDensity_in_cgs;
  dt *= UnitTime_in_s / HubbleParam;
  fac *= UnitMass_in_g / UnitEnergy_in_cgs;

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho * fac;

  m = m_old;
  m_lower = m;
  m_upper = m;

  LambdaNet = CoolingRateFromU(u, rho, ne_guess);

  /* bracketing */

  if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
    {
      m_upper *= sqrt(1.1);
      m_lower /= sqrt(1.1);
      while(m_upper - m_old -
	    m_upper * m_upper / m_old * ratefact * CoolingRateFromU(u, rho * m_upper / m_old,
								    ne_guess) * dt < 0)
	{
	  m_upper *= 1.1;
	  m_lower *= 1.1;
	}
    }

  if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
    {
      m_lower /= sqrt(1.1);
      m_upper *= sqrt(1.1);
      while(m_lower - m_old -
	    m_lower * m_lower / m_old * ratefact * CoolingRateFromU(u, rho * m_lower / m_old,
								    ne_guess) * dt > 0)
	{
	  m_upper /= 1.1;
	  m_lower /= 1.1;
	}
    }

  do
    {
      m = 0.5 * (m_lower + m_upper);

      LambdaNet = CoolingRateFromU(u, rho * m / m_old, ne_guess);

      if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0)
	{
	  m_upper = m;
	}
      else
	{
	  m_lower = m;
	}

      dm = m_upper - m_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("m= %g\n", m);
    }
  while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      printf("m_old= %g\n", m_old);
      exit(11);
    }

  return m;
}








void cool_test(void)
{
  double uin, rhoin, tempin, muin, nein;

  uin = 6.01329e+09;
  rhoin = 7.85767e-29;
  tempin = 34.0025;
  muin = 0.691955;

  nein = (1 + 4 * yhelium) / muin - (1 + yhelium);

  printf("%g\n", convert_u_to_temp(uin, rhoin, &nein));
}


/* this function determines the electron fraction, and hence the mean 
 * molecular weight. With it arrives at a self-consistent temperature.
 * Element abundances and the rates for the emission are also computed
 */
double convert_u_to_temp(double u, double rho, double *ne_guess)
{
  double temp, temp_old, temp_new, max = 0, ne_old;
  double mu;
  int iter = 0;

  double u_input, rho_input, ne_input;

  u_input = u;
  rho_input = rho;
  ne_input = *ne_guess;

  mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  do
    {
      ne_old = *ne_guess;

      find_abundances_and_rates(log10(temp), rho, ne_guess);
      temp_old = temp;

      mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);

      temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max =
	dmax(max,
	     temp_new / (1 + yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;

      if(iter > (MAXITER - 10))
	printf("-> temp= %g ne=%g\n", temp, *ne_guess);
    }
  while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in convert_u_to_temp()\n");
      printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);

      exit(12);
    }

  return temp;
}


/* this function computes the actual abundance ratios 
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess)
{
  double neold, nenew;
  int j, niter;
  double Tlow, Thi, flow, fhi, t;

  double logT_input, rho_input, ne_input;

  logT_input = logT;
  rho_input = rho;
  ne_input = *ne_guess;

  if(logT <= Tmin)		/* everything neutral */
    {
      nH0 = 1.0;
      nHe0 = yhelium;
      nHp = 0;
      nHep = 0;
      nHepp = 0;
      ne = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Tmax)		/* everything is ionized */
    {
      nH0 = 0;
      nHe0 = 0;
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *ne_guess = ne;		/* note: in units of the hydrogen number density */
      return;
    }

  t = (logT - Tmin) / deltaT;
  j = (int) t;
  Tlow = Tmin + deltaT * j;
  Thi = Tlow + deltaT;
  fhi = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  ne = *ne_guess;
  neold = ne;
  niter = 0;
  necgs = ne * nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;

      aHp = flow * AlphaHp[j] + fhi * AlphaHp[j + 1];
      aHep = flow * AlphaHep[j] + fhi * AlphaHep[j + 1];
      aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
      ad = flow * Alphad[j] + fhi * Alphad[j + 1];
      geH0 = flow * GammaeH0[j] + fhi * GammaeH0[j + 1];
      geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
      geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];

      if(necgs <= 1.e-25 || J_UV == 0)
	{
	  gJH0ne = gJHe0ne = gJHepne = 0;
	}
      else
	{
	  gJH0ne = gJH0 / necgs;
	  gJHe0ne = gJHe0 / necgs;
	  gJHepne = gJHep / necgs;
	}

      nH0 = aHp / (aHp + geH0 + gJH0ne);	/* eqn (33) */
      nHp = 1.0 - nH0;		/* eqn (34) */

      if((gJHe0ne + geHe0) <= SMALLNUM)	/* no ionization at all */
	{
	  nHep = 0.0;
	  nHepp = 0.0;
	  nHe0 = yhelium;
	}
      else
	{
	  nHep = yhelium / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
	  nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
	  nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
	}

      neold = ne;

      ne = nHp + nHep + 2 * nHepp;	/* eqn (38) */
      necgs = ne * nHcgs;

      if(J_UV == 0)
	break;

      nenew = 0.5 * (ne + neold);
      ne = nenew;
      necgs = ne * nHcgs;

      if(fabs(ne - neold) < 1.0e-4)
	break;

      if(niter > (MAXITER - 10))
	printf("ne= %g  niter=%d\n", ne, niter);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("no convergence reached in find_abundances_and_rates()\n");
      printf("logT_input= %g  rho_input= %g  ne_input= %g\n", logT_input, rho_input, ne_input);
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      exit(13);
    }

  bH0 = flow * BetaH0[j] + fhi * BetaH0[j + 1];
  bHep = flow * BetaHep[j] + fhi * BetaHep[j + 1];
  bff = flow * Betaff[j] + fhi * Betaff[j + 1];

  *ne_guess = ne;
}



/*  this function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates 
 *  (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRateFromU(double u, double rho, double *ne_guess)
{
  double temp;

  temp = convert_u_to_temp(u, rho, ne_guess);

  return CoolingRate(log10(temp), rho, ne_guess);
}


/*  this function computes the self-consistent temperature
 *  and abundance ratios 
 */
double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer)
{
  double temp;

  rho *= UnitDensity_in_cgs * HubbleParam * HubbleParam;	/* convert to physical cgs units */
  u *= UnitPressure_in_cgs / UnitDensity_in_cgs;

  temp = convert_u_to_temp(u, rho, ne_guess);

  *nH0_pointer = nH0;

  return temp;
}




extern FILE *fd;

/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRate(double logT, double rho, double *nelec)
{
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double T;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT;	/* floor at Tmin */


  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */


  if(logT < Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec);

      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */

      T = pow(10.0, logT);

      LambdaExcH0 = bH0 * ne * nH0;
      LambdaExcHep = bHep * ne * nHep;
      LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */

      LambdaIonH0 = 2.18e-11 * geH0 * ne * nH0;
      LambdaIonHe0 = 3.94e-11 * geHe0 * ne * nHe0;
      LambdaIonHep = 8.72e-11 * geHep * ne * nHep;
      LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */

      LambdaRecHp = 1.036e-16 * T * ne * (aHp * nHp);
      LambdaRecHep = 1.036e-16 * T * ne * (aHep * nHep);
      LambdaRecHepp = 1.036e-16 * T * ne * (aHepp * nHepp);
      LambdaRecHepd = 6.526e-11 * ad * ne * nHep;
      LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;

      LambdaFF = bff * (nHp + nHep + 4 * nHepp) * ne;

      Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

      LambdaCmptn = 0;

      Heat = 0;
      if(J_UV != 0)
	Heat += (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs;
    }
  else				/* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present.  
         Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
	LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

      /* very hot: H and He both fully ionized */
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *nelec = ne;		/* note: in units of the hydrogen number density */

      T = pow(10.0, logT);
      LambdaFF =
	1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * ne;

      LambdaCmptn = 0;

      Lambda = LambdaFF + LambdaCmptn;
    }

  return (Heat - Lambda);
}





double INLINE_FUNC LogTemp(double u, double ne)	/* ne= electron density in terms of hydrogen density */
{
  double T;

  if(u < ethmin)
    u = ethmin;

  T = log10(GAMMA_MINUS1 * u * mhboltz * (1 + 4 * yhelium) / (1 + ne + yhelium));

  return T;
}


void *mymalloc(size_t size)
{
  void *ptr;

  ptr = malloc(size);

  if(!(ptr))
    {
      printf("failed to allocate %d bytes of memory.\n", (int) size);
      exit(14);
    }

  return ptr;
}

void InitCoolMemory(void)
{
  BetaH0 = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  BetaHep = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  AlphaHp = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  AlphaHep = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  Alphad = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  AlphaHepp = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  GammaeH0 = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  GammaeHe0 = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  GammaeHep = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
  Betaff = (double *) mymalloc((NCOOLTAB + 1) * sizeof(double));
}


void MakeCoolingTable(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19 */
{
  int i;
  double T;
  double Tfact;

  yhelium = (1 - XH) / (4 * XH);

  mhboltz = PROTONMASS / BOLTZMANN;

  MinGasTemp= 20.0;

  Tmin = log10(0.1 * MinGasTemp);
  

  deltaT = (Tmax - Tmin) / NCOOLTAB;

  ethmin = pow(10.0, Tmin) * (1. + yhelium) / ((1. + 4. * yhelium) * mhboltz * GAMMA_MINUS1);
  /* minimum internal energy for neutral gas */

  for(i = 0; i <= NCOOLTAB; i++)
    {
      BetaH0[i] =
	BetaHep[i] =
	Betaff[i] =
	AlphaHp[i] = AlphaHep[i] = AlphaHepp[i] = Alphad[i] = GammaeH0[i] = GammaeHe0[i] = GammaeHep[i] = 0;

      T = pow(10.0, Tmin + deltaT * i);

      Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

      if(118348 / T < 70)
	BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;

      if(473638 / T < 70)
	BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;

      Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));

      AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);

      AlphaHep[i] = 1.5e-10 * pow(T, -0.6353);

      AlphaHepp[i] = 4. * AlphaHp[i];

      if(470000 / T < 70)
	Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));

      if(157809.1 / T < 70)
	GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;

      if(285335.4 / T < 70)
	GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;

      if(631515.0 / T < 70)
	GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
    }


}





/* table input (from file TREECOOL) for ionizing parameters */

#define JAMPL	1.0		/* amplitude factor relative to input table */
#define TABLESIZE 200		/* Max # of lines in TREECOOL */

static float inlogz[TABLESIZE];
static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
static int nheattab;		/* length of table */


void ReadIonizeParams(char *fname)
{
  int i;
  FILE *fdcool;

  if(!(fdcool = fopen(fname, "r")))
    {
      printf(" Cannot read ionization table in file `%s'\n", fname);
      exit(456);
    }

  for(i = 0; i < TABLESIZE; i++)
    gH0[i] = 0;

  for(i = 0; i < TABLESIZE; i++)
    if(fscanf(fdcool, "%g %g %g %g %g %g %g",
	      &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF)
      break;

  fclose(fdcool);

  /*  nheattab is the number of entries in the table */

  for(i = 0, nheattab = 0; i < TABLESIZE; i++)
    if(gH0[i] != 0.0)
      nheattab++;
    else
      break;


  printf("\n\nread ionization table with %d entries in file `%s'.\n\n", nheattab, fname);
}


void IonizeParams(void)
{
#ifdef IONIZETABLE
  IonizeParamsTable();
#else
  IonizeParamsFunction();
#endif
}



void IonizeParamsTable(void)
{

  gJHe0 = gJHep = gJH0 = 0;
  epsHe0 = epsHep = epsH0 = 0;
  J_UV = 0;
  return;


}


void SetZeroIonization(void)
{
  gJHe0 = gJHep = gJH0 = 0;
  epsHe0 = epsHep = epsH0 = 0;
  J_UV = 0;
}


void IonizeParamsFunction(void)
{

  J_UV = 0.;
  gJHe0 = gJHep = gJH0 = 0.;
  epsHe0 = epsHep = epsH0 = 0.;

}



void InitCool(void)
{
  InitCoolMemory();
  MakeCoolingTable();

#ifdef IONIZETABLE
  ReadIonizeParams("TREECOOL");
#endif


  IonizeParams();
}

