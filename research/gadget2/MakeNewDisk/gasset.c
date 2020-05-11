#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"


static double R, z;



void set_gas_velocities(void)
{
  int i;
  long dum;
  int iz, ir;
  double ur, uz;
  double vstream_phi;
  double vr, vphi;
  double vx, vy, vz;



  if(N_GAS == 0)
    return;

  dum = drand48() * 1e8;

  printf("set gas velocities...");
  fflush(stdout);

  for(i = 1; i <= N_GAS; i++)
    {
      R = sqrt(xp_gas[i] * xp_gas[i] + yp_gas[i] * yp_gas[i]);
      z = fabs(zp_gas[i]);

      if(R < Baselen)
	ir = 0;
      else
	ir = (log(R) - log(Baselen)) / (log(LL) - log(Baselen)) * (RSIZE - 1) + 1;

      ur = (R - list_R[ir]) / (list_R[ir + 1] - list_R[ir]);

      if(z < Baselen)
	iz = 0;
      else
	iz = (log(z) - log(Baselen)) / (log(LL) - log(Baselen)) * (ZSIZE - 1) + 1;

      uz = (z - list_z[iz]) / (list_z[iz + 1] - list_z[iz]);

      if(ir < 0 || ir >= RSIZE)
	{
	  printf("ir=%d out of range\n", ir);
	}

      if(iz < 0 || iz >= ZSIZE)
	{
	  printf("iz=%d out of range\n", iz);
	}


      vstream_phi = VelStreamGas[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamGas[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamGas[ir][iz + 1] * (1 - ur) * (uz) + VelStreamGas[ir + 1][iz + 1] * (ur) * (uz);



      vr = 0;
      vz = 0;

      vphi = vstream_phi;

      vx = vr * xp_gas[i] / R - vphi * yp_gas[i] / R;
      vy = vr * yp_gas[i] / R + vphi * xp_gas[i] / R;

      vxp_gas[i] = vx;
      vyp_gas[i] = vy;
      vzp_gas[i] = vz;

      /*
      u_gas[i] = U4/100.0;
      */
      u_gas[i] = U4;

    }

  printf("done.\n");
  fflush(stdout);
}





void set_gas_positions(void)
{
  int i, n_disk, n_HI;
  double q, R, f, f_, Rold, phi, u;
  int rbin, zbin;
  int iter;

  if(N_GAS == 0)
    return;

  printf("set gas positions...\n");

  if(HI_GasMassFraction > 0)
    {
      n_disk = (1 - HI_GasMassFraction) * N_GAS;
      n_HI = N_GAS - n_disk;
    }
  else
    {
      n_disk = N_GAS;
      n_HI = 0;
    }

  for(i = 1; i <= n_disk; i++)
    {
      do
        q = drand48();
      while(q>0.995);

      R = 1.0;
      iter=0;
      do
	{
	  f = (-1 + (1 + R) * exp(-R)) / ((1 - (1 + NumGasScaleLengths) * exp(-NumGasScaleLengths))) + q;
	  f_ = -R * exp(-R) / ((1 - (1 + NumGasScaleLengths) * exp(-NumGasScaleLengths)));

	  Rold = R;
	  R = R - f / f_;
	  if(iter > 200)
	    printf("R=%g q=%g rold=%g iter=%d\n", R, q, Rold, iter);
	  if(iter > 210)
	    break;
	  iter++;
	}
      while(fabs(R - Rold) / R > 1e-5);

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
	    Zcumul[zbin] = Zcumul[zbin - 1] + 0.5 * (Zrho[zbin] + Zrho[zbin - 1]) * (list_z[zbin] - list_z[zbin - 1]);

	  if(Zcumul[ZSIZE] == 0)
	    rbin--;
	}
      while(Zcumul[ZSIZE] == 0 && rbin>0);


      for(zbin = 0; zbin <= ZSIZE; zbin++)
	Zcumul[zbin] = Zcumul[zbin] / Zcumul[ZSIZE];

      do
        q = drand48();
      while(q>0.995);
      
      for(zbin = 0; zbin < ZSIZE - 1; zbin++)
	{
	  if(q >= Zcumul[zbin] && q < Zcumul[zbin + 1])
	    break;
	}

      u = (q - Zcumul[zbin]) / (Zcumul[zbin + 1] - Zcumul[zbin]);

      z = (1 - u) * list_z[zbin] + u * list_z[zbin + 1];

      if(drand48() > 0.5)
	z = -z;

      q = drand48();
      phi = drand48() * PI * 2;

      xp_gas[i] = R * cos(phi);
      yp_gas[i] = R * sin(phi);
      zp_gas[i] = z;
    }



  for(i = 1 + n_disk; i <= N_GAS;)
    {
      q = drand48();
      
      zp_gas[i] = Z0 / 2 * log(q / (1 - q));

      q = drand48();

      R = H * sqrt(q) * HI_GasDiskScaleLength;

      phi = drand48() * PI * 2;

      xp_gas[i] = R * cos(phi);
      yp_gas[i] = R * sin(phi);


      if(fabs(zp_gas[i]) < 10*Z0)
	i++;
    }



  for(i = 1; i <= N_GAS; i++)
    mp_gas[i] = M_GAS / N_GAS;
}
