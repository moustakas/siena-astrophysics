#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "prototypes.h"
#include "globvars.h"


static double R, z;


double set_disk_velocities(void)
{
  int i;
  long dum;

  int iz, ir;
  double ur, uz;
  double vdisp_rz, vdisp_phi, vstream_phi;
  double vr, vphi;
  double vx, vy, vz;


  if(N_DISK == 0)
    return 0;

  dum = drand48() * 1e8;


  printf("set disk velocities...");
  fflush(stdout);

  for(i = 1; i <= N_DISK; i++)
    {
      R = sqrt(xp_disk[i] * xp_disk[i] + yp_disk[i] * yp_disk[i]);
      z = fabs(zp_disk[i]);

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

      if(ir < 0 || ir >= RSIZE || ur < 0 || ur > 1)
	{
	  printf("ir=%d out of range ur=%g\n", ir, ur);
	  exit(0);
	}

      if(iz < 0 || iz >= ZSIZE || uz < 0 || uz > 1)
	{
	  printf("iz=%d out of range uz=%g\n", iz, uz);
	  exit(0);
	}


      vdisp_rz = VelDispRz_disk[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispRz_disk[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispRz_disk[ir][iz + 1] * (1 - ur) * (uz) + VelDispRz_disk[ir + 1][iz + 1] * (ur) * (uz);

      vdisp_phi = VelDispPhi_disk[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispPhi_disk[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispPhi_disk[ir][iz + 1] * (1 - ur) * (uz) + VelDispPhi_disk[ir + 1][iz + 1] * (ur) * (uz);

      vstream_phi = VelStreamPhi_disk[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamPhi_disk[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamPhi_disk[ir][iz + 1] * (1 - ur) * (uz) + VelStreamPhi_disk[ir + 1][iz + 1] * (ur) * (uz);


      if(vdisp_rz < 0)
	{
	  printf("in disk: vdisp_rz:%g   %g %g %d %d \n", vdisp_rz, ur, uz, ir, iz);
	  vdisp_rz = -vdisp_rz;
	}
      if(vdisp_phi < 0)
	{
	  printf("in disk: vdisp_phi:%g  %g %g %d %d\n", vdisp_phi, ur, uz, ir, iz);

	  vdisp_phi = -vdisp_phi;
	}


      vr = gasdev(&dum) * sqrt(vdisp_rz);  /*  * 1.5; */
      vz = gasdev(&dum) * sqrt(vdisp_rz);

      vphi = vstream_phi + gasdev(&dum) * sqrt(vdisp_phi);


      vx = vr * xp_disk[i] / R - vphi * yp_disk[i] / R;
      vy = vr * yp_disk[i] / R + vphi * xp_disk[i] / R;

      vxp_disk[i] = vx;
      vyp_disk[i] = vy;
      vzp_disk[i] = vz;

    }

  printf("done.\n");
  fflush(stdout);

  return 0;
}




double set_disk_positions(void)
{
  int i, countr, countz;
  double q, f, f_, Rold, phi;


  if(N_DISK == 0)
    return 0;

  srand48(222);


  printf("set disk positions...\n");


  for(i = 1, countr = countz = 0; i <= N_DISK;)
    {
      q = drand48();

      zp_disk[i] = Z0 / 2 * log(q / (1 - q));

      q = drand48();

      R = 1.0;
      do
	{
	  f = (1 + R) * exp(-R) + q - 1;
	  f_ = -R * exp(-R);

	  Rold = R;
	  R = R - f / f_;
	}
      while(fabs(R - Rold) / R > 1e-7);

      R *= H;

      phi = drand48() * PI * 2;

      xp_disk[i] = R * cos(phi);
      yp_disk[i] = R * sin(phi);


      if(R > LL || fabs(zp_disk[i]) > LL)
	countr++;
      else
	i++;
    }

  for(i = 1; i <= N_DISK; i++)
    mp_disk[i] = (M_DISK - M_GAS) / N_DISK;

  return R;
}
