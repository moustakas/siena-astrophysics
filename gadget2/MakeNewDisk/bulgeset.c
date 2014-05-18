#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



static double R, z;


double set_bulge_velocities(void)
{
  int i;
  long dum;
  int iz, ir;
  double ur, uz;
  double vdisp_rz, vdisp_phi, vstream_phi;
  double vr, vphi;
  double vx, vy, vz;

  if(N_BULGE == 0)
    return 0;

  dum = drand48() * 1e8;

  printf("set bulge velocities...");
  fflush(stdout);

  for(i = 1; i <= N_BULGE; i++)
    {
      R = sqrt(xp_bulge[i] * xp_bulge[i] + yp_bulge[i] * yp_bulge[i]);
      z = fabs(zp_bulge[i]);

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




      vdisp_rz = VelDispRz_bulge[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispRz_bulge[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispRz_bulge[ir][iz + 1] * (1 - ur) * (uz) + VelDispRz_bulge[ir + 1][iz + 1] * (ur) * (uz);

      vdisp_phi = VelDispPhi_bulge[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispPhi_bulge[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispPhi_bulge[ir][iz + 1] * (1 - ur) * (uz) + VelDispPhi_bulge[ir + 1][iz + 1] * (ur) * (uz);

      vstream_phi = VelStreamPhi_bulge[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamPhi_bulge[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamPhi_bulge[ir][iz + 1] * (1 - ur) * (uz) + VelStreamPhi_bulge[ir + 1][iz + 1] * (ur) * (uz);


      if(vdisp_rz < 0)
	{
	  printf("in bulge: vdisp_rz:%g   %g %g %d %d \n", vdisp_rz, ur, uz, ir, iz);
	  vdisp_rz = -vdisp_rz;
	}
      if(vdisp_phi < 0)
	{
	  printf("in bulge: vdisp_phi:%g  %g %g %d %d\n", vdisp_phi, ur, uz, ir, iz);

	  vdisp_phi = -vdisp_phi;
	}

      vr = gasdev(&dum) * sqrt(vdisp_rz);
      vz = gasdev(&dum) * sqrt(vdisp_rz);

      vphi = vstream_phi + gasdev(&dum) * sqrt(vdisp_phi);


      vx = vr * xp_bulge[i] / R - vphi * yp_bulge[i] / R;
      vy = vr * yp_bulge[i] / R + vphi * xp_bulge[i] / R;

      vxp_bulge[i] = vx;
      vyp_bulge[i] = vy;
      vzp_bulge[i] = vz;
    }

  printf("done.\n");
  fflush(stdout);

  return 0;
}





double set_bulge_positions(void)
{
  int i, countr, countz;
  double q, R, phi, theta;

  if(N_BULGE == 0)
    return 0;

  srand48(222);


  for(i = 1, countr = countz = 0; i <= N_BULGE;)
    {
      q = drand48();
      R = 1 / (1 - sqrt(q)) - 1;

      phi = drand48() * PI * 2;
      theta = acos(drand48() * 2 - 1);

      xp_bulge[i] = R * sin(theta) * cos(phi);
      yp_bulge[i] = R * sin(theta) * sin(phi);
      zp_bulge[i] = R * cos(theta);

      /* now stretch */

      xp_bulge[i] *= A;
      yp_bulge[i] *= A;
      zp_bulge[i] *= A;


      mp_bulge[i] = M_BULGE / N_BULGE;

      if((R * A) > LL)
	countr++;
      else
	i++;
    }

  /*
     printf("bulge discarded:  %d  \n",countr);
   */

  return 0;
}
