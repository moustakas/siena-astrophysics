#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



static double R, z;



double set_halo_velocities(void)
{
  int i, rangecount=0;
  long dum;
  int iz, ir;
  double ur, uz;
  double vdisp_rz, vdisp_phi, vstream_phi;
  double vr, vphi;
  double vx, vy, vz;


  if(N_HALO == 0)
    return 0;

  dum = drand48() * 1e8;

  printf("set halo velocities...\n");
  fflush(stdout);

  for(i = 1; i <= N_HALO; i++)
    {
      R = sqrt(xp_halo[i] * xp_halo[i] + yp_halo[i] * yp_halo[i]);
      z = fabs(zp_halo[i]);

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
	  rangecount++;
	  /* 
	     printf("ir=%d out of range R=%g\n", ir, R);
	  */
	  vxp_halo[i] = 0;
	  vyp_halo[i] = 0;
	  vzp_halo[i] = 0;
	  continue;
	}

      if(iz < 0 || iz >= ZSIZE || uz < 0 || uz > 1)
	{
	  rangecount++;
	  /*
	  printf("iz=%d out of range Z=%g\n", iz, z);
	  */
	  vxp_halo[i] = 0;
	  vyp_halo[i] = 0;
	  vzp_halo[i] = 0;
	  continue;
	}


      vdisp_rz = VelDispRz_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispRz_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispRz_halo[ir][iz + 1] * (1 - ur) * (uz) + VelDispRz_halo[ir + 1][iz + 1] * (ur) * (uz);

      vdisp_phi = VelDispPhi_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispPhi_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispPhi_halo[ir][iz + 1] * (1 - ur) * (uz) + VelDispPhi_halo[ir + 1][iz + 1] * (ur) * (uz);
if((i==104317) || (i==104320))
{
	printf("ir= %d  iz= %d  ur= %g  uz= %g\n",ir,iz,ur,uz);
	printf("VelDispPhi_halo[ir][iz]= %g\n",VelDispPhi_halo[ir][iz]);
	printf("VelDispPhi_halo[ir+1][iz]= %g\n",VelDispPhi_halo[ir+1][iz]);
	printf("VelDispPhi_halo[ir][iz+1]= %g\n",VelDispPhi_halo[ir][iz+1]);
	printf("VelDispPhi_halo[ir+1][iz+1]= %g\n",VelDispPhi_halo[ir+1][iz+1]);
	printf("---\n");
}

      vstream_phi = VelStreamPhi_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamPhi_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamPhi_halo[ir][iz + 1] * (1 - ur) * (uz) + VelStreamPhi_halo[ir + 1][iz + 1] * (ur) * (uz);
if((i==104317) || (i==104320))
{
        printf("VelStreamPhi_halo[ir][iz]= %g\n",VelStreamPhi_halo[ir][iz]);
        printf("VelStreamPhi_halo[ir+1][iz]= %g\n",VelStreamPhi_halo[ir+1][iz]);
        printf("VelStreamPhi_halo[ir][iz+1]= %g\n",VelStreamPhi_halo[ir][iz+1]);
        printf("VelStreamPhi_halo[ir+1][iz+1]= %g\n",VelStreamPhi_halo[ir+1][iz+1]);
        fflush(stdout);
}

      if(vdisp_rz < 0)
	{
	  printf("in halo: vdisp_rz:%g   %g %g %d %d \n", vdisp_rz, ur, uz, ir, iz);
	  vdisp_rz = -vdisp_rz;
	}
      if(vdisp_phi < 0)
	{
	  printf("in halo: vdisp_phi:%g  %g %g %d %d\n", vdisp_phi, ur, uz, ir, iz);

	  vdisp_phi = -vdisp_phi;
	}

      vr = gasdev(&dum) * sqrt(vdisp_rz);
      vz = gasdev(&dum) * sqrt(vdisp_rz);

      vphi = vstream_phi + gasdev(&dum) * sqrt(vdisp_phi);
if(i==104317)
{
	printf("vphi= %g    vstream_phi= %g  vdisp_phi= %g  gasdev(&dum)= %g\n", vphi, vstream_phi, vdisp_phi, gasdev(&dum));
	fflush(stdout);
}


      vx = vr * xp_halo[i] / R - vphi * yp_halo[i] / R;
      vy = vr * yp_halo[i] / R + vphi * xp_halo[i] / R;

      vxp_halo[i] = vx;
      vyp_halo[i] = vy;
      vzp_halo[i] = vz;
if((i > 104312) && (i < 104322))
{
	printf("i= %d   R= %g vr= %g  vz= %g  vphi= %g  xyz= %g|%g|%g  v_xyz= %g|%g|%g\n",i,R,vr,vz,vphi,xp_halo[i],yp_halo[i],zp_halo[i],vx,vy,vz); fflush(stdout);
}
    }

  printf("rangecount=%d\n", rangecount);
  printf("done.\n");
  fflush(stdout);

  return 0;
}






double set_halo_positions(void)
{
  int i, countr, countz;
  double q, R, phi, theta;
  double halo_q_to_r(double q);

  if(N_HALO == 0)
    return 0;

  srand48(22);

  printf("set halo positions...\n");

  for(i = 1, countr = countz = 0; i <= N_HALO;)
    {

      do
	{
	  q = drand48();

	  R = halo_q_to_r(q);
	}
      while(R > 50 * R200);

      phi = drand48() * PI * 2;
      theta = acos(drand48() * 2 - 1);

      xp_halo[i] = R * sin(theta) * cos(phi);
      yp_halo[i] = R * sin(theta) * sin(phi);
      zp_halo[i] = R * cos(theta);

      i++;
    }

  for(i = 1; i <= N_HALO; i++)
    mp_halo[i] = M_HALO / N_HALO;

  return 0;
}


