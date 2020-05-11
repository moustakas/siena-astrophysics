#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"

static double R, z;

void dump_veldisp_field(void)
{
  FILE *fd;
  int si;

  fd = fopen("deldisp.dat", "w");

  si = RSIZE + 1;
  fwrite(&si, sizeof(int), 1, fd);
  si = ZSIZE + 1;
  fwrite(&si, sizeof(int), 1, fd);

  fwrite(list_R, sizeof(double), RSIZE + 1, fd);
  fwrite(list_z, sizeof(double), ZSIZE + 1, fd);

  fwrite(&VelDispRz_disk[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
  fwrite(&VelDispPhi_disk[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
  fwrite(&VelStreamPhi_disk[0][0], sizeof(double), (ZSIZE + 1) * (RSIZE + 1), fd);
  fwrite(&epi_gamma2[0], sizeof(double), RSIZE + 1, fd);
  fwrite(&epi_kappa2[0], sizeof(double), RSIZE + 1, fd);

  fclose(fd);
}

void compute_velocity_dispersions_disk(void)
{
  int i, j;
  double z, R;
  double rho;

  if(N_DISK == 0)
    return;

  printf("disk velocity dispersion field...\n");
  fflush(stdout);

  for(i = 0; i <= RSIZE; i++)
    {
      for(j = 0; j <= ZSIZE; j++)
	{

	  xl[j + 1] = list_z[j];
	  yl[j + 1] = Dphi_z[i][j] * comp_rho_disk(list_R[i], list_z[j]);
	}

      spline(xl, yl, ZSIZE + 1, 1e40, 1e40, D2yl);

      for(j = ZSIZE - 1, VelDispRz_disk[i][ZSIZE] = 0; j >= 0; j--)
	{
	  VelDispRz_disk[i][j] = VelDispRz_disk[i][j + 1];
	  if(fabs(yl[j + 2]) > 1e-100 && fabs(yl[j + 1]) > 1e-100)
	    VelDispRz_disk[i][j] += qromb(splint_xl_yl_D2yl, list_z[j], list_z[j + 1]);
	}
    }

  for(i = 0; i <= RSIZE; i++)
    {
      for(j = 0; j <= ZSIZE; j++)
	{
	  xl[j + 1] = list_z[j];
	  yl[j + 1] = Dphi_z_dR[i][j] * comp_rho_disk(list_RplusdR[i], list_z[j]);

	}

      spline(xl, yl, ZSIZE + 1, 1e40, 1e40, D2yl);

      for(j = ZSIZE - 1, VelDispRz_dR_disk[i][ZSIZE] = 0; j >= 0; j--)
	{
	  VelDispRz_dR_disk[i][j] = VelDispRz_dR_disk[i][j + 1];
	  if(fabs(yl[j + 2]) > 1e-100 && fabs(yl[j + 1]) > 1e-100)
	    VelDispRz_dR_disk[i][j] += qromb(splint_xl_yl_D2yl, list_z[j], list_z[j + 1]);
	}

    }


  for(i = 0; i <= RSIZE; i++)
    {
      for(j = 0; j <= ZSIZE; j++)
	{

	  R = list_R[i];
	  z = list_z[j];

	  rho = comp_rho_disk(R, z);

	  if(rho > 0)
	    {
	      if(i > 0)
		VelDispPhi_disk[i][j] =
		  R / rho * (VelDispRz_dR_disk[i][j] - VelDispRz_disk[i][j]) / (list_RplusdR[i] - list_R[i]);
	      else
		VelDispPhi_disk[i][j] = 0;

	      VelDispRz_disk[i][j] /= rho;
	    }
	  else
	    VelDispRz_disk[i][j] = VelDispPhi_disk[i][j] = 0;

	  VelVc2_disk[i][j] = R * Dphi_R[i][j];

	  VelDispPhi_disk[i][j] *= RadialDispersionFactor;

	  VelDispPhi_disk[i][j] += VelVc2_disk[i][j] + VelDispRz_disk[i][j] * RadialDispersionFactor;


	  if(VelDispPhi_disk[i][j] > VelDispRz_disk[i][j] * RadialDispersionFactor / epi_gamma2[i])
	    {
	      VelStreamPhi_disk[i][j] =
		sqrt(VelDispPhi_disk[i][j] - VelDispRz_disk[i][j] * RadialDispersionFactor / epi_gamma2[i]);
	      VelDispPhi_disk[i][j] = VelDispRz_disk[i][j] * RadialDispersionFactor / epi_gamma2[i];
	    }
	  else
	    {
	      VelStreamPhi_disk[i][j] = 0;
	      VelDispPhi_disk[i][j] = VelDispRz_disk[i][j] * RadialDispersionFactor / epi_gamma2[i];
	    }


	  if(VelDispRz_disk[i][j] < 0)
	    VelDispRz_disk[i][j] = 0;

	  if(VelDispPhi_disk[i][j] < 0)
	    VelDispPhi_disk[i][j] = 0;
	}
    }

  printf("done.\n");
  fflush(stdout);

}







double intz_di(double k);
double intz_di_abs(double);

double intR_di(double k);
double intR_di_abs(double k);


double comp_Dphi_Z_disk_tree(double RR, double zz)
{
  double pos[3], acc[3];
  int nint;

  if(M_DISK > 0)
    {
      pos[0] = RR;
      pos[1] = 0;
      pos[2] = zz;
      
      nint = force_treeevaluate(pos, 0.01 * H, acc);

      return -G * acc[2];
    }
  else
    return 0;
}


double comp_Dphi_R_disk_tree(double RR, double zz)
{
  double pos[3], acc[3];
  int nint;

  if(M_DISK > 0)
    {
      pos[0] = RR;
      pos[1] = 0;
      pos[2] = zz;
      
      nint = force_treeevaluate(pos, 0.01 * H, acc);
  
      return -G * acc[0];
    }
  else
    return 0;
}


double comp_Dphi_z_disk(double RR, double zz)
{
  double pos[3], acc[3];

  if(M_DISK > 0)
    {
      pos[0] = RR;
      pos[1] = 0;
      pos[2] = zz;
      
      force_treeevaluate(pos, 0.01 * H, acc);
      
      return -G * acc[2];
    }
  else
    return 0;
}

/*

double comp_Dphi_z_disk(double RR, double zz)
{
  double comp_Dphi_z_disk_sph(double RR, double zz);
  double comp_Dphi_z_disk_exact(double RR, double zz);

  if(sqrt(RR * RR + zz * zz) > 10 * H)
    return comp_Dphi_z_disk_sph(RR, zz);
  else
    return comp_Dphi_z_disk_exact(RR, zz);

}

*/





double comp_Dphi_R_disk(double RR, double zz)
{
  double comp_Dphi_R_disk_sph(double RR, double zz);
  double comp_Dphi_R_disk_exact(double RR, double zz);

  if(RR > 0)
    {
      if(sqrt(RR * RR + zz * zz) > 10 * H)
	return comp_Dphi_R_disk_sph(RR, zz);
      else
	return comp_Dphi_R_disk_exact(RR, zz);
    }
  else
    return 0;
}






double comp_Dphi_z_disk_sph(double RR, double zz)
{
  double m;
  double r;


  r = sqrt(RR * RR + zz * zz);

  m = M_DISK * (1 - (1 + r / H) * exp(-r / H));

  return G * zz / (r * r * r) * m;
}



double comp_Dphi_z_disk_exact(double RR, double zz)
{
  int i;
  double dphiz;
  double Sigma0;
  double in1, in2, in3, bb;
  double deltaz, zpos;


  if(N_DISK == 0 && N_GAS == 0)
    return 0;

  if(fabs(zz) < 4 * Z0)
    {
      deltaz = (6.0 * Z0) / NSHEETS;

      dphiz = 0;

      for(i = 0; i < NSHEETS; i++)
	{
	  zpos = -3.0 * Z0 + (i + 0.5) * Z0 * 6.0 / NSHEETS;

	  R = RR;
	  z = zz - zpos;

	  Sigma0 =
	    (M_DISK) / (2 * PI * H * H) * deltaz / (2 * Z0) * pow(2 / (exp(zpos / Z0) + exp(-zpos / Z0)), 2);

	  in1 = qromb(intz_di, 0, 2 / H);

	  bb = 2;
	  do
	    {
	      in2 = qromb(intz_di, bb / H, (bb + 2) / H);
	      in3 = qromb(intz_di_abs, bb / H, (bb + 2) / H);
	      in1 += in2;
	      bb += 2;
	    }
	  while(fabs(in3 / in1) > 1e-2);

	  dphiz += 2 * PI * G * Sigma0 * H * H * (in1);
	}
      return dphiz;
    }
  else
    {
      R = RR;
      z = zz;
      Sigma0 = (M_DISK) / (2 * PI * H * H);

      in1 = qromb(intz_di, 0, 2 / H);

      bb = 2;
      do
	{
	  in2 = qromb(intz_di, bb / H, (bb + 2) / H);
	  in3 = qromb(intz_di_abs, bb / H, (bb + 2) / H);
	  in1 += in2;
	  bb += 2;
	}
      while(fabs(in3 / in1) > 1e-2);

      dphiz = 2 * PI * G * Sigma0 * H * H * (in1);

      return dphiz;

    }
}








double comp_Dphi_R_disk_sph(double RR, double zz)
{
  double m;
  double r;


  r = sqrt(RR * RR + zz * zz);

  m = M_DISK * (1 - (1 + r / H) * exp(-r / H));

  return G * RR / (r * r * r) * m;
}



double comp_Dphi_R_disk_exact(double RR, double zz)
{
  int i;
  double dphiR;
  double Sigma0;
  double in1, in2, in3, bb;
  double deltaz, zpos;


  if(N_DISK == 0 && N_GAS == 0)
    return 0;

  if(fabs(zz) < 4 * Z0)
    {
      deltaz = (6.0 * Z0) / NSHEETS;

      dphiR = 0;

      for(i = 0; i < NSHEETS; i++)
	{
	  zpos = -3.0 * Z0 + (i + 0.5) * Z0 * 6.0 / NSHEETS;

	  R = RR;
	  z = zz - zpos;

	  Sigma0 =
	    (M_DISK) / (2 * PI * H * H) * deltaz / (2 * Z0) * pow(2 / (exp(zpos / Z0) + exp(-zpos / Z0)), 2);

	  in1 = qromb(intR_di, 0, 2 / H);

	  bb = 2;
	  do
	    {
	      in2 = qromb(intR_di, bb / H, (bb + 2) / H);
	      in3 = qromb(intR_di_abs, bb / H, (bb + 2) / H);
	      in1 += in2;
	      bb += 2;
	    }
	  while(fabs(in3 / in1) > 1e-2);

	  dphiR += 2 * PI * G * Sigma0 * H * H * (in1);
	}
      return dphiR;
    }
  else
    {
      R = RR;
      z = zz;
      Sigma0 = (M_DISK) / (2 * PI * H * H);

      in1 = qromb(intR_di, 0, 2 / H);

      bb = 2;
      do
	{
	  in2 = qromb(intR_di, bb / H, (bb + 2) / H);
	  in3 = qromb(intR_di_abs, bb / H, (bb + 2) / H);
	  in1 += in2;
	  bb += 2;
	}
      while(fabs(in3 / in1) > 1e-2);

      dphiR = 2 * PI * G * Sigma0 * H * H * (in1);

      return dphiR;

    }
}






double comp_Dphi_R_disk_razorthin(double RR, double zz)
{
  double Sigma0, y;
  double dphidR;

  if(RR > 0)
    {

      Sigma0 = (M_DISK) / (2 * PI * H * H);
      y = RR / (2 * H);

      if(y > 1e-4)
	dphidR = 2 * PI * G * Sigma0 * y * (bessi0(y) * bessk0(y) - bessi1(y) * bessk1(y));
      else
	dphidR = 0;

      return dphidR;
    }
  else
    return 0;
}









double comp_rho_disk(double R, double z)
{
  double x;

  x =
    (1 - GasFraction) * (M_DISK) / (4 * PI * H * H * Z0) * exp(-R / H) * pow(2 / (exp(z / Z0) + exp(-z / Z0)),
									     2);

  return x;
}



double intz_di(double k)
{
  if(z > 0)
    return (bessj0(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5));
  else
    return (-bessj0(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5));
}


double intz_di_abs(double k)
{
  if(z > 0)
    return fabs(bessj0(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5));
  else
    return fabs(-bessj0(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5));
}






double intR_di(double k)
{
  if(z >= 0)
    return bessj1(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5);
  else
    return bessj1(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5);
}

double intR_di_abs(double k)
{
  if(z >= 0)
    return fabs(bessj1(k * R) * k * exp(-z * k) / pow(1 + k * k * H * H, 1.5));
  else
    return fabs(bessj1(k * R) * k * exp(z * k) / pow(1 + k * k * H * H, 1.5));
}





double mass_cumulative_disk(double R)
{
  return M_DISK * (1 - (1 + R / H) * exp(-R / H));
}
