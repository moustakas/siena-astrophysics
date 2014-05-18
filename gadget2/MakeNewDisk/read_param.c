#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "globvars.h"
#include "prototypes.h"


void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "CC");
  addr[nt] = &CC;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "V200");
  addr[nt] = &V200;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "LAMBDA");
  addr[nt] = &LAMBDA;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "MD");
  addr[nt] = &MD;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "MBH");
  addr[nt] = &MBH;
  id[nt++] = FLOAT;


  strcpy(tag[nt], "MB");
  addr[nt] = &MB;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "JD");
  addr[nt] = &JD;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "GasFraction");
  addr[nt] = &GasFraction;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "DiskHeight");
  addr[nt] = &DiskHeight;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BulgeSize");
  addr[nt] = &BulgeSize;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "N_HALO");
  addr[nt] = &N_HALO;
  id[nt++] = INT;

  strcpy(tag[nt], "N_DISK");
  addr[nt] = &N_DISK;
  id[nt++] = INT;

  strcpy(tag[nt], "N_GAS");
  addr[nt] = &N_GAS;
  id[nt++] = INT;

  strcpy(tag[nt], "N_BULGE");
  addr[nt] = &N_BULGE;
  id[nt++] = INT;


  strcpy(tag[nt], "RadialDispersionFactor");
  addr[nt] = &RadialDispersionFactor;
  id[nt++] = FLOAT;


  strcpy(tag[nt], "HI_GasMassFraction");
  addr[nt] = &HI_GasMassFraction;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "HI_GasDiskScaleLength");
  addr[nt] = &HI_GasDiskScaleLength;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "MaxGasDiskHeight");
  addr[nt] = &MaxGasDiskHeight;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputFile");
  addr[nt] = OutputFile;
  id[nt++] = STRING;


  strcpy(tag[nt], "FactorEVP");
  addr[nt] = &FactorEVP;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "FactorSN");
  addr[nt] = &FactorSN;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "MaxSfrTimescale");
  addr[nt] = &MaxSfrTimescale;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "TempSupernova");
  addr[nt] = &TempSupernova;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "TempClouds");
  addr[nt] = &TempClouds;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "FactorForSofterEQS");
  addr[nt] = &FactorForSofterEQS;
  id[nt++] = FLOAT;

#ifdef REDSHIFT_SCALING
  strcpy(tag[nt], "REDSHIFT");
  addr[nt] = &REDSHIFT;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omega_m0");
  addr[nt] = &Omega_m0;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Omega_L0");
  addr[nt] = &Omega_L0;
  id[nt++] = FLOAT;
#endif

  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  buf[0] = 0;
	  fgets(buf, 200, fd);

	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy(addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
	      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
	      errorFlag = 1;
	    }
	}
      fclose(fd);

    }
  else
    {
      fprintf(stdout, "Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	  fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }

  if(errorFlag)
    {
      exit(1);
    }



#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
