/*** code to combine two galaxies such that they collide
     on a parabolic encounter ***/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "globvars.h"




int main(int argc, char *argv[])
{
  double theta1, phi1, theta2, phi2;
  double rmin, rstart, ecc=1.0;
  char gal_fname1[100], gal_fname2[100], gal_output[100];

#ifdef ECC_NONZERO
  if(argc != 11)
#else
  if(argc != 10)
#endif
    {
      printf("\nwrong number of arguments\n");
      printf
#ifdef ECC_NONZERO
	("call with:\n\n<fname_gal1> <theta1> <phi1>\n<fname_gal2> <theta2> <phi2>\n<rmin> <rstart>\n<eccentricity>\n<fname_galout>\n\n");
#else
	("call with:\n\n<fname_gal1> <theta1> <phi1>\n<fname_gal2> <theta2> <phi2>\n<rmin> <rstart>\n<fname_galout>\n\n");
#endif
      printf("(angles in degrees.)\n\n");
      exit(0);
    }

#ifdef ECC_NONZERO
#ifdef ORBITCORRECTION
  printf(" can't have ECC_NONZERO and ORBITCORRECTION turned on at the same time.\n");
  printf(" please fix.\n");
  exit(0);
#endif
#endif

  strcpy(gal_fname1, argv[1]);
  theta1 = atof(argv[2]);
  phi1 = atof(argv[3]);

  strcpy(gal_fname2, argv[4]);
  theta2 = atof(argv[5]);
  phi2 = atof(argv[6]);

  rmin = atof(argv[7]);
  rstart = atof(argv[8]);

#ifdef ECC_NONZERO
  ecc= atof(argv[9]);
  strcpy(gal_output, argv[10]);
#else
  strcpy(gal_output, argv[9]);
#endif

  load_particles(gal_fname1, &gal[0], &header[0]);
  load_particles(gal_fname2, &gal[1], &header[1]);

  turn_galaxy(&gal[0], theta1, phi1);
  turn_galaxy(&gal[1], theta2, phi2);


  move_galaxies(&gal[0], &gal[1], rmin, rstart, ecc);

  save_combined(gal_output, &gal[0], &gal[1]);

  return 0;
}
