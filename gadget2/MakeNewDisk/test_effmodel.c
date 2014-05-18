#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "prototypes.h"
#include "globvars.h"



int main(int argc, char *argv[])
{
  char parameter_filename[100];
  char buf[1000];
  FILE *fd;
  int rep;

  if(argc != 2)
    {
      fprintf(stderr, "\n\nwrong argument(s).  Specify a parameterfile.\n\n");
      exit(0);
    }
  strcpy(parameter_filename, argv[1]);

  read_parameterfile(parameter_filename);

  init_units();			/* set system of units */

  init_clouds();

  printf("done testing effmodel.\n");


  return 0;
}




