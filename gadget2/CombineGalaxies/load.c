#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "globvars.h"
#include "nrsrc/nrutil.h"




void load_particles(char *fname, struct galaxy_data *g, struct io_header *header)
{
  FILE *fd;
  int i, type, pc, blksize;
  int ntot_withmasses;

#define SKIP  fread(&blksize,sizeof(int),1,fd);

  if(!(fd = fopen(fname, "r")))
    {
      fprintf(stderr, "File '%s' not found.\n", fname);
      exit(0);
    }

  SKIP;
  fread(header, sizeof(struct io_header), 1, fd);
  SKIP;



  g->Ngas = header->npart[0];
  g->Ndm = header->npart[1];

  for(i = 0, g->Ntot = 0; i < 6; i++)
    g->Ntot += header->npart[i];


  g->pos = matrix(1, g->Ntot, 1, 3);
  g->vel = matrix(1, g->Ntot, 1, 3);
  g->id = ivector(1, g->Ntot);
  g->m = vector(1, g->Ntot);

  if(g->Ngas)
    g->u = vector(1, g->Ngas);


  SKIP;
  fread(&g->pos[1][1], sizeof(float), 3 * g->Ntot, fd);
  SKIP;

  SKIP;
  fread(&g->vel[1][1], sizeof(float), 3 * g->Ntot, fd);
  SKIP;

  SKIP;
  fread(&g->id[1], sizeof(int), g->Ntot, fd);
  SKIP;


  for(i = 0, ntot_withmasses = 0; i < 6; i++)
    {
      if(header->mass[i] == 0)
	ntot_withmasses += header->npart[i];
    }

  if(ntot_withmasses)
    SKIP;
  for(type = 0, pc = 1; type < 6; type++)
    {
      for(i = 0; i < header->npart[type]; i++)
	{
	  if(header->mass[type] == 0)
	    fread(&g->m[pc], sizeof(float), 1, fd);
	  else
	    g->m[pc] = header->mass[type];

	  pc++;
	}
    }
  if(ntot_withmasses)
    SKIP;

  if(g->Ngas)
    {
      SKIP;
      fread(&g->u[1], sizeof(float), g->Ngas, fd);
      SKIP;
    }

  fclose(fd);


  for(i = 1, g->Mtot = 0; i <= g->Ntot; i++)
    g->Mtot += g->m[i];

  for(i = 1 + g->Ngas, g->Mdm = 0; i <= (g->Ngas + g->Ndm); i++)
    g->Mdm += g->m[i];
}
