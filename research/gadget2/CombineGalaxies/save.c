#include <stdio.h>
#include <math.h>
#include <stdlib.h>



#include "globvars.h"




void save_combined(char *fname, struct galaxy_data *g1, struct galaxy_data *g2)
{
  FILE *fd;
  int i, type, pc1, pc2, blklen, ntot_withmasses;
  struct io_header new_header;

#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);


  if(!(fd = fopen(fname, "w")))
    {
      fprintf(stderr, "Error in writing to '%s'\n", fname);
      exit(0);
    }

  new_header = header[0];

  for(i = 0; i < 6; i++)
    {
      new_header.npart[i] += header[1].npart[i];
      new_header.npartTotal[i] += header[1].npartTotal[i];

      if(header[0].mass[i] != header[1].mass[i])
	new_header.mass[i] = 0;
    }

  for(i = 1; i <= g2->Ntot; i++)
    g2->id[i] += g1->Ntot;




  blklen = sizeof(struct io_header);
  BLKLEN;
  fwrite(&new_header, sizeof(struct io_header), 1, fd);
  BLKLEN;


  blklen = 3 * (g1->Ntot + g2->Ntot) * sizeof(float);
  BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	fwrite(&g1->pos[pc1++][1], sizeof(float), 3, fd);

      for(i = 1; i <= header[1].npart[type]; i++)
	fwrite(&g2->pos[pc2++][1], sizeof(float), 3, fd);
    }
  BLKLEN;


  blklen = 3 * (g1->Ntot + g2->Ntot) * sizeof(float);
  BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	fwrite(&g1->vel[pc1++][1], sizeof(float), 3, fd);

      for(i = 1; i <= header[1].npart[type]; i++)
	fwrite(&g2->vel[pc2++][1], sizeof(float), 3, fd);
    }
  BLKLEN;


  blklen = (g1->Ntot + g2->Ntot) * sizeof(int);
  BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	fwrite(&g1->id[pc1++], sizeof(int), 1, fd);

      for(i = 1; i <= header[1].npart[type]; i++)
	fwrite(&g2->id[pc2++], sizeof(int), 1, fd);
    }
  BLKLEN;



  for(i = 0, ntot_withmasses = 0; i < 6; i++)
    {
      if(new_header.mass[i] == 0)
	ntot_withmasses += new_header.npart[i];
    }

  blklen = ntot_withmasses * sizeof(float);
  if(ntot_withmasses)
    BLKLEN;
  for(type = 0, pc1 = pc2 = 1; type < 6; type++)
    {
      for(i = 1; i <= header[0].npart[type]; i++)
	if(new_header.mass[type] == 0)
	  fwrite(&g1->m[pc1++], sizeof(float), 1, fd);
	else
	  pc1++;

      for(i = 1; i <= header[1].npart[type]; i++)
	if(new_header.mass[type] == 0)
	  fwrite(&g2->m[pc2++], sizeof(float), 1, fd);
	else
	  pc2++;
    }
  if(ntot_withmasses)
    BLKLEN;


  if(new_header.npart[0] > 0)
    {
      blklen = new_header.npartTotal[0] * sizeof(float);
      BLKLEN;
      for(i = 1; i <= header[0].npart[0]; i++)
	fwrite(&g1->u[i], sizeof(float), 1, fd);
      for(i = 1; i <= header[1].npart[0]; i++)
	fwrite(&g2->u[i], sizeof(float), 1, fd);
      BLKLEN;
    }

  fclose(fd);
}
