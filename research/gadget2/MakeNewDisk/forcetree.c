#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "globvars.h"
#include "prototypes.h"
#include "forcetree.h"


static struct NODE
{
  float len;			/* sidelength of treenode */
  float center[3];		/* geometrical center */
  union
  {
    int suns[8];		/* ancestor nodes */
    struct
    {
      float s[3];
      float mass;
      int sibling;
      int nextnode;
    }
    d;
  }
  u;
}
 *Nodes, *Nodes_base;



static int MaxNodes;		/* maximum allowed number of internal nodes */
static int MaxPart;
static int *Nextnode;		/* gives next node in tree walk  (nodes array) */


static int last;		/* auxialiary variable used to set-up non-recursive walk */


int force_treebuild(void)
{
  int i, j, subnode = 0, parent = -1, numnodes;
  int nfree, th, nn;
  double lenhalf;
  struct NODE *nfreep;
  double len, xmin[3], xmax[3];


  Nodes = Nodes_base - MaxPart;

  /* select first node */
  nfree = MaxPart;		/* index */
  nfreep = &Nodes[nfree];

  /* create an empty  root node  */

  for(j = 0; j < 3; j++)	/* find enclosing rectangle */
    xmin[j] = xmax[j] = P[0].Pos[j];

  for(i = 1; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	if(P[i].Pos[j] > xmax[j])
	  xmax[j] = P[i].Pos[j];
	if(P[i].Pos[j] < xmin[j])
	  xmin[j] = P[i].Pos[j];
      }

  for(j = 1, len = xmax[0] - xmin[0]; j < 3; j++)	/* determine maxmimum extension */
    if((xmax[j] - xmin[j]) > len)
      len = xmax[j] - xmin[j];

  for(j = 0; j < 3; j++)
    nfreep->center[j] = (xmin[j] + xmax[j]) / 2;
  nfreep->len = len;

  for(i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;



  nfreep = &Nodes[nfree];

  for(i = 0; i < NumPart; i++)	/* insert all particles */
    {
      th = MaxPart;

      while(1)
	{
	  if(th >= MaxPart)	/* we are dealing with an internal node */
	    {
	      subnode = 0;
	      if(P[i].Pos[0] > Nodes[th].center[0])
		subnode += 1;
	      if(P[i].Pos[1] > Nodes[th].center[1])
		subnode += 2;
	      if(P[i].Pos[2] > Nodes[th].center[2])
		subnode += 4;

	      nn = Nodes[th].u.suns[subnode];

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = th;	/* note: subnode can still be used in the next step of the walk */
		  th = nn;
		}
	      else
		{
		  /* here we have found an empty slot where we can 
		   * attach the new particle as a leaf 
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* we try to insert into a leaf with a single particle
	       * need to generate a new internal node at this point 
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;

	      subnode = 0;
	      if(P[th].Pos[0] > nfreep->center[0])
		subnode += 1;
	      if(P[th].Pos[1] > nfreep->center[1])
		subnode += 2;
	      if(P[th].Pos[2] > nfreep->center[2])
		subnode += 4;

	      nfreep->u.suns[subnode] = th;

	      th = nfree;	/* resume trying to insert the new particle at 
				   the newly created internal node */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  printf("maximum number %d of tree-nodes reached.\n", MaxNodes);
		  printf("for particle %d\n", i);
		  exit(1);
		}
	    }
	}
    }



  /* now compute the multipole moments recursively */
  last = -1;
  force_update_node_recursive(MaxPart, -1);

  if(last >= MaxPart)
    Nodes[last].u.d.nextnode = -1;
  else
    Nextnode[last] = -1;

  return numnodes;
}



void force_update_node_recursive(int no, int sib)
{
  int j, jj, p, pp = 0, nextsib, suns[8];
  double s[3], mass;

  if(no >= MaxPart)		/* internal node */
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= MaxPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib);

	      if(p >= MaxPart)	/* an internal node or pseudo particle */
		{
		  mass += Nodes[p].u.d.mass;	/* we assume a fixed particle mass */
		  s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		  s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		  s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		}
	      else		/* a particle */
		{
		  mass += P[p].Mass;
		  s[0] += P[p].Mass * P[p].Pos[0];
		  s[1] += P[p].Mass * P[p].Pos[1];
		  s[2] += P[p].Mass * P[p].Pos[2];
		}
	    }
	}

      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;

      Nodes[no].u.d.sibling = sib;
    }
  else				/* single particle */
    {
      if(last >= 0)
	{
	  if(last >= MaxPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    Nextnode[last] = no;
	}

      last = no;
    }
}




int force_treeevaluate(double *pos, double softening, double *acc)
{
  struct NODE *nop = 0;
  int no, ninteractions;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z;

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;


  pos_x = pos[0];
  pos_y = pos[1];
  pos_z = pos[2];

  h = 2.8 * softening;
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;

  no = MaxPart;			/* first node */

  while(no >= 0)
    {
      if(no < MaxPart)
	{
	  /* the index of the node is the index of the particle */
	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;

	  r2 = dx * dx + dy * dy + dz * dz;

	  mass = P[no].Mass;

	  no = Nextnode[no];
	}
      else			/* we have an  internal node */
	{
	  nop = &Nodes[no];

	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;

	  r2 = dx * dx + dy * dy + dz * dz;

	  mass = nop->u.d.mass;

/*
	  printf("nop->len=%g  r=%g  %g %g\n",
		 nop->len, sqrt(r2), nop->len * nop->len , r2 * ErrTolTheta * ErrTolTheta);
*/

	  if(nop->len * nop->len > r2 * ErrTolTheta * ErrTolTheta)
	    {
	      /* open cell */
	      no = nop->u.d.nextnode;
	      continue;
	    }

	  /* check in addition whether we lie inside the cell */

	  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
	    {
	      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
		    {
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}
	    }

	  no = nop->u.d.sibling;	/* ok, node can be used */
	}

      r = sqrt(r2);

      if(r >= h)
	fac = mass / (r2 * r);
      else
	{
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

      ninteractions++;
    }

  /* store result at the proper place */

  acc[0] = acc_x;
  acc[1] = acc_y;
  acc[2] = acc_z;

  return ninteractions;
}

















int force_treeevaluate_direct(double *pos, double softening, double *acc)
{
  int no;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z;

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;

  pos_x = pos[0];
  pos_y = pos[1];
  pos_z = pos[2];

  h = 2.8 * softening;
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;

  for(no = 0; no < NumPart; no++)
    {
      /* the index of the node is the index of the particle */
      dx = P[no].Pos[0] - pos_x;
      dy = P[no].Pos[1] - pos_y;
      dz = P[no].Pos[2] - pos_z;

      r2 = dx * dx + dy * dy + dz * dz;

      mass = P[no].Mass;

      r = sqrt(r2);

      if(r >= h)
	fac = mass / (r2 * r);
      else
	{
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

    }

  /* store result at the proper place */

  acc[0] = acc_x;
  acc[1] = acc_y;
  acc[2] = acc_z;

  return 0;
}







/* this function allocates memory used for storage of the tree
 * and auxiliary arrays for tree-walk and link-lists.
 */
void force_treeallocate(int maxnodes, int maxpart)	/* usually maxnodes=0.7*maxpart is sufficient */
{
  size_t bytes, allbytes = 0;

  MaxNodes = maxnodes;
  MaxPart = maxpart;

  if(!(Nodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
      exit(3);
    }
  allbytes += bytes;

  if(!(Nextnode = malloc(bytes = (maxpart) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%d bytes)\n", maxpart, bytes);
      exit(0);
    }
  allbytes += bytes;

  printf("\nUse %g MByte for BH-tree.\n\n", allbytes / (1024.0 * 1024.0));
}


/* free the allocated memory
 */
void force_treefree(void)
{
  free(Nextnode);
  free(Nodes_base);
}
