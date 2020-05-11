
#include <stdio.h>
#include <math.h>
#include <stdlib.h>



#include "globvars.h"




void turn_galaxy(struct galaxy_data *g, double theta, double phi)
{
  int i;
  void turn_1(float xyz[3], double phi);
  void turn_2(float xyz[3], double phi);


  for(i = 1; i <= g->Ntot; i++)
    {
      turn_2(&g->pos[i][1], theta);
      turn_1(&g->pos[i][1], phi);

      turn_2(&g->vel[i][1], theta);
      turn_1(&g->vel[i][1], phi);
    }
}


void turn_1(float xyz[3], double phi)	/* turns around z-axes by phi */
{
  double xx, yy, zz, co, si;

  co = cos(phi / 180.0 * PI);
  si = sin(phi / 180.0 * PI);

  xx = co * xyz[0] - si * xyz[1];

  yy = +si * xyz[0] + co * xyz[1];
  zz = xyz[2];

  xyz[0] = xx;
  xyz[1] = yy;
  xyz[2] = zz;
}



void turn_2(float xyz[3], double phi)	/* turns around y-axes by phi */
{
  double xx, yy, zz, co, si;

  co = cos(phi / 180.0 * PI);
  si = sin(phi / 180.0 * PI);

  xx = co * xyz[0] + si * xyz[2];
  yy = xyz[1];
  zz = -si * xyz[0] + co * xyz[2];


  xyz[0] = xx;
  xyz[1] = yy;
  xyz[2] = zz;
}
