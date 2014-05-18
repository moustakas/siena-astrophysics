#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"



void compute_vstream_gas(void)
{
  int i, j;
  double P1, P2, vphipress, vphi;


  printf("comp vstream gas\n");

  for(i = 0; i < RSIZE; i++)
    {
      for(j = 0; j <= ZSIZE; j++)
	{
	  if(i == 0 || i == RSIZE)
	    VelStreamGas[i][j] = 0;
	  else
	    {
	      VelStreamGas[i][j] = 0;

	      vphi = Dphi_R[i][j];

	      P2 = eqn_of_state(RhoGas[i][j]);
	      P1 = eqn_of_state(1.1*RhoGas[i][j]);

	      if(RhoGas[i][j] > 0 && RhoGas[i+1][j]> 0 )
		{
		  vphipress = - log(P1 / P2) / log(1.1) * P2 / RhoGas[i][j] / H;

                  if(vphi+vphipress<0) 
                    {
                      /*
                      printf("r=%g vphi=%g vphipress=%g j=%d  %g %g %g %g  %g %g %g %g\n", list_R[i], sqrt(vphi), sqrt(-vphipress), j, 
                             P1,  P2, RhoGas[i][j], RhoGas[i+1][j],
                             log(P1 / P2) / log(1.1),
                             log(RhoGas[i+1][j] / RhoGas[i][j]), log(list_R[i+1]/list_R[i]),
                             list_R[i]);
                      */
                      vphi = 0;
                    }
                  else
                    vphi += vphipress;
		}

	      if(vphi > 0)
		{
		  VelStreamGas[i][j] = sqrt(vphi * list_R[i]);
		      
		  /*
		    if(j==0)
		    printf("R=%g  vphi=%g  %g P2=%g P1=%g rho2=%g rho1=%g  delta=%g\n", list_R[i], 
		    VelStreamGas[i][j], sqrt(Dphi_R[i][j] * list_R[i]),
		    P2, P1, RhoGas[i][j], RhoGas[i-1][j],
		    log(P2/P1)/log(list_R[i]/list_R[i-1]) * P2/ RhoGas[i][j]/list_R[i]);
		  */
		}
	    }
	}
    }

  for(j = 0; j <= ZSIZE; j++)
    {
      VelStreamGas[1][j] = VelStreamGas[2][j];
    }

  printf("done\n");
}
