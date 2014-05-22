
#define PI 3.1415926
#define G  43007.0                   /* gravitational constant */




extern struct galaxy_data {
  double Mdark, Scalefac;  
  int Ntot, Ngas, Ndm;
  double Mtot, Mdm;
  float **pos,**vel,*m,*u;
  int   *id;
} gal[2];

/* Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header[2];



void move_galaxies(struct galaxy_data *g1,
		   struct galaxy_data *g2,double rmin,double rstart,double ecc);


void load_particles(char *fname, struct galaxy_data *g, 
		    struct io_header *header);

void turn_galaxy(struct galaxy_data *g,double omega,double incl);

void save_combined(char *fname,struct galaxy_data *g1,struct galaxy_data *g2);



