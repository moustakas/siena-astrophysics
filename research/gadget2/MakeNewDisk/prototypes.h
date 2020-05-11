
#include <stdio.h>

#include "forcetree.h"

void init_clouds(void);
double halo_mass(double r);

void integrate_gasdensity(void);
void integrate_and_adjust(void);

void read_parameterfile(char *fname);

void determine_cumulative_gasdensity(void);
void init_central_densities(void);
double surface_density_disk(double r);

void dump_gas_density(void);
void compute_vertical_force_field(void);

void integrate_surfacedensity(void);

void dump_eqn_of_state(void);
double comp_Dphi_Z_disk_tree(double RR, double zz);

double eqn_of_state(double rho);

void  set_dummy_particles(void);

double comp_Dphi_R_disk_tree(double RR, double zz);

void compute_vstream_gas(void);

void dump_veldisp_field(void);

double set_halo_positions(void);
double set_disk_positions(void);
double set_bulge_positions(void);

double set_halo_velocities(void);
double set_disk_velocities(void);
double set_bulge_velocities(void);

void compute_velocity_dispersions_halo(void);
void compute_velocity_dispersions_disk(void);
void compute_velocity_dispersions_bulge(void);
void compute_local_escape_speed(void);
void plot_toomre_stability(FILE *fd);
void allocate_memory(int NumPart);

void set_gas_velocities(void);
void plot_circular_speeds(FILE *fd);

void compute_radial_force_field(void);

void save_particles(char *fname);
void set_gas_positions(void);
void init_units(void);
void init(void);

void structure(void);
void write_header(FILE *fd);

void prepare_cumlative_profile(void);
double mass_cumulative_disk(double);
double mass_cumulative_bulge(double);
double disk_angmomentum(void);
double additional_mass_in_halo_cutoff(void);
double fc(double);
double gc(double);
void solve_mass_shells(void);
void setup_massprofile(void);


double comp_Dphi_z_halo(double R,double z);
double comp_Dphi_R_halo(double R,double z);
double comp_rho_halo(double R,double z);

double comp_Dphi_R_disk_razorthin(double RR,double zz);

double comp_Dphi_z_bulge(double R,double z);
double comp_Dphi_R_bulge(double R,double z);
double comp_rho_bulge(double R,double z);


double comp_Dphi_z_disk(double R,double z);
double comp_Dphi_R_disk(double R,double z);
double comp_rho_disk(double R,double z);


double comp_Dphi_z(double R,double z);
double comp_Dphi_R(double R,double z);

double epicyclic_kappa2(double R);


double splint_zs_y_D2y(double z);
double splint_xl_yl_D2yl(double t);


double   drand48(void);

void write_cumulative_mass(void);



double mass_cumulative_disk(double);
double mass_cumulative_bulge(double);
/* returns the maximum of two double
 */
double dmax(double x, double y);
double dmin(double x, double y);
int imax(int x, int y);
int imin(int x, int y);
