/**********************************************************************/

/*******************Taylor-bubble 2D sensitive analysis******************************

First determine the standard N_f = 50, Eo = 50, rho_r = 1000, mu_r = 100
Including:
# mesh analysis
# effect of initial bubble length
# effect of rho_r
# effect of mu_r
# effect of surface tension (Eo)
# effect of inverse viscosity number (N_f)

Longkai GUO
*/

#define REDUCED 1
#define ADAPT 1

#include "axi.h" 
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
//#include "navier-stokes/conserving.h"
#include "tension.h"

#if REDUCED
  #include "reduced.h"
#endif

#include "view.h"
#include "navier-stokes/perfs.h"

/******************** N_f &Eo &RHOR &MUR &l_taylor  ***********************/
double N_f = 50;
double Eo = 50;
double RHOR = 10;
double MUR = 100;
double l_taylor = 6; //bubble length ratio

#define length_ratio l_taylor

//#define N_f 50 (discuss the effect of N_f)
#define ratio 1

/**********************************************************/
#define Fr_0 1
#define _g 1 //in this case, consider U_0 = sqrt(gD_0)--->Fr_0^2 = U_0^2/(gD_0)

/* R_0 characteristic length (radias of lower tube) */
#define lambda_1  (2*(sqrt(1+2.44*pow(N_f,0.667))-1)/(2.44*pow(N_f,0.667))) 
// dimensionless liquid film thickness------refernce from chengsi's thesis tto guarantee the consistency with experiment
#define R_0 0.5
#define lambda  lambda_1*R_0 // liquid film
#define r_0 (R_0 - lambda)// or (R_0 - lambda)
#define R_ex R_0*ratio
#define L_lower 8
#define L_upper 8
#define L0 L_lower+L_upper

/* grid infromation */
#define maxlevel 11
#define minlevel 9

// boundary conditions
bid top_lower;
bid top_upper;

u.t[top_lower] = dirichlet(0);
u.n[top_lower] = dirichlet(0);
u.t[top_upper] = dirichlet(0);
u.n[top_upper] = dirichlet(0);

uf.n[bottom] = 0.;
// dirichlet boundary condition for left and right walls
u.n[right] = dirichlet(0);
u.t[right]  = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[left]  = dirichlet(0);

#define TIME 40


// geometry of Taylor-bubble
double sphere(double x, double y, coord center, double radius) {

  return ( sq(x ) + sq (y ) - sq (radius));

}

double geometry(double x, double y) {

  coord center;
  foreach_dimension()
  center.x = 0;
  double distance = 2;
  double s = sphere (x-(length_ratio+distance)*r_0, y, center, r_0);
  double left = x - distance*r_0;
  double right = -x + (length_ratio+distance)*r_0;
  double top = -y + r_0;
  double bottom = 0;
  double c = -min(min(min(top,left),right),bottom);
  double geom = min(s, c);

  return geom;
}

int main() {

// domain size
  size (L0);
  init_grid(1<<minlevel);
  
// dimensionless properties
  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./N_f;
  mu2 = mu1/MUR;
  f.sigma = 0;


  TOLERANCE = 1e-4;

#if REDUCED
  G.x -= _g;
  Z.x = _g;
#endif

run();  
}

void createFolder(const char* folder)

{
    char *str = NULL;
    str = (char *) malloc(sizeof(char) *1024);
    sprintf(str, "mkdir %s", folder); 
    system(str);   
}

event init (t = 0) {
// creat folders
createFolder("pressure");
createFolder("interface");
createFolder("tip");
createFolder("gfsview");
createFolder("velocity");
createFolder("field");
createFolder("shear");
createFolder("curvature");
createFolder("log_1");
createFolder("movie");

// computational domain
mask((x<=L_lower&&y>R_0)?top_lower:(x>L_lower&&y>R_ex)?top_upper:none);

// mesh refine initallly for interested region
refine ((x<=8&&(y-R_0<0)) && level < maxlevel);
refine ((x>=8&&(y-R_0*ratio<0)) && level < maxlevel);

fraction (f, geometry(x,y));
unrefine(y>R_0);
}


/**
We add the acceleration of gravity. */
#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= _g;
}
#endif


event logfile (t = 0; t <= TIME; t += 0.01) {

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 1024);
  sprintf (name, "log_1/out");
  static FILE * fp = fopen (name, "w");
  double xb = 0., yb = 0, vx = 0., vy = 0,  sb = 0., ke = 0;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:vx)  reduction(+:vy) reduction(+:sb) reduction(+:ke)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv;
    ke += dv()*rho(f[])*(sq(u.x[]) + sq(u.y[]));
  }
  
  fprintf (fp,"%g %g %g %g %g %g %g \n", t, sb, xb/sb, yb/sb, vx/sb, vy/sb, ke);
  fflush (fp);

}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){1e-3}, maxlevel,minlevel);
}
#endif

/**********************************************************/
#include "output_taylor.h"

// movie
event movies (t = 0; t <= TIME; t += 0.05) {

  char legend_time[80];
  double Mo=pow(Eo,3)/pow(N_f,4);
  sprintf(legend_time, "t* = %g\nN_f = %g\nEo=%g\nMo = %g", t,N_f,Eo,Mo);
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  clear();
// interface
  view (fov = 20, tx=-0.5,ty=0);
  draw_vof ("f",fc = {1,0,0}, lw = 2);
  squares("f",linear = true);
  draw_string(legend_time, 1, size = 30., lw = 2.);
  mirror({0,1}){
  draw_vof ("f",fc = {1,0,0}, lw = 2);
  squares("f",linear = true);

    }

char name[1000];
sprintf (name, "movie/Taylor-bubble.mp4");
save(name);

}