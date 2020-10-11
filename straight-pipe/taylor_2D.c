/**********************************************************************/

/*******************Taylor-bubble 2D test cases, compared with experimental data from chengsi******************************

Including:
# taylor bubble rising in straight vertical pipe
    1) 60%, 70%, 80%, 90%, 95%
    2) change ratio: 1.72, 1.45, 1.33, 1.24, 1.12, 0.9, 0.81, 0.69
# taylor bubble rising in sudden expansion/contraction vertical pipe
# taylor bubble rising in gradual expansion/contraction vertical pipe
# Popinet's trick (bview+dump)  ------ this makes possible in 3D cases, but interface not smooth
# heat transfer ----- diffusion.h

Three main dimensionless numbers are used in the code:

Eo = Buoyancy_force/surface_tension = (\Delta \rho g D^2)/\sigma
Nf = sqrt(\Delta\rho g D^3)/\mu_l // to describe the properties of continuous phase
Mo = g\mu_l^4\Delta\rho / (rho_l^2 \sigma^3) = Eo^3/Nf^4

Other dimensionless numbers:
\rho_r = \rho_l/\rho_g
\mu_r = \mu_l/\mu_g

Longkai GUO
*/
#define AXIS 1
#define REDUCED 1
#define ADAPT 1
#define HEAT 0  // heat transfer
#define VERTICAL 0 //vertical tube
#define GRADUAL 0 //gradual tube
#define no_change 1 //vertical tube no-change
#define non_stick 0
#define MASK 1
#define popinet_trick  0
#define DEGREE 0
#define REFERENCE 1

#define percent_60 0
#define percent_70 0
#define percent_80 0
#define percent_90 0
#define percent_95 0


#define case_1 0
#define case_2 0
#define case_3 0
#define case_4 0
#define case_5 0
#define case_6 0
#define case_7 0
#define case_8 0

// different change-ratio
#if case_1
#define ratio 1.72
#endif

#if case_2
#define ratio 1.45
#endif

#if case_3
#define ratio 1.33
#endif

#if case_4
#define ratio 1.24
#endif

#if case_5
#define ratio 1.12
#endif

#if case_6
#define ratio 0.9
#endif

#if case_7
#define ratio 0.81
#endif

#if case_8
#define ratio 0.69
#endif

/********************************************/
// fluid peroperties
/*
length_ratio(bubble tail to the bottom(4r0)+bubble body length(not including head)
This is compared with chengsi's experiment
*/

#if REFERENCE
#define RHOR 1000
#define MUR 100
#define N_f 4.2835144447
#define Mo 1.6e-2
#define Eo 9.8051217678  //pow(Mo*pow(N_f,4),1/3) 
#define length_ratio 10
#define ratio 1
#endif

#if percent_60
#define RHOR 910
#define MUR 1655
#define N_f 230.9
#define Eo 40.55

#if case_1
#define length_ratio 7.2 
#endif 
#if case_2
#define length_ratio 7.84 
#endif 
#if case_3
#define length_ratio 7.6 
#endif 
#if case_4
#define length_ratio 7.66 
#endif 
#if case_5
#define length_ratio 7.93 
#endif 
#if case_6
#define length_ratio 6.7 
#endif 
#if case_7
#define length_ratio 7.84 
#endif 
#if case_8
#define length_ratio 6.7 
#endif 

#endif

#if percent_70
#define RHOR 934
#define MUR 2906
#define N_f 134.93
#define Eo 41.86

#if case_1
#define length_ratio 7.38 
#endif 
#if case_2
#define length_ratio 7.6 
#endif 
#if case_3
#define length_ratio 8.87 
#endif 
#if case_4
#define length_ratio 7.1 
#endif 
#if case_5
#define length_ratio 8.1 
#endif 
#if case_6
#define length_ratio 10.15 
#endif 
#if case_7
#define length_ratio 9.42 
#endif 
#if case_8
#define length_ratio 6.27 
#endif 

#endif

#if percent_80
#define RHOR 1044
#define MUR 5350
#define N_f 81.91
#define Eo 47.52

#if case_1
#define length_ratio 7.9 
#endif 
#if case_2
#define length_ratio 8.07 
#endif 
#if case_3
#define length_ratio 8.69 
#endif 
#if case_4
#define length_ratio 10.36
#endif 
#if case_5
#define length_ratio 8.23 
#endif 
#if case_6
#define length_ratio 7.58 
#endif 
#if case_7
#define length_ratio 7.68 
#endif 
#if case_8
#define length_ratio 7.9 
#endif 
#endif

#if percent_90
#define RHOR 1067
#define MUR 16344
#define N_f 27.4
#define Eo 49.1

#if case_1
#define length_ratio 11.85 
#endif 
#if case_2
#define length_ratio 9.66 
#endif 
#if case_3
#define length_ratio 9.43 
#endif 
#if case_4
#define length_ratio 10.36 
#endif 
#if case_5
#define length_ratio 12.15 
#endif 
#if case_6
#define length_ratio 8.2 
#endif 
#if case_7
#define length_ratio 7.91 
#endif 
#if case_8
#define length_ratio 7.69 
#endif 
#endif

#if percent_95
#define RHOR 1075
#define MUR 29809
#define N_f 15.1
#define Eo 49.8

#if case_1
#define length_ratio 8.38 
#endif 
#if case_2
#define length_ratio 10.41 
#endif 
#if case_3
#define length_ratio 9.88
#endif 
#if case_4
#define length_ratio 8.38 
#endif 
#if case_5
#define length_ratio 10.74 
#endif 
#if case_6
#define length_ratio 9.88 
#endif 
#if case_7
#define length_ratio 12 
#endif 
#if case_8
#define length_ratio 10 
#endif 
#endif 



#if AXIS
  #include "axi.h" // fixme: does not run with -catch
#endif

#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
//#include "navier-stokes/conserving.h"
#include "tension.h"

#if HEAT
#include "diffusion.h"
#include "tracer.h"
#endif


#if REDUCED
  #include "reduced.h"
#endif

#include "vtk.h"
#include "view.h"
#include "output.h"

#include "navier-stokes/perfs.h"

/**********************************************************/

#if HEAT
scalar T[];
scalar *tracers = {T};
face vector D[];
#endif

#define Fr_0 1
#define _g 1 //in this case, consider U_0 = sqrt(gD_0)--->Fr_0^2 = U_0^2/(gD_0)

/* R_0 characteristic length (radias of lower tube) */
#define lambda_1  (2*(sqrt(1+2.44*pow(N_f,0.667))-1)/(2.44*pow(N_f,0.667))) // dimensionless liquid film thickness------refernce from chengsi's thesis tto guarantee the consistency with experiment
#define R_0 0.5
#define lambda  lambda_1*R_0 // liquid film
#define r_0 (R_0-lambda)
#define R_ex R_0*ratio
#define L_lower 8
#define L_upper 8
#define L0 L_lower+L_upper
/* grid infromation */
#define maxlevel 11
#define minlevel 9
#define delta_ratio (ratio - 1)
#define delta_length delta_ratio*R_0

#if DEGREE
#define degree 45
#define theta degree*pi/180
#define delta_x delta_length/tan(theta)
#endif  

#if !DEGREE
#define delta_x 1
#define theta arctan(delta_length/delta_x)  //gradual degree
#endif

/*-----------Froud numbers-------------------*/
/* predict the terminal velocity
double Froud(N_f, Eo){
    double m;
    double Fr;
    if(N_f<18)
        m = 25;
    else if(N_f>18&&N_f<258)
        m = 69*pow(N_f,-0.35);
    else
        m = 10;
    Fr = 0.345(1-pow(e,-0.01*N_f/0.345)(1-pow(e,(3.37-Eo)/m))); 
    return Fr;
}
// reference from "Numerical study of hydrodynamic characteristics of gas-liquid slug flow in vertical pipes"
*/


#if HEAT
#define T_low 293
#define T_top 293
#define T_other 293
#define T_taylor 293
#define T_liquid 293
#define D1 0.01
#define D2 0.1
#define D(f)  (clamp(f,0.,1.)*(D1 - D2) + D2)


mgstats mgd;
#endif


// boundary conditions
#if VERTICAL

bid top_lower;
bid top_upper;
scalar f1[];

uf.n[bottom] = 0.;
uf.n[top_lower] = 0.;
uf.n[top_upper] = 0.;
u.t[bottom] = neumann(0);
u.t[top_lower] = dirichlet(0);
//u.t[bottom] = dirichlet(0.35);
//u.t[top_lower] = dirichlet(0.35);
u.t[top_upper] = dirichlet(0);


// right side for outflow
u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
//u.n[right] = dirichlet(0.35);

#endif

#if popinet_trick
f1[top] = 0;
#endif

#if GRADUAL
bid top_lower;
bid top_upper;
bid gra;

uf.n[bottom] = 0.;
uf.n[top_lower] = 0.;
uf.n[top_upper] = 0.;
uf.n[gra] = 0.;
u.t[bottom] = dirichlet(0);
u.t[top_lower] = dirichlet(0);
//u.t[bottom] = dirichlet(0.35);
//u.t[top_lower] = dirichlet(0.35);
u.t[top_upper] = dirichlet(0);
f[right]   = neumann(0.); //liquid

u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
u.t[gra] = dirichlet(0);
//u.n[right] = dirichlet(0.35);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
#endif


#if no_change
uf.n[bottom] = 0.;
uf.n[top] = 0.;
u.t[bottom] = neumann(0);

u.t[top] = dirichlet(0);
//f[right]   = neumann(0.); //liquid
u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
#endif





#if HEAT
T[top_lower] = neumann(5);//heat flux
T[top_upper] = neumann(5);
T[bottom] = neumann(0);
T[left] = dirichlet(T_other);
T[right] = dirichlet(T_other);
#endif

#define TIME 30


// geometry of Taylor-bubble

double sphere(double x, double y, coord center, double radius) {

  return ( sq(x ) + sq (y ) - sq (radius));

}


double geometry(double x, double y) {

  coord center;
  foreach_dimension()
  center.x = 0;
  double s = sphere (x-length_ratio*r_0, y, center, r_0);
  double left = x - 4*r_0;
  double right = -x + length_ratio*r_0;
  double top = -y + r_0;
  double bottom=0;
  double c = -min(min(min(top,left),right),bottom);

  double geom = min(s, c);

  return geom;
}



int main() {
  size (L0);
  #if GRADUAL
  init_grid (1 << 9);
  #endif

  #if !GRADUAL
  init_grid(1<<minlevel);
  #endif
  
  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./N_f;
  mu2 = 1./(MUR*N_f);
  f.sigma = 1./Eo;
  DT = 1e-3;
  TOLERANCE = 1e-4;

  #if REDUCED
    G.x -= _g;
    Z.x = _g;
  #endif

  run();
}

void createFolder(const char* folder)
{
    char str[128];
    sprintf(str, "mkdir %s", folder); 
    system(str);   
}



#if GRADUAL
#define gradual_geo(x,y,R) (y-R)-(x-L_lower)*tan(theta)
#define dR_refine (2.*L0/(1 << minlevel))
#endif

event init (t = 0) {
// creat folders
createFolder("pressure");
createFolder("interface");
createFolder("cell");
createFolder("interfac-1");
createFolder("tip");
createFolder("gfsview");
createFolder("velocity");
createFolder("field");
createFolder("bview");

#if MASK
if (!restore (file = "dump")) {
#if no_change
mask(y>R_0?top:none);
#endif

#if VERTICAL
  mask((x<=L_lower&&y>R_0)?top_lower:
  (x>L_lower&&y>R_ex)?top_upper:
  none);
#endif

#if GRADUAL
  mask((x<=L_lower&&y>R_0)?top_lower:
  (L_lower<x&&x<=(L_lower+delta_x)&&(y-R_0)>(x-L_lower)*tan(theta))?gra:
  (x>(L_lower+delta_x)&&y>(R_ex))?top_upper: none);
  refine(x<=L_lower&&(y-R_ex<0)&&level<maxlevel);
  refine (L_lower<x&&x<=(L_lower+delta_x)&&(y>((R_0-dR_refine)+(x-L_lower)*tan(theta))) && (y<((R_0+dR_refine)+(x-L_lower)*tan(theta))) && level < maxlevel);
  refine(x>(L_lower+delta_x)&&(y-R_ex<0)&&level<maxlevel);
 #endif

 #if !GRADUAL
//refine ((x<=6&&(y-R_0<0)) && level < maxlevel);
  refine (y-R_ex<0 && level < maxlevel);
 #endif
 }
 #endif

 #if popinet_trick
  fraction (f1, difference(x,x<6? y>R_0 :y>R_ex));
  fraction (f, min(x<6? R_0 - y :R_ex - y,geometry(x,y))); // popinet's trick needs more grids due to consider entire domain.
  for (scalar s in {f1,f})
    s.refine = s.prolongation = fraction_refine;
  foreach() 
        u.x[] = 0;
  boundary({f,f1,u.x});
  #endif

  
  
  #if !popinet_trick
  //refine ((x<=6&&(y-R_0<0)) && level < maxlevel);
  fraction (f, geometry(x,y));
  #endif
  
}

//Popinet's trick
#if popinet_trick
event velocity (i++) {
  foreach()
    foreach_dimension()
      u.x[] = (1. - f1[])*u.x[];
  boundary ((scalar *){u});
}
#endif

#if HEAT
event init_T (i = 0) {
  foreach() {
	  T[] = T_liquid;
  }

  for (scalar f in interfaces){
	  foreach(){
		if ( f[] < 1e-6  ){ 	  
    	  T[] = T_taylor;
	    }
      }
   }
  boundary ({T});
}


// Time integration
event integration (i++) {

  foreach_face(){
	 D.x[] = D(f[]); 
  }
  
  mgd = diffusion (T, DT, D);
}
#endif


/**
We add the acceleration of gravity. */

#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= _g;
}
#endif

/**
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

/*
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (t = 0; t <= TIME; t += 0.01) {
  double xb = 0., yb = 0, vx = 0., vy = 0,  sb = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:vx)  reduction(+:vy) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv;
  }
  
  printf ("%g %g %g %g %g %g %g %g %g %g ", 
	  t, sb, -1., xb/sb, yb/sb, vx/sb, vy/sb, dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);

}


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-3,1e-2,1e-2}, maxlevel,minlevel);
  //adapt_wavelet ({f}, (double[]){5e-4}, LEVEL);
}
#endif

/*

void my_output_facets (struct My_OutputFacets p)
{
  scalar c = p.c;
  face vector s = p.s;
  if (!p.fp) p.fp = stdout;
  if (!s.x.i) s.x.i = -1;

  foreach()
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (c[], n);
#if dimension == 2      
      coord segment[2];
      if (facets (n, alpha, segment) == 2)
	fprintf (p.fp, "%g %g %g %g\n %g %g %g %g\n\n", 
		 x + segment[0].x*Delta, y + segment[0].y*Delta, 
		 interpolate(u.x, x + segment[0].x*Delta, y + segment[0].y*Delta),
		 interpolate(u.y, x + segment[0].x*Delta, y + segment[0].y*Delta),
		 x + segment[1].x*Delta, y + segment[1].y*Delta, 
		 interpolate(u.x, x + segment[1].x*Delta, y + segment[1].y*Delta),
		 interpolate(u.y, x + segment[1].x*Delta, y + segment[1].y*Delta));		 
		 //x + segment[1].x*Delta, y + segment[1].y*Delta);
#else // dimension == 3
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
	fprintf (p.fp, "%g %g %g\n",
		 x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
	fputc ('\n', p.fp);
#endif
    }

  fflush (p.fp);
}


*/


/******************************** output **********************************************/

event interface (t += 1) {

// interface
  char *outfile1 = NULL;
  outfile1 = (char *) malloc(sizeof(char) * 256);
  sprintf(outfile1, "interface/interface-%g", t);
  FILE * fp_interface = fopen (outfile1, "w");
  output_facets (f, fp_interface);
  fclose(fp_interface);

// cell 

  char *outfile = NULL;
  outfile = (char *) malloc(sizeof(char) * 1024);
  sprintf(outfile, "cell/cell-%g.dat", t);  
  FILE * fp = fopen (outfile, "w");
  output_cells (fp);
  fclose (fp);
}

/* This blocked the bubble go through the sudden changed pipe

event initial_interface (t += 0.1){
  scalar pos[];
  position (f, pos, {1,0});
  double max = statsf(pos).max;

if ((max-L_lower) >= 0){
  char *outfile1 = NULL;
  outfile1 = (char *) malloc(sizeof(char) * 256);
  sprintf(outfile1, "interface-1/interface-%g", t);
  FILE * fp_interface = fopen (outfile1, "w");
  
  output_facets (f, fp_interface);
  fclose(fp_interface);
}
}
*/

#if popinet_trick
event bview_output (t += 1) {
  char name[80];
  sprintf (name, "bview/dump-%d", (int) t);
  p.nodump = false;
  dump (file = name, list = {f,p,u.x,u.y,rho,mu});
}
#endif

// Delta h
event head (t += 0.01) {

  scalar pos[];
  position (f, pos, {1,0});
  double max = statsf(pos).max;

  //Point point = locate(11*r_0,0,0);
  
//interpolate(u.x, max,0),interpolate(u.y, max,0)
  char name[80];
  sprintf (name, "tip/head");
  static FILE * fp = fopen (name, "w");
  if ((max-6) >= 0)
    fprintf(fp,"%g %g\n",t,(max-6));
  fflush (fp);
}

// tip velocity
event tip (t += 0.01) {

  scalar pos[];
  position (f, pos, {1,0});
  double max = statsf(pos).max;

  //Point point = locate(11*r_0,0,0);
  
//interpolate(u.x, max,0),interpolate(u.y, max,0)
  char name[80];
  sprintf (name, "tip/tip");
  static FILE * fp = fopen (name, "w");
  fprintf(fp,"%g %g\n",t,max);
  fflush (fp);
}


//tail velocity
event tail (t += 0.01) {

  scalar pos[];
  position (f, pos, {1,0});
  double min = statsf(pos).min;

  //Point point = locate(11*r_0,0,0);
  
//interpolate(u.x, max,0),interpolate(u.y, max,0)
  char name[80];
  sprintf (name, "tip/tail");
  static FILE * fp = fopen (name, "w");
  fprintf(fp,"%g %g\n",t,min);
  fflush (fp);
}

// bubble length
event bubble_length (t += 0.01) {

  scalar pos_tip[], pos_tail[];
  position (f, pos_tip, {1,0});
  position (f, pos_tail, {1,0});
  double max = statsf(pos_tip).max;
  double min = statsf(pos_tail).min;
  //Point point = locate(11*r_0,0,0);
  char name[80];
  sprintf (name, "tip/length"); 
  static FILE * fp = fopen (name, "w");
  fprintf(fp,"%g %g \n",t,(max-min));
  fflush (fp);
  
}

// interface trajectory
event trajectory (t += 0.01) {

  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  double min = statsf(pos).min;
  //Point point = locate(11*r_0,0,0);
  
//interpolate(u.x, max,0),interpolate(u.y, max,0)
  char name[80];
  sprintf (name, "tip/trajectory");
  static FILE * fp = fopen (name, "w");
  fprintf(fp,"%g %g %g\n",t,max,min);
  fflush (fp);
}



// output velocity distribution
event velocity (t+=0.1){

  char name[100];
  sprintf (name, "velocity/vprof-%.4f", t);
  FILE*fp1 = fopen (name, "w");

  for (double y = 0; y <= R_ex; y += 0.02)
  for (double x = 0; x <= L0; x += 0.2)
    fprintf (fp1, " %g %g %g %g\n", x,y,
	    interpolate (u.x, x, y), interpolate (u.y, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
    fflush(fp1);
  }

// gfsview
event gfsview (t = 0; t <= TIME; t += 0.5)
{
  char name[80];
  sprintf (name, "gfsview/snapshot-%.4f.gfs",t);
  output_gfs (file = name, t = t);
}

// output vtk with field
event save_vtu_surface (t = 0; t <=  TIME; t += 1){
foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

char *name = NULL;
name = (char *) malloc(sizeof(char) * 256);
FILE *fp;
sprintf (name, "field/paraview--%d-%.6f.vtk",pid(),t);
fp = fopen (name, "w");
output_vtk ({f,p,rho,mu,u.x,u.y},N, fp, true);//draw streamlines in paraview-------(generate velocity vector first)
fclose(fp);

}

// pressure --- point at expansion center
event vprofile1 (t += 0.5){

  char name[80];
  sprintf (name, "pressure/pressure");
  FILE*fp = fopen (name, "w");
  fprintf(fp,"Grid"     "Time"      "pressure"     "Vx"      "Vy\n");
 
  float p_point = interpolate (p, 6, 0.1);
  float u_point = interpolate (u.x, 6, 0.1);
  float v_point = interpolate (u.y, 6, 0.1);

  fprintf(fp,"%g  %g  %g  %g  \n", 
  t, p_point, u_point,v_point);
  fprintf(ferr,"%g %g  %g  %g  \n", 
   t, p_point, u_point,v_point);

}

// movie
event movies (t = 0; t <= TIME; t += 0.05) {
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  char legend_time[80];
  sprintf(legend_time, "t = %g s", t);

  //char legend_parameter[80];
  //sprintf(legend_parameter, "Eo = 49.2\n Nf = 27.4\n Ratio = 1.24\n Mo = 2.1e-1");


  clear();
// interface
  view (tx = -0.5,ty=0);
  draw_vof ("f",fc = {1,0,0}, lw = 2);
  p.nodump = false;
  #if !HEAT
  squares("p",linear = true);
  #endif
  #if HEAT
 squares("T",linear = true);
  #endif
  draw_string(legend_time, 1, size = 30., lw = 2.);
  //draw_string(legend_parameter, 2, size = 30., {0.99,0.5,0.99});
  mirror({0,1}){
  draw_vof ("f",fc = {1,0,0}, lw = 2);
  p.nodump = false;
  #if !HEAT
  squares("p",linear = true);
  #endif
   #if HEAT
 squares("T",linear = true);
  #endif
    }

  save("interface.mp4");
  }
  
  /************************gnuplot*******************************/
/* -----Get the bubble shape

#!/bin/bash

for i in interface-*;do
echo '
reset
set size ratio -1
set terminal pdf
set output " '$i'.pdf"
unset tics
unset border
set object 1 rectangle from 6,-0.5*1.33 to 14,0.5*1.33 lw 1 fc rgb "light-blue"
set object 2 rectangle from 0,0.5 to 6,-0.5 lw 1 fc rgb "light-blue"
p [0:14][-0.5*1.33:0.5*1.33] "'$i'" u 1:2 w l lw 1 lc rgb "red" notitle,\
"'$i'" u 1:(-$2) w l lw 1 lc rgb "red" notitle
' | gnuplot
done

-------Get the bubble tip velocity
#!/bin/bash

awk '{print $1,$2-6}' tip>tip2

awk 'NR == 1{p = $2;next}
	{print $1, ($2 - p)/0.1; p = $2}
	END {print $1,p}' tip > tip-1

join tip-1 tip2 > tip-velocity

awk '{print $1,$2-6}' taill>taill2

awk 'NR == 1{p = $2;next}
	{print $1, ($2 - p)/0.1; p = $2}
	END {print $1,p}' taill > taill-1

join taill-1 taill2 > tail-velocity
*/

/************************run the code*******************************/
/*

First method: make $file.tst
This method linked the run commands in Basilisk, to generate one $file directory.
..........It should be noted that this doesn't work with 'mask' function mixed MPI, so if you want to use it with mask, just one processor.

Second method: tap the code by (Openmpi)
export OMP_NUM_THREADS=num (num means numbers of processors you want)
qcc -Wall -O3 -fopenmp $file.c -o $file  -L$BASILISK/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm //run this code
./taylor >> log 2> out  // log file and other output data file by function in Basilisk

*/

/* post-processing****************
1) Froude number evolution before and after change pipe
2) mass-center velocity
3) tip-velocity, tail-velocity
4) streamlines and velocity vector field
5) liquid film thickness 
6) bubble body length
*/


/************************python drawing pictures*******************************/
// cut experiment pictures
/*

from PIL import Image
import os.path, sys

path = "C:\\Users\\glk12\\Desktop\\nnn\\png\\80_3"
dirs = os.listdir(path)
path_new = "C:\\Users\\glk12\\Desktop\\nnn\\png\\80_3\\new"

def crop():
    for item in dirs:
        fullpath = os.path.join(path,item)         #corrected
        fullpath_new = os.path.join(path_new,item)
        if os.path.isfile(fullpath):
            im = Image.open(fullpath)
            f, e = os.path.splitext(fullpath_new)
            imCrop = im.crop((200, 537, 670, 713)) #corrected
            imCrop.save(f  + '.png', "PNG", quality=100)
            

crop()
*/

/*************streamlines&&velocity vectors***********************/

/*
face vector alpha[];
vector alpha_tmp[];
tensor alpha_new[];

foreach_face()
    alpha_tmp.x[] = (alpha.x[1] + alpha.x[]) / 2;
    boundary((scalar *){alpha_tmp});
foreach(){
    foreach_dimension()
        alpha_new.x.x[] = (alpha_tmp.x[1] - alpha_tmp.x[-1])/Delta;
        alpha_new.y.x[] = (alpha_tmp.y[1] - alpha_tmp.y[-1])/Delta;
        }
    boundary((vector *){alpha_new});
    */