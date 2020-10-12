/*
The output file for 2D cylinderical coordinate Taylor bubble rising, as follows:
1). bubble interface and cell evolution with time
2). the distance between bubble nose and expansion after bubble through expansion
3). coordinate of bubble nose and tail evolution with time
4). bubble frontal radius 
5). velocity field 
6). wall shear stress
7). pressure ---- middle liquid film (for bubble body)
8). gfsview file
9). output_vtk (for other scalar or vector field)
*/


/******************************** output **********************************************/
#include "fracface.h"


event interface (t += 0.1) {
    p.nodump = false;

// interface
    char *outfile1 = NULL;
    outfile1 = (char *) malloc(sizeof(char) * 256);
    sprintf(outfile1, "interface/interface-%.1f",t);
    FILE * fp1 = fopen (outfile1, "w");
    output_facets (f, fp1);
    fclose(fp1);
// cell 

    char *outfile = NULL;
    outfile = (char *) malloc(sizeof(char) * 1024);
    sprintf(outfile, "interface/cell-%.1f",t);  
    FILE * fp = fopen (outfile, "w");
    output_cells (fp);
    fclose (fp);
}


// liquid film
event liquid (t += 0.01) {
    scalar pos[];
    position (f, pos, {1,0},{8,0});
    double max = statsf(pos).max;
    double min = statsf(pos).min;
    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "tip/liquid");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g\n", t, max,min);
    
    fflush (fp);
}

// frontal radius
event front_radius(t += 0.1){
    vector h[];
    heights (f, h); // return height 
    
    scalar pos[];
    position (f, pos, {1,0});
    double max = statsf(pos).max;
    scalar pos1[];
    position (f, pos1, {0,1});
    double max_1 = statsf(pos1).max;

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "curvature/curvature");
    static FILE * fp = fopen (name, "w");
    Point point = locate (max, 0);
    double a = height_curvature(point, f, h); // return kappa
    fprintf (fp, "%g %g %g %g %g\n", t, 1/fabs(a),max,max_1, 2/(fabs(a)*max_1));
    fflush (fp);
}




// output velocity distribution
event velocity (t+=0.1){

    char name[1000];
    sprintf (name, "velocity/vprof-%.1f",t);
    FILE*fp1 = fopen (name, "w");


    for (double y = 0; y < R_ex; y += 0.02)
    for (double x = 0; x < L0; x += 0.1)
        fprintf (fp1, " %g %g %g %g %g\n", t, x,y,
            interpolate (u.x, x, y), interpolate (u.y, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
        fflush(fp1);
  }


event velocity_center (t+=0.1){

    char name[1000];
    sprintf (name, "velocity/vcenter-%.1f",t);
    FILE*fp1 = fopen (name, "w");

    double vx = 0., vy = 0,  sb = 0.;
  foreach(reduction(+:vx)  reduction(+:vy) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv;
  }

vector uu[];
foreach_face(){
  uu.x[] = u.x[] - vx/sb*fm.x[];
  uu.y[] = u.y[] - vy/sb*fm.y[];
}
boundary((scalar *){uu});


    for (double y = 0; y < R_0; y += 0.02)
    for (double x = 0; x < L0; x += 0.1)
        fprintf (fp1, " %g %g %g %g %g\n", t, x,y,
            interpolate (uu.x, x, y), interpolate (uu.y, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
        fflush(fp1);
  }


// gfsview
event gfsview (t = 0; t <= TIME; t += 0.5)
{
    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    FILE *fp;
    sprintf (name, "gfsview/snapshot-%.1f.gfs",t);
    fp = fopen (name, "w");
    output_gfs (fp);
}



event shear (t+=0.1){

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "shear/shear-%.1f",t);
  FILE*fp1 = fopen (name, "w");

  scalar shear[];
  foreach()
    shear[] = mu.y[]*(u.x[0,1] - u.x[0,-1])/(2*Delta);
  boundary ({shear});

  for (double y = (R_ex-0.002); y < R_ex; y += 0.002)
  for (double x = 0; x <= L0; x += 0.1)
    fprintf (fp1, " %g %g %g %g %g\n", x,y,
	    interpolate (u.x, x, y), interpolate (u.y, x, y),
      interpolate (shear, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
    fflush(fp1);
  }


// pressure --- middle line of liquid film

event pressure (t+=0.1){

scalar pos[];
    position (f, pos, {0,1});
    double max = statsf(pos).max;

    scalar pos_1[];
    position (f, pos_1, {1,0});
    double max_1 = statsf(pos_1).max;
    double min_1 = statsf(pos_1).min;

    char name[1000];
    sprintf (name, "pressure/pressure-%.1f",t);
    FILE*fp = fopen (name, "w");
    
    for (double y = (max+(R_0 - max)/2); y < (max+(R_0 - max)/2) + 0.0001; y += 0.0001)
    for (double x = min_1; x <= max_1; x += 0.1)
        fprintf (fp, " %g %g %g\n", x,y,interpolate (p, x, y)); 
    fflush(fp);
  }



event mid_expansion(t+=0.1){
    scalar shear[];
    foreach()
        shear[] = mu.y[]*(u.x[0,1] - u.x[0,-1])/(Delta);
    boundary ({shear});
    char name[1000];
    sprintf (name, "pressure/mid-prop");
    FILE*fp = fopen (name, "a");
    fprintf (fp, " %g %g %g %g\n", t,interpolate (p, 8, 0), interpolate(u.x,8,0),interpolate(shear,8,0));
    fflush(fp);


}