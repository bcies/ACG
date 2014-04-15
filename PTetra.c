#include <FPT.h>
#include <D3d_matrix.h>

double Mat[15][4][4], Imat[15][4][4];
double rgb[15][3];
int numobjects;
double light_in_eye_space[3];
double AMBIENT;
double MAX_DIFFUSE;
double SPECPOW;
double degrees_of_half_angle;
double z_buffer[600][600];

int sphere(double p[3], double u, double v) {
  p[0] = cos(v)*cos(u);
  p[1] = sin(u);
  p[2] = sin(v)*cos(u);
}

void path (int frame_number, double path_xyz[3]) {
  double u,v,r ;
  double x,y,z ;
  u = 5*frame_number*M_PI/180 ;
  v = 0.3*u ;
  r = 2.0 + 1.4*sin(u) ;
  x = r*cos(u)*cos(v) ;
  y = r*sin(u) ;
  z = r*cos(u)*sin(v) ;
  path_xyz[0] = x ;
  path_xyz[1] = y ;
  path_xyz[2] = z ;
}

int init_scene (int frame_number) {
  // model variables
  double xcen[4],ycen[4],zcen[4],brad ; // four nodes of tetrahedron
  double ccx,ccy,ccz,ccr ; // location of center of center sphere and radius
  
  double eye[3],coi[3],up[3] ;
  double light_position[3], amb, diff, spow ;
  int k;
  double theta;
  double brad, ccr;
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  // build a ball and stick model of a tetrahedron
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  // 3 equally spaced pts around unit circle in the xz-plane 
  // form the base 
  for (k = 0 ; k < 3 ; k++) {
    theta = 2*M_PI*k/3 ;
    xcen[k] = cos(theta) ;
    ycen[k] = 0 ;
    zcen[k] = sin(theta) ;
  }
  
  // you figure where the 4th node of the regular tetrahedron
  xcen[3] = 0 ; ycen[3] = 1.5 ; zcen[3] = 0 ;
  
  // also, figure out location of the 5th node of the model
  // which is at the center of mass of the tetrahedron
  ccx = 0 ; ccy = 1 ; ccz = 0 ;
  
  brad = 0.08 ; // radius of the 4 verts of the tetrahedron
  ccr  = 0.20 ; // the radius of the center node of the model
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  degrees_of_half_angle = 25 ;
  
  path (frame_number, eye) ;
  
  coi[0] = ccx ;
  coi[1] = ccy ;
  coi[2] = ccz ;
  
  path (frame_number + 1, up) ;
  
  printf("eye = %lf %lf %lf\n",eye[0],eye[1],eye[2]) ;
  printf("coi = %lf %lf %lf\n",coi[0],coi[1],coi[2]) ;
  printf("up  = %lf %lf %lf\n",up[0],up[1],up[2]) ;
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  path (frame_number + 10, light_position) ;
  AMBIENT  = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 80 ;
  //////////////////////////////////////////////
  // SET UP MATRICES
  //////////////////////////////////////////////
  int Tn;
  int Ttypelist[100];
  double Tvlist[100];
  
  numobjects = 0;
  
  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[0]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[0]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[0]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;
  
  D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				    Tn, Ttypelist, Tvlist);
  //Make the spheres blue
  rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
  numobjects++;
  
  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[1]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[1]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[1]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;
  
  D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				    Tn, Ttypelist, Tvlist);
  //Make the spheres blue
  rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
  numobjects++;

  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[2]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[2]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[2]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;
  
  D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				    Tn, Ttypelist, Tvlist);
  //Make the spheres blue
  rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
  numobjects++;
  
  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[3]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[3]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[3]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;
  
  D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				    Tn, Ttypelist, Tvlist);
  //Make the spheres blue
  rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
  numobjects++;

  //Center Sphere
  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = ccx; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ccy; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = ccz; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = ccr; Tn++;
  
  D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				    Tn, Ttypelist, Tvlist);
  rgb[numobjects][0] = 1; rgb[numobjects][1] = 0; rgb[numobjects][2] = 0;
  numobjects++;
  ///////////////////////////////////////////
  //TRANSLATE EYE TO ORIGIN AND ADJUST ALL OBJECTS 
  //AND LIGHT/REFERENCE PTS
  ///////////////////////////////////////////
  double view[4][4], view_inverse[4][4];
  
  if(!D3d_view(view, view_inverse, eye, coi, up)) {
    printf("Building view matrices went wrong.\n");
  }

  for(k = 0; k < numobjects; k++) {
    D3d_mat_mult(Mat[numobjects], view, Mat[numobjects]);
    D3d_mat_mult(Imat[numobjects], Imat[numobjects], view_inverse);
  } 

  //TODO make sure your top sphere and center sphere are positioned correctly
}
