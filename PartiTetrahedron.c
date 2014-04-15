#include <FPT.h>
#include <D3d_matrix.h>

double light_in_eye_space[3];
double degrees_of_half_angle;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 80 ;

void path (int frame_number, double path_xyz[3])
{
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

int sphere(double p[3], double u, double v) {
  p[0] = cos(v)*cos(u);
  p[1] = sin(u);
  p[2] = sin(v)*cos(u);
}

void nu_light_model (double irgb[3],
                     double s[3],
                     double p[3],
                     double n[3],
                     double argb[3])
// s,p,n in eyespace
// irgb == inherent color of object (input to this function)
// s = location of start of ray
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// assume eye is at the origin
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]
{
  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;
  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;
  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;
  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below
  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ;
    NdotE *= (-1.0) ;    NdotL *= (-1.0) ;
  }
  // Blinn's variant
  double H[3] ;
  H[0] = E[0] + L[0] ;
  H[1] = E[1] + L[1] ;  
  H[2] = E[2] + L[2] ;
  len = sqrt(H[0]*H[0] + H[1]*H[1] + H[2]*H[2]) ;
  H[0] /= len ;  H[1] /= len ;  H[2] /= len ;
  double NdotH = N[0]*H[0] + N[1]*H[1] + N[2]*H[2] ;
  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }
  double specular ;
  if (NdotH <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(NdotH,SPECPOW) ;}
  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;
 LLL : ;
  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }
}

void draw(int (*func)(double[3], double, double), double Mat[15][4][4], 
	  double Imat[15][4][4], double rgb[15][3], int numobjects, 
	  double eye[3]) {
  double p[3], pu[3], pv[3], normvec[3];
  double u, v, h;
  int xtemp, ytemp;
  double xproj, yproj;
  double argb[3];
  int i;
  h = tan(degrees_of_half_angle * M_PI / 180);
  for(i = 0; i < numobjects; i++) {
    for(u = 0; u < M_PI*2; u += 0.005) {
      for(v = 0; v < M_PI*2; v += 0.005) {
	(*func)(p, u, v);
	(*func)(pu, u+0.001, v);
	(*func)(pv, u, v+0.001);
	D3d_mat_mult_pt(p, Mat[i], p);
	xproj = (300.0 / h) * (p[0] / p[2]) + 300;
	yproj = (300.0 / h) * (p[1] / p[2]) + 300;
	xtemp = (int)xproj;
	ytemp = (int)yproj;
	if(xtemp < 600 && ytemp < 600 && xtemp >= 0 && ytemp >= 0) {
	  if(z_buffer[xtemp][ytemp] > p[2]) {
	    z_buffer[xtemp][ytemp] = p[2];
	    D3d_mat_mult_pt(pu, Mat[i], pu);
	    D3d_mat_mult_pt(pv, Mat[i], pv);
	    findNormal(p, pu, pv, normvec);
	    nu_light_model(rgb[i], eye, p, normvec, argb);
	    G_rgb(argb[0], argb[1], argb[2]);
	    G_point(xtemp, ytemp);
	  }
	}
      }
    }
  }
}

int init_scene (int frame_number, double eye[3], double coi[3], double up[3],
		double Mat[15][3][3], double Imat[15][3][3])
{
  // model variables
  double xcen[4],ycen[4],zcen[4],brad ; // four nodes of tetrahedron
  double ccx,ccy,ccz,ccr ; // location of center of center sphere and radius

  int Tn;
  int Ttypelist[100];
  double Tvlist[100];


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


  //Path for eye and the up point for the eye.

  degrees_of_half_angle = 25 ;

  path (frame_number, eye) ;

  coi[0] = ccx ;
  coi[1] = ccy ;
  coi[2] = ccz ;

  path (frame_number + 1, up) ;

  printf("eye = %lf %lf %lf\n",eye[0],eye[1],eye[2]) ;
  printf("coi = %lf %lf %lf\n",coi[0],coi[1],coi[2]) ;
  printf("up  = %lf %lf %lf\n",up[0],up[1],up[2]) ;


  //Path for Light

  path (frame_number + 10, light_position) ;

  int i ;
  for(i = 0; i < 3; i++) {
    light_in_eye_space[i] = light_position[i];
  }

  AMBIENT  = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 80 ;

  //Make matrices for sphere positioning

  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[0]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[0]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[0]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;

  D3d_make_movement_sequence_matrix(Mat[0], Imat[0], Tn, Ttypelist, Tvlist);

  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[1]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[1]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[1]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;

  D3d_make_movement_sequence_matrix(Mat[1], Imat[1], Tn, Ttypelist, Tvlist);

  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[2]; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[2]; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[2]; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;

  D3d_make_movement_sequence_matrix(Mat[2], Imat[2], Tn, Ttypelist, Tvlist);

  Tn = 0;
  Ttypelist[Tn] = TX; Tvlist[Tn] = ccx; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ccy; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = ccz; Tn++;
  Ttypelist[Tn] = SX; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = ccr; Tn++;

  D3d_make_movement_sequence_matrix(Mat[3], Imat[3], Tn, Ttypelist, Tvlist);

}

int main() {
  double Mat[15][4][4], Imat[15][4][4];
  int numobjects;
  double eye[3], coi[3], up[3];

  double rgb[15][3]; 

  int frame_number;
  frame_number = 0;

  light_in_eye_space[0] = 100;
  light_in_eye_space[1] = 200;
  light_in_eye_space[2] = -50;

  init_scene(frame_number, eye, coi, up, Mat, Imat);

  numobjects = 5;
  int i;
  for(i = 0; i < numobjects; i++) {
    rgb[i][0] = 0; rgb[i][1] = 0; rgb[i][2] = 1;
  }

  draw(&sphere, Mat, Imat, rgb, numobjects, eye;
}
