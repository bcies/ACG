#include <FPT.h>
#include <D3d_matrix.h>

double Mat[15][4][4], Imat[15][4][4];
double Drawmat[15][4][4], Idrawmat[15][4][4];
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

int hyperboloid(double p[3], double u, double v) {
  //p[0] = 0.25*cosh(u)*cos(v);
  //p[1] = 0.25*cosh(u)*sin(v);
  //p[2] = 0.25*sinh(u);
  //p[0] = .25*sqrt(1+pow(u,2))*cos(v);
  //p[1] = .25*sqrt(1+pow(u,2))*sin(v);
  //p[2] = .25*u;

  double r = .25;
  double H = 1;
  double R = 1;

  double b = sqrt(r*r*H*H / (R*R - r*r));

  double y = u;
  double d = r/b*sqrt(y*y + b*b);

  double X = d* cos(v);
  double Z = d * sin(v);

  p[0] = X;
  p[1] = y;
  p[2] = Z;
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

int init_start() {
  // model variables
  double xcen[4], ycen[4], zcen[4], brad ; // four nodes of tetrahedron
  double hxcen[10], hycen[10], hzcen[10];
  double ccx, ccy, ccz, ccr ; // location of center of center sphere and radius

  int k;
  double theta;
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
  xcen[3] = 0 ; ycen[3] = sqrt(2) ; zcen[3] = 0 ;
  
  // also, figure out location of the 5th node of the model
  // which is at the center of mass of the tetrahedron
  ccx = 0 ; ccy = ycen[3] / 4 ; ccz = 0 ;
  
  //Finding the centers of all hyperboloids

  for(k = 0; k < 3; k++) {
    hxcen[k] = ((xcen[(k + 1)%3] - xcen[k])/2) + xcen[k];
    hycen[k] = 0;
    hzcen[k] = ((zcen[(k + 1)%3] - zcen[k])/2) + zcen[k];
    hxcen[k+3] = ((xcen[3]-xcen[k])/2) + xcen[k];
    hycen[k+3] = ((ycen[3]-ycen[k])/2) + ycen[k];
    hzcen[k+3] = ((zcen[3]-zcen[k])/2) + zcen[k];
  }
  //Middle, connecting Hyperboloids
  for(k = 6; k < 10; k++) {
    hxcen[k] = ((ccx - xcen[k - 6])/2) + xcen[k - 6];
    hycen[k] = ((ccy - ycen[k - 6])/2) + ycen[k - 6];
    hzcen[k] = ((ccz - zcen[k - 6])/2) + zcen[k - 6];
  }

  brad = 0.08 ; // radius of the 4 verts of the tetrahedron
  ccr  = 0.20 ; // the radius of the center node of the model
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  degrees_of_half_angle = 45 ;

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
  
  for(k = 0; k < 4; k++) {
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[k]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[k]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[k]; Tn++;
    
    D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				      Tn, Ttypelist, Tvlist);
    //Make the spheres blue
    rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
    numobjects++;
  }

  //Center Sphere
  Tn = 0;
  Ttypelist[Tn] = SX; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = TX; Tvlist[Tn] = ccx; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ccy; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = ccz; Tn++;
  
  D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects], 
				    Tn, Ttypelist, Tvlist);
  rgb[numobjects][0] = 1; rgb[numobjects][1] = 0; rgb[numobjects][2] = 0;
  numobjects++;
  int m;
  //Hyperboloids
  for(k = 0; k < 3; k++) {
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad*0.55; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad*10; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad*0.55; Tn++;
    Ttypelist[Tn] = RZ; Tvlist[Tn] = 90; Tn++;
    Ttypelist[Tn] = RY; 
    if(k == 0) { Tvlist[Tn] = 30; }
    else if(k == 1) { Tvlist[Tn] = 90; }
    else { Tvlist[Tn] = -30; }
    Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = hxcen[k]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = hycen[k]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = hzcen[k]; Tn++;
    
    D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects],
				      Tn, Ttypelist, Tvlist);
    rgb[numobjects][0] = 0.95; rgb[numobjects][1] = 0.95; 
    rgb[numobjects][2] = 0.95;
    numobjects++;
  }
  for(k = 0; k < 3; k++) {
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad*0.55; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad*10; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad*0.55; Tn++;
    Ttypelist[Tn] = RX;
    Tvlist[Tn] = -atan2(ycen[3], 1) * 180/M_PI + 20;
    printf("%lf\n", Tvlist[Tn]);
    Tn++;
    Ttypelist[Tn] = RY;
    if(k == 0) { Tvlist[Tn] = 90; }
    else if(k == 1) { Tvlist[Tn] = -30; }
    else { Tvlist[Tn] = 210; }
    Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = hxcen[k+3]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = hycen[k+3]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = hzcen[k+3]; Tn++;
    
    D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects],
				      Tn, Ttypelist, Tvlist);
    rgb[numobjects][0] = 0.95; rgb[numobjects][1] = 0.95; 
    rgb[numobjects][2] = 0.95;
    numobjects++;
  }
  //Inner hyperboloids
  for(k = 0; k < 4; k++) {
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad*.5; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = .5; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad*.5; Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = hxcen[k + 6]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = hycen[k + 6]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = hzcen[k + 6]; Tn++;
    
    D3d_make_movement_sequence_matrix(Mat[numobjects], Imat[numobjects],
				      Tn, Ttypelist, Tvlist);
    rgb[numobjects][0] = 0.95; rgb[numobjects][1] = 0.95; 
    rgb[numobjects][2] = 0.95;
    numobjects++;
  }

}

int init_scene (int frame_number) {
  int k;
  double ccx, ccy, ccz;
  double eye[3],coi[3],up[3] ;
  double light_position[3], amb, diff, spow ;
  ccx = 0 ; ccy = sqrt(2) / 4 ; ccz = 0 ;
  
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
 
  ///////////////////////////////////////////
  //TRANSLATE EYE TO ORIGIN AND ADJUST ALL OBJECTS 
  //AND LIGHT/REFERENCE PTS
  ///////////////////////////////////////////
  double view[4][4], view_inverse[4][4];
  
  if(!D3d_view(view, view_inverse, eye, coi, up)) {
    printf("Building view matrices went wrong.\n");
  }

  D3d_mat_mult_pt(light_in_eye_space, view, light_position);

  //Draw mat's are so I don't change the object matrices when I adjust for eyespace.
  for(k = 0; k < numobjects; k++) {
    D3d_copy_mat(Drawmat[k], Mat[k]);
    D3d_copy_mat(Idrawmat[k], Imat[k]);
    D3d_mat_mult(Drawmat[k], view, Drawmat[k]);
    D3d_mat_mult(Idrawmat[k], Idrawmat[k], view_inverse);
  } 

  //TODO make sure your top sphere and center sphere are positioned correctly
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

void findUnitVector(double vec[3], double returnvec[3]) {
  int i;
  double divider;
  divider = sqrt(pow(vec[0],2)+pow(vec[1],2)+pow(vec[2],2));
  for(i = 0; i < 3; i++) {
    returnvec[i] = vec[i]/divider;
  }
}

int findNormal(double p[3], double pu[3], double pv[3], double rvec[3]) {
  double vec1[3];
  double vec2[3];
  double nvec[3];
  vec1[0] = pv[0] - p[0];
  vec1[1] = pv[1] - p[1];
  vec1[2] = pv[2] - p[2];
  vec2[0] = pu[0] - p[0];
  vec2[1] = pu[1] - p[1];
  vec2[2] = pu[2] - p[2];
  D3d_x_product(nvec, vec1, vec2);
  findUnitVector(nvec, nvec);
  rvec[0] = nvec[0]; rvec[1] = nvec[1]; rvec[2] = nvec[2];
}

void draw(double matrix[15][4][4], 
	  double I_matrix[15][4][4], double rgbval[15][3], int numobjects) {
  double p[3], pu[3], pv[3], normvec[3];
  double u, v, h;
  int xtemp, ytemp;
  double xproj, yproj;
  double argb[3];
  int i;
  double eye[3];
  double start;
  double limit;
  eye[0] = 0; eye[1] = 0; eye[2] = 0;
  h = tan(degrees_of_half_angle * M_PI / 180);
  for(i = 0; i < numobjects; i++) {
    if(i < 5) {
      start = -M_PI / 2;
      limit = M_PI / 2;
    } else {
      start = -1;
      limit = 1;
    }
    for(u = start; u < limit; u += 0.005) {
      for(v = 0; v < M_PI*2; v += 0.005) {
	if(i < 5) {
	  sphere(p, u, v);
	  sphere(pu, u+0.001, v);
	  sphere(pv, u, v+0.001);
	} else {
	  hyperboloid(p, u, v);
	  hyperboloid(pu, u+0.001, v);
	  hyperboloid(pv, u, v+0.001);
	}
	
	D3d_mat_mult_pt(p, matrix[i], p);
	//printf("px: %lf py:%lf pz:%lf\n", p[0], p[1], p[2]);
	xproj = (300.0 / h) * (p[0] / p[2]) + 300;
	yproj = (300.0 / h) * (p[1] / p[2]) + 300;
	xtemp = (int)xproj;
	ytemp = (int)yproj;
	//printf("xtemp: %d, ytemp: %d, xproj: %lf, yproj: %lf\n", xtemp, ytemp, xproj, yproj);
	if(xtemp < 600 && ytemp < 600 && xtemp >= 0 && ytemp >= 0) {
	  if(z_buffer[xtemp][ytemp] > p[2]) {
	    //printf("Do I get here?\n");
	    z_buffer[xtemp][ytemp] = p[2];
	    D3d_mat_mult_pt(pu, matrix[i], pu);
	    D3d_mat_mult_pt(pv, matrix[i], pv);
	    findNormal(p, pu, pv, normvec);
	    nu_light_model(rgbval[i], eye, p, normvec, argb);
	    G_rgb(argb[0], argb[1], argb[2]);
	    G_point(xtemp, ytemp);
	  }
	}
      }
    }
  }
}

int main() {
  int frame_number;
  int i, j;

  init_start();

  G_init_graphics(600,600);

  for(frame_number = 0; frame_number < 40; frame_number++) {
    init_scene(frame_number);
    
    for(i = 0; i < 600; i++) {
      for(j = 0; j < 600; j++) {
	z_buffer[i][j] = 1000;
      }
    }
    
    
    G_rgb(0,0,0);
    G_clear();
    
    draw(Drawmat, Idrawmat, rgb, numobjects);
    
    G_wait_key();
  }
}
