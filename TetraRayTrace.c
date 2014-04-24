#include <FPT.h>
#include <D3d_matrix.h>

double degrees_of_half_angle;
int numobjects;
int objectType[15];
double M[15][4][4], invM[15][4][4];
double dM[15][4][4], dinvM[15][4][4];
double rgb[15][3];
double eye[3], light_in_eye_space[3];
double AMBIENT;
double MAX_DIFFUSE;
double SPECPOW;

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

int init_scene() {
  double xcen[4], ycen[4], zcen[4], brad ;
  double hxcen[10], hycen[10], hzcen[10];
  double ccx, ccy, ccz, ccr ;
  int k;
  double theta;
  int Tn;
  int Ttypelist[100];
  double Tvlist[100];
  AMBIENT = 0.2;
  MAX_DIFFUSE = 0.5;
  SPECPOW = 80;
  degrees_of_half_angle = 25;

  //Centers of all spheres in scene.
  for (k = 0 ; k < 3 ; k++) {
    theta = 2*M_PI*k/3 ;
    xcen[k] = cos(theta) ;
    ycen[k] = 0 ;
    zcen[k] = sin(theta) ;
  }
  xcen[3] = 0 ; ycen[3] = sqrt(2) ; zcen[3] = 0 ;
  ccx = 0 ; ccy = ycen[3] / 4 ; ccz = 0 ;
 
  //Centers of all hyperboloids
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
  
  //Making the movement matrices
  // Object Type Indicators:
  // 0 = Sphere
  // 1 = Hyperboloid

  numobjects = 0; 
  for(k = 0; k < 4; k++) {
    objectType[numobjects] = 0;
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad; Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = xcen[k]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = ycen[k]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = zcen[k]; Tn++;
    D3d_make_movement_sequence_matrix(M[numobjects], invM[numobjects], 
				      Tn, Ttypelist, Tvlist);
    //Make the spheres blue
    rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
    numobjects++;
  }
  //Center Sphere
  objectType[numobjects] = 0;
  Tn = 0;
  Ttypelist[Tn] = SX; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SY; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = SZ; Tvlist[Tn] = ccr; Tn++;
  Ttypelist[Tn] = TX; Tvlist[Tn] = ccx; Tn++;
  Ttypelist[Tn] = TY; Tvlist[Tn] = ccy; Tn++;
  Ttypelist[Tn] = TZ; Tvlist[Tn] = ccz; Tn++;
  D3d_make_movement_sequence_matrix(M[numobjects], invM[numobjects], 
				    Tn, Ttypelist, Tvlist);
  rgb[numobjects][0] = 1; rgb[numobjects][1] = 0; rgb[numobjects][2] = 0;
  numobjects++;
  //Hyperboloids
  for(k = 0; k < 3; k++) {
    objectType[numobjects] = 1;
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad*.25; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad*.25; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad*5; Tn++;
    Ttypelist[Tn] = RY; 
    if(k == 0) { Tvlist[Tn] = -60; }
    else if(k == 1) { Tvlist[Tn] = 0; }
    else { Tvlist[Tn] = 60; }
    Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = hxcen[k]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = hycen[k]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = hzcen[k]; Tn++; 
    D3d_make_movement_sequence_matrix(M[numobjects], invM[numobjects],
				      Tn, Ttypelist, Tvlist);
    rgb[numobjects][0] = 0.95; rgb[numobjects][1] = 0.95; 
    rgb[numobjects][2] = 0.95;
    numobjects++;
  }
  for(k = 0; k < 3; k++) {
    objectType[numobjects] = 1;
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad*.25; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad*.25; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad*5; Tn++;
    Ttypelist[Tn] = RX;
    Tvlist[Tn] = atan2(ycen[3], 1) * 180/M_PI;
    //printf("%lf\n", Tvlist[Tn]);
    Tn++;
    Ttypelist[Tn] = RY;
    if(k == 0) { Tvlist[Tn] = 90; }
    else if(k == 1) { Tvlist[Tn] = -30; }
    else { Tvlist[Tn] = 210; }
    Tn++;
    Ttypelist[Tn] = TX; Tvlist[Tn] = hxcen[k+3]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = hycen[k+3]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = hzcen[k+3]; Tn++;
    D3d_make_movement_sequence_matrix(M[numobjects], invM[numobjects],
				      Tn, Ttypelist, Tvlist);
    rgb[numobjects][0] = 0.95; rgb[numobjects][1] = 0.95; 
    rgb[numobjects][2] = 0.95;
    numobjects++;
  }
  //Inner hyperboloids
  for(k = 0; k < 4; k++) {
    objectType[numobjects] = 1;
    Tn = 0;
    Ttypelist[Tn] = SX; Tvlist[Tn] = brad*.2; Tn++;
    Ttypelist[Tn] = SY; Tvlist[Tn] = brad*.2; Tn++;
    Ttypelist[Tn] = SZ; Tvlist[Tn] = brad*2.75; Tn++;
    Ttypelist[Tn] = RX;
    if(k == 3) {
      Tvlist[Tn] = 90;
    } else {
      Tvlist[Tn] = 20;
    }
    Tn++;
    if(k != 3) {
      Ttypelist[Tn] = RY;
      if(k == 0) { Tvlist[Tn] = 90; }
      else if (k == 1) { Tvlist[Tn] = -30; }
      else { Tvlist[Tn] = 210; }
      Tn++;
    }
    Ttypelist[Tn] = TX; Tvlist[Tn] = hxcen[k + 6]; Tn++;
    Ttypelist[Tn] = TY; Tvlist[Tn] = hycen[k + 6]; Tn++;
    Ttypelist[Tn] = TZ; Tvlist[Tn] = hzcen[k + 6]; Tn++;
    
    D3d_make_movement_sequence_matrix(M[numobjects], invM[numobjects],
				      Tn, Ttypelist, Tvlist);
    rgb[numobjects][0] = 0.95; rgb[numobjects][1] = 0.95; 
    rgb[numobjects][2] = 0.95;
    numobjects++;
  }
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

int makePaths(int frame_number) {
  int k;
  double ccx, ccy, ccz;
  double tempeye[3],coi[3],up[3] ;
  double light_position[3], amb, diff, spow ;
  ccx = 0 ; ccy = sqrt(2) / 4 ; ccz = 0 ;
  path (frame_number, tempeye) ;
  coi[0] = ccx ;
  coi[1] = ccy ;
  coi[2] = ccz ;
  path (frame_number + 1, up) ;
  printf("eye = %lf %lf %lf\n",tempeye[0],tempeye[1],tempeye[2]) ;
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
  if(!D3d_view(view, view_inverse, tempeye, coi, up)) {
    printf("Building view matrices went wrong.\n");
  }
  D3d_mat_mult_pt(light_in_eye_space, view, light_position);
  //Draw mat's are so I don't change the object matrices when I adjust for eyespace.
  for(k = 0; k < numobjects; k++) {
    D3d_copy_mat(dM[k], M[k]);
    D3d_copy_mat(dinvM[k], invM[k]);
    D3d_mat_mult(dM[k], view, dM[k]);
    D3d_mat_mult(dinvM[k], dinvM[k], view_inverse);
  } 
}

int quadEquation(double intersects[2], double a, double b, double c) {
  double discriminant;
  discriminant = b*b - 4*a*c;
  if(a == 0) {
    if(b == 0) {
      intersects[0] = -1;
      intersects[1] = -1;
    } else if(c == 0) {
      intersects[0] = 0;
      intersects[1] = -1;
    } else {
      intersects[0] = -c / b;
      intersects[1] = -1;
    }
  } else if(discriminant < 0) {
    intersects[0] = -1;
    intersects[1] = -1;
  } else {
    intersects[0] = (-b - sqrt(discriminant))/(2*a);
    intersects[1] = (-b + sqrt(discriminant))/(2*a);
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

double dotProduct(double v1[3], double v2[3]) {
  double product;
  int i;
  product = 0;
  for(i = 0; i < 3; i++) {
    product += v1[i] * v2[i];
  }
  return product;
}

int generateRay(int x, int y, double argb[3]) {
  double rayScreen[3], rayVec[15][3];
  double h;
  double nx, ny;
  int i, j;
  double eye_obj_space[15][3], ray_obj_space[15][3];
  double a, b, c;
  double intersects[15][2];
  double P[3], N[3];
  double testP;
  h = tan(degrees_of_half_angle * M_PI / 180);
  nx = (x - 300) * (h / 300.0);
  ny = (y - 300) * (h / 300.0);
  rayScreen[0] = nx;
  rayScreen[1] = ny;
  rayScreen[2] = 1;

  //find all intersections
  for(i = 0; i < numobjects; i++) {
    D3d_mat_mult_pt(eye_obj_space[i], dinvM[i], eye);
    D3d_mat_mult_pt(ray_obj_space[i], dinvM[i], rayScreen);
    a = 0; b = 0; c = 0;
    if(objectType[i] == 0) {
      for(j = 0; j < 3; j++) {
	rayVec[i][j] = ray_obj_space[i][j] - eye_obj_space[i][j];
	a += rayVec[i][j] * rayVec[i][j];
	b += eye_obj_space[i][j] * rayVec[i][j];
	c += eye_obj_space[i][j] * eye_obj_space[i][j];
      }
      b *= 2;
      c -= 1;
    } else if(objectType[i] == 1) {
      for(j = 0; j < 2; j++) {
	rayVec[i][j] = ray_obj_space[i][j] - eye_obj_space[i][j];
	a -= rayVec[i][j] * rayVec[i][j];
	b -= eye_obj_space[i][j] * rayVec[i][j];
	c -= eye_obj_space[i][j] * eye_obj_space[i][j];
      }
      rayVec[i][2] = ray_obj_space[i][2] - eye_obj_space[i][2];
      a += rayVec[i][j] * rayVec[i][j];
      b += eye_obj_space[i][2] * rayVec[i][2];
      c += eye_obj_space[i][2] * eye_obj_space[i][2];
      b *= 2;
      c += 1;
    }
    quadEquation(intersects[i], a, b, c);
    //Cutting hyperboloid
    for(j = 0; j < 2; j++) {
      if(intersects[i][j] >= 0) {
    	testP = eye_obj_space[i][2] + rayVec[i][2] * intersects[i][j];
	if(abs(testP) > 1) {
	  intersects[i][j] = -1;
	}
     }
    }
  }

  //find closest intersect
  int objnum = -1;
  double closeT = 1000;
  for(i = 0; i < numobjects; i++) {
    if(intersects[i][0] >= 0 && intersects[i][0] < closeT) {
      objnum = i;
      closeT = intersects[i][0];
    }
    if(intersects[i][1] >= 0 && intersects[i][1] < closeT) {
      objnum = i;
      closeT = intersects[i][1];
    }
  }

  //if(objnum != -1) {
  // printf("x: %d, y: %d\n", x, y);
  // printf("nx: %lf, ny: %lf\n", nx, ny);
  // printf("Objnum is %d.\n", objnum);
  //  printf("T is %lf\n", closeT);
  //}
  for(i = 0; i < 3; i++) {
    P[i] = eye_obj_space[objnum][i] + rayVec[objnum][i] * closeT;
  }

  //find normal
  if(objectType[objnum] == 0) {
    double transposedM[4][4];
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
	transposedM[i][j] = dinvM[objnum][j][i];
      }
      transposedM[3][i] = 0;
      transposedM[i][3] = 0;
    }
    transposedM[3][3] = 1;
    D3d_mat_mult_pt(N, transposedM, P);
  } else if(objectType[objnum] == 1) {
    for(i = 0; i < 2; i++) {
      N[i] = -2*P[i];
    }
    N[2] = 2*P[i];
    findUnitVector(N, N);
    if(dotProduct(N, P) < 0) {
      for(i = 0; i < 3; i++) {
	N[i] = -N[i];
      }
    }
  }
  
  //Need to move P back to eyespace
  D3d_mat_mult_pt(P, dM[objnum], P);

  nu_light_model(rgb[objnum], eye, P, N, argb);

}

int main() {
  
  init_scene();

  //This is just because I never actually adjust eye,
  //Only use a temp eye in makePaths, which adjusts
  //everything so that it's in eyespace. (Eye at origin.)
  eye[0] = 0; eye[1] = 0; eye[2] = 0;

  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();
  int ix, iy, i;
  double argb[3];
  for(i = 0; i <= 40; i++) {
    makePaths(i);
    for(ix = 0; ix < 600; ix++) {
      for(iy = 0; iy < 600; iy++) {
	generateRay(ix, iy, argb);
	G_rgb(argb[0], argb[1], argb[2]);
	G_point(ix, iy);
      }
    }
    G_display_image();
    G_wait_key();
  }
}
