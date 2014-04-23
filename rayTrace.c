#include <FPT.h>
#include <D3d_matrix.h>

double degrees_of_half_angle;
int numobjects;
int objectType[15];
double M[15][4][4], invM[15][4][4];
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
  int Tn;
  int Ttypelist[100];
  double Tvlist[100];
  AMBIENT = 0.2;
  MAX_DIFFUSE = 0.5;
  SPECPOW = 80;
  light_in_eye_space[0] = 100;
  light_in_eye_space[1] = 200;
  light_in_eye_space[2] = -50;

  eye[0] = 0; eye[1] = 0; eye[2] = 0;

  //degrees_of_half_angle = 45;

  // Object Type Indicators:
  // 0 = Sphere
  // 1 = Hyperboloid

  numobjects = 0;
  
  //Build the four pancakes
  objectType[numobjects] = 0;
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   10    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
  
  D3d_make_movement_sequence_matrix (M[numobjects], invM[numobjects],
				     Tn, Ttypelist, Tvlist) ;
  rgb[numobjects][0] = 0; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;

  numobjects++;

  objectType[numobjects] = 0;
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   10    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
  
  D3d_make_movement_sequence_matrix (M[numobjects], invM[numobjects],
				     Tn, Ttypelist, Tvlist) ;
  rgb[numobjects][0] = 0; rgb[numobjects][1] = 1; rgb[numobjects][2] = .25;
  
  numobjects++;
  
  objectType[numobjects] = 0;
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   10    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =   60    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
  
  D3d_make_movement_sequence_matrix (M[numobjects], invM[numobjects],
				     Tn, Ttypelist, Tvlist) ;
  rgb[numobjects][0] = 1; rgb[numobjects][1] = 1; rgb[numobjects][2] = 1;
  
  numobjects++;
  
  objectType[numobjects] = 0;
  Tn = 0 ; // number of transformations
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   10    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =  120    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
  
  D3d_make_movement_sequence_matrix (M[numobjects], invM[numobjects],
				     Tn, Ttypelist, Tvlist) ;
  rgb[numobjects][0] = 1; rgb[numobjects][1] = 0; rgb[numobjects][2] = 1;
  
  numobjects++;

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

int generateRay(int x, int y, double argb[3]) {
  double rayScreen[3], rayVec[15][3];
  double h;
  double nx, ny;
  int i, j;
  double eye_obj_space[15][3], ray_obj_space[15][3];
  double a, b, c;
  double intersects[15][2];
  double P[3], N[3];
  h = tan(degrees_of_half_angle * M_PI / 180);
  nx = (x - 300) * (h / 300.0);
  ny = (y - 300) * (h / 300.0);
  rayScreen[0] = nx;
  rayScreen[1] = ny;
  rayScreen[2] = 1;

  //find all intersections
  for(i = 0; i < numobjects; i++) {
    D3d_mat_mult_pt(eye_obj_space[i], invM[i], eye);
    D3d_mat_mult_pt(ray_obj_space[i], invM[i], rayScreen);
    //for(j = 0; j < 3; j++) {
    // printf("Eye obj: %lf, Ray Obj: %lf\n", eye_obj_space[i][j], ray_obj_space[i][j]);
    //}
    a = 0; b = 0; c = 0;
    for(j = 0; j < 3; j++) {
      rayVec[i][j] = ray_obj_space[i][j] - eye_obj_space[i][j];
      a += rayVec[i][j] * rayVec[i][j];
      b += eye_obj_space[i][j] * rayVec[i][j];
      c += eye_obj_space[i][j] * eye_obj_space[i][j];
    }
    b *= 2;
    c -= 1;
    quadEquation(intersects[i], a, b, c);
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
  double transposedM[4][4];
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      transposedM[i][j] = invM[objnum][j][i];
    }
    transposedM[3][i] = 0;
    transposedM[i][3] = 0;
  }
  transposedM[3][3] = 1;
  D3d_mat_mult_pt(N, transposedM, P);

  //Need to move P back to eyespace
  D3d_mat_mult_pt(P, M[objnum], P);

  nu_light_model(rgb[objnum], eye, P, N, argb);

}

int main() {
  
  init_scene();

  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();
  int ix, iy, i;
  double argb[3];
  for(i = 0; i <= 24; i++) {
    degrees_of_half_angle = 80 - 2*i;
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
