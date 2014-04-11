#include <FPT.h>
#include <D3d_matrix.h>

double z_buffer[600][600];
double light_in_eye_space[3];
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 80 ;
double degrees_of_half_angle;

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

void crossProduct(double v1[3], double v2[3], double returnvec[3]) {
  double result[3];
  int i;
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  for(i = 0; i < 3; i++) {
    returnvec[i] = result[i];
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
  crossProduct(vec1, vec2, nvec);
  findUnitVector(nvec, nvec);
  rvec[0] = nvec[0]; rvec[1] = nvec[1]; rvec[2] = nvec[2];
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

int main() {
  double Mat[15][4][4], Imat[15][4][4];
  int numobjects;
  numobjects = 0;
  int Tn;
  int Ttypelist[100];
  double Tvlist[100];
  double eye[3], coi[3], up[3];
  int frame_number;
  double rgb[15][3];
  int i, j;
  
  eye[0] = 0 ;
  eye[1] = 0 ;
  eye[2] = 0 ;
  
  coi[0] = 0 ;
  coi[1] = 0 ;
  coi[2] = 1 ;
  
  up[0] =  0 ;
  up[1] =  1 ; 
  up[2] =  0 ;
  
  // position of light
  light_in_eye_space[0] = 100 ; 
  light_in_eye_space[1] = 200 ; 
  light_in_eye_space[2] = -50 ;
  //amb  = 0.2 ;
  //diff = 0.5 ;
  //spow = 80 ;
  

    
    //Build the four pancakes.
    Tn = 0 ; 
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   10    ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
    
    D3d_make_movement_sequence_matrix (Mat[numobjects], Imat[numobjects],
				       Tn, Ttypelist, Tvlist) ;
    // color 0.00, 0.00, 1.00
    rgb[numobjects][0] = 0;
    rgb[numobjects][1] = 0;
    rgb[numobjects][2] = 1;
    numobjects++;
    
    
    Tn = 0 ; // number of transformations
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   10    ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
    
    D3d_make_movement_sequence_matrix (
				       Mat[numobjects],
				       Imat[numobjects],
				       Tn,
				       Ttypelist,
				       Tvlist) ;
    
    // color 0.00, 1.00, 0.25
    rgb[numobjects][0] = 0;
    rgb[numobjects][1] = 1;
    rgb[numobjects][2] = .25;
    
    numobjects++;
    
    Tn = 0 ; // number of transformations
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   10    ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = RY ; Tvlist[Tn] =   60    ; Tn++ ;
    Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
    
    D3d_make_movement_sequence_matrix (
				       Mat[numobjects],
				       Imat[numobjects],
				       Tn,
				       Ttypelist,
				       Tvlist) ;
    
    // color 1.00, 1.00, 1.00
    rgb[numobjects][0] = 1;
    rgb[numobjects][1] = 1;
    rgb[numobjects][2] = 1;
    
    numobjects++;
    
    Tn = 0 ; // number of transformations
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   10    ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  100    ; Tn++ ;
    Ttypelist[Tn] = RY ; Tvlist[Tn] =  120    ; Tn++ ;
    Ttypelist[Tn] = RX ; Tvlist[Tn] =  -70    ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  400    ; Tn++ ;
    
    D3d_make_movement_sequence_matrix (
				       Mat[numobjects],
				       Imat[numobjects],
				       Tn,
				       Ttypelist,
				       Tvlist) ;
    
    
    // color 1.00, 0.00, 1.00
    rgb[numobjects][0] = 1;
    rgb[numobjects][1] = 0;
    rgb[numobjects][2] = 1;
    
    numobjects++;
    
  G_init_graphics(600,600);
  
  for(frame_number = 0; frame_number <= 24; frame_number++) { 
    printf("Frame %d\n", frame_number);
    for(i = 0; i < 600; i++) {
      for(j = 0; j < 600; j++) {
	z_buffer[i][j] = 1000;
      }
    }
    
    // eyespace 
    degrees_of_half_angle = 80-2*frame_number ;
    G_rgb(0,0,0);
    G_clear();

    draw(&sphere, Mat, Imat, rgb, numobjects, eye);
    //G_save_image_to_file("bpancakeimg0001.xwd");
    G_wait_key();
  }

}
