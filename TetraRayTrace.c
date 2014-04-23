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
