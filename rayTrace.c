#include <FPT.h>
#include <D3d_matrix.h>

double degrees_of_half_angle;

int sphereEquation(double rayStart[3], double t[3]) {

}

int generateRay(int x, int y, double rgb[3]) {
  double tvals[3];
  double rayStart[3];
  double h;
  double nx, ny;
  h = tan(degrees_of_half_angle * M_PI / 180);

  nx = x - 300 * (h / 300.0);
  ny = y - 300 * (h / 300.0);

  tvals[0] = nx;
  tvals[1] = ny;
  tvals[2] = 1;

  rayStart[0] = 0; rayStart[1] = 0; rayStart[2] = 0;



}

int main() {
  degrees_of_half_angle = 45;

  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();
  int ix, iy;
  double rgb[3];
  for(ix = 0; ix < 600; ix++) {
    for(iy = 0; iy < 600; iy++) {
      generateRay(ix, iy, rgb);
      G_rgb(rgb[0], rgb[1], rgb[2]);
      G_point(ix, iy);
    }
  }
}
