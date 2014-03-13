#include <FPT.h>
#include <D3d_matrix.h>

int main() {
  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();
  int ix, iy;
  double r, g, b;
  for(ix = 0; ix < 600; ix++) {
    for(iy = 0; iy < 600; iy++) {
      
      G_rgb(r, g, b);
      G_point(ix, iy);
    }
  }

}
