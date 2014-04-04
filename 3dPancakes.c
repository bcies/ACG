#include <FPT.h>
#include <D3d_matrix.h>

int main() {
  double Mat[4][4], Imat[4][4];
  int Tn;
  int TtypeList[100];
  double Tvlist[100];
  double eye[3], coi[3], up[3];
  double z_buffer[600][600];

  // eyespace 
  degrees_of_half_angle = 80-2*frame_number ;
  
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
  pos[0] = 100 ; pos[1] = 200 ; pos[2] = -50 ;
  amb  = 0.2 ;
  diff = 0.5 ;
  spow = 80 ;

  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();

  //Build the four pancakes.

}
