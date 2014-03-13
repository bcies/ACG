#include <FPT.h>
#include <D3d_matrix.h>

double P_INCREMENT = 0.001;
  
int circle (double p[2], double u) {
  p[0] = cos(u); 
  p[1] = sin(u);
}

int  sum4 (double p[2], double u) {
  if(u <= M_PI * 0.5) {
    p[0] = sqrt(fabs(cos(u)));
    p[1] = sqrt(fabs(sin(u))); 
  } else if(u <= M_PI) {
    p[0] = -sqrt(fabs(cos(u)));
    p[1] = sqrt(fabs(sin(u)));
  } else if(u <= 1.5 * M_PI) {
    p[0] = -sqrt(fabs(cos(u)));
    p[1] = -sqrt(fabs(sin(u)));
  } else if(u <= 2 * M_PI) {
    p[0] = sqrt(fabs(cos(u)));
    p[1] = -sqrt(fabs(sin(u)));
  }
}

  

int square (double p[2], double u) {
	if(u <= M_PI * 0.5) {
		p[0] = pow(cos(u), 2);
 		p[1] = pow(sin(u), 2);
	} else if(u <= M_PI) {
		p[0] = -pow(cos(u), 2);
		p[1] = pow(sin(u), 2);
	} else if(u <= 1.5 * M_PI) {
		p[0] = -pow(cos(u), 2);
 		p[1] = -pow(sin(u), 2);
	} else if(u <= 2 * M_PI) {
		p[0] = pow(cos(u), 2);
 		p[1] = -pow(sin(u), 2);
	}
}

int asteroid (double p[2], double u) { 
  //sqrt(|x|) + sqrt(|y|) = 1  
  if(u <= M_PI * 0.5) {
  	p[0] = pow(cos(u), 4);
  	p[1] = pow(sin(u), 4);
  } else if(u <= M_PI) {
  	p[0] = -pow(cos(u), 4);
  	p[1] = pow(sin(u), 4);   
  } else if(u <= 1.5 * M_PI) {
   	p[0] = -pow(cos(u), 4);
  	p[1] = -pow(sin(u), 4);   
  } else if(u <= 2 * M_PI) {
   	p[0] = pow(cos(u), 4);
  	p[1] = -pow(sin(u), 4);   
  }
 }


int hyperbola (double p[2], double u) {
  p[0] = cosh(u);
  p[1] = sinh(u);
  //Faces right.
}

int parabola (double p[2], double u) {
  p[0] = u;
  p[1] = pow(u, 2);
  //Faces up.
}

int lemon (double p[2], double u) {
  //x^2 - (1-y^2)^3 = 0
  p[0] = pow(cos(u), 3);
  p[1] = sin(u);
}

int drawNormal(int (*func)(double[2], double), double u, double padd[2]) {
  //Find the tangent line at a point and then find the perpendicular line
  double newu;
  double xslope;
  double yslope;
  double p1[2], p2[2];
  double m;
  double h;
  newu = u + 0.0001;
  (*func)(p1, u);
  (*func)(p2, newu);
  xslope = p2[0] - p1[0];
  yslope = p2[1] - p1[1];
  if(yslope == 0) {
    padd[0] = p1[0];
    padd[1] = p1[1] + 10;
  } else if(xslope == 0) {
    padd[0] = p1[0] + 10;
    padd[1] = p1[1];
  } else {
    m = -1/(yslope/xslope);
    h = sqrt(1 + pow(m, 2));
    padd[0] = (1);
    padd[1] = (m);
  }
}   

void graph(int (*func)(double[2], double), double m[4][4], double minv[4][4],
	double start, double end, double increment) {
    double p[2];
    double p2[2];
    double u;
    int i;
    printf("Graphing. Start is %lf, end is %lf, increment is %lf.\n", start, end, increment);
    i = 0;
    for(u = start; u <= end; u+= increment) {
    	(*func)(p, u);
    	D3d_mat_mult_pt(p, m, p);
    	G_point(p[0], p[1]);
	if(i % 150 == 0) {
	  drawNormal(&(*func),u, p2);
	  //D3d_mat_mult_pt(p2, m, p2);
	  G_line(p[0], p[1], p[0] + p2[0], p[1]+ p2[1]);
	}
	i++;
    }
}
    
int main() {
    double m[4][4], minv[4][4];
    int num_movements;
    int mt_list[6];
    double p_list[6]; 

    G_init_graphics(600,600);
    G_rgb(0,0,0);
    G_clear();

	printf("Circle\n");
    //Making the Circle
    num_movements = 4;
    mt_list[0] = SX; p_list[0] = 50.0;
    mt_list[1] = SY; p_list[1] = 100.0;
    mt_list[2] = TX; p_list[2] = 300.0;
    mt_list[3] = TY; p_list[3] = 500.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&circle, m, minv, 0.25 * M_PI, 1.50 * M_PI, P_INCREMENT);

    G_wait_key();
    G_rgb(0,0,0);
    G_clear();

printf("Sum4\n");
    //Making the sum4
    num_movements = 4;
    mt_list[0] = SX; p_list[0] = 30.0;
    mt_list[1] = SY; p_list[1] = 60.0;
    mt_list[2] = TX; p_list[2] = 300.0;
    mt_list[3] = TY; p_list[3] = 300.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&sum4, m, minv, 0.50 * M_PI, 1.75 * M_PI, P_INCREMENT);

    G_wait_key();
    G_rgb(0,0,0);
    G_clear();

printf("Square\n");
    //Making the Square
    num_movements = 4;
    mt_list[0] = SX; p_list[0] = 150.0;
    mt_list[1] = SY; p_list[1] = 100.0;
    mt_list[2] = TX; p_list[2] = 500.0;
    mt_list[3] = TY; p_list[3] = 500.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&square, m, minv, 0, 2.0 * M_PI, P_INCREMENT);

    G_wait_key();
    G_rgb(0,0,0);
    G_clear();

printf("Asteroid\n");
    //Making the Asteroid
    num_movements = 5;
    mt_list[0] = SX; p_list[0] = 80.0;
    mt_list[1] = SY; p_list[1] = 40.0;
    mt_list[2] = RZ; p_list[2] = 45.0;
    mt_list[3] = TX; p_list[3] = 500.0;
    mt_list[4] = TY; p_list[4] = 300.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&asteroid, m, minv, 0, 2.0 * M_PI, P_INCREMENT);

    G_wait_key();
    G_rgb(0,0,0);
    G_clear();

printf("Hyperbola\n");
    //Making the Hyperbola
    num_movements = 5;
    mt_list[0] = NY; p_list[0] = 0.0;
    mt_list[1] = SX; p_list[1] = 100.0;
    mt_list[2] = SY; p_list[2] = 100.0;
    mt_list[3] = TX; p_list[3] = 250.0;
    mt_list[4] = TY; p_list[4] = 250.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&hyperbola, m, minv, -1, 2, P_INCREMENT);

    G_wait_key();
    G_rgb(0,0,0);
    G_clear();

printf("Parabola\n");
    //Making the Parabola
    num_movements = 5;
    mt_list[0] = SX; p_list[0] = 150.0;
    mt_list[1] = SY; p_list[1] = 50.0;
    mt_list[2] = RZ; p_list[2] = 60.0;
    mt_list[3] = TX; p_list[3] = 250.0;
    mt_list[4] = TY; p_list[4] = 250.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&parabola, m, minv, -1, 2, P_INCREMENT);

    G_wait_key();
    G_rgb(0,0,0);
    G_clear();

printf("Lemon\n");
    //Making of Lemon
    num_movements = 5;
    mt_list[0] = SX; p_list[0] = 150.0;
    mt_list[1] = SY; p_list[1] = 150.0;
    mt_list[2] = RZ; p_list[2] = 60.0;
    mt_list[3] = TX; p_list[3] = 600.0;
    mt_list[4] = TY; p_list[4] = 150.0;
    D3d_make_movement_sequence_matrix(m, minv, num_movements, mt_list, p_list);
    G_rgb(1,0,0);
    graph(&lemon, m, minv, 0, 2.0 * M_PI, P_INCREMENT);

    G_wait_key();
}
