#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "graphics/graphics.h"
#include "file_operations/file_operations.h"

const float circleRadius = 0.025, circleColor = 0;
const int windowWidth = 800;

typedef struct particule {
    double pos_x;
    double pos_y;
    double mass;
    double vel_x;
    double vel_y;
}particule;


int main (int argc, char *argv[]){

    if (argc != 6){
        printf("Use %s nbr_of_star filename nsteps delta_t graphics_0/1", argv[0]);
        return -1;
    }
    int N = atoi(argv[1]);
    const char* fileName = argv[2];
    int nsteps = atoi(argv[3]);
    int delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    double d;

    printf("file open \n");

    double buf[5*N];
    if (read_doubles_from_file(5*N, buf, fileName) != 0){
        printf("Error reading file \n");
        return -1;
    }
    int i;
    for(i = 0; i<N; i++){
    	particules[i].pos_x = buf[i*5 + 0];
    	particules[i].pos_y = buf[i*5 + 1];     
    	particules[i].mass  = buf[i*5 + 2];     
    	particules[i].vel_x = buf[i*5 + 3];     
    	particules[i].vel_y = buf[i*5 + 4];     
    }
    for (i=0; i<nsteps; i++) {
    for (p=0;p<N;p++) {
      sum=0;
      for (j=0;j<N-1;j++) {
        if (j != p) {
          distancex = particules[p].pos_x - particules[j].pos_x;
          distancey = particules[p].pos_y - particules[j].pos_y;
          sum += particules[j].mass / pow(  pow(distancex,2) - pow(distancey, 2) ,2);
        }
      }
      F = (-100/N) * particules[p].mass * sum;
      A = F/particules[p].mass;

      particules[p].pos_x += delta_t*particules[p].vel_x + pow(delta_t,2)*A;
      particules[p].pos_y += delta_t*particules[p].vel_y + pow(delta_t,2)*A;
      particules[p].vel_x += delta_t*A;
      particules[p].vel_y += delta_t*A;
    }
  }
  for (i=0;i<N;i++) {
    printf("---------------------------------\n");
    printf("x:%f\ny:%f\nmass:%f\nvel_x:%f\nvel_y:%f\n", particules[i].pos_x, particules[i].pos_y, particules[i].mass, particules[i].vel_x, particules[i].vel_y);
  }
  return 0;
}


    if (graphics == 1){
        float L=1, W=1;
        int i, j;
	double x, y;

        InitializeGraphics(argv[0], windowWidth, windowWidth);
        SetCAxes(0,1);
        for (i=0; i<nsteps; i++){
            ClearScreen();
            for (j=0; j<N; j++){
                /*set each particule*/
                x = 0;
                j = 1;
                DrawCircle(x, y, L, W, circleRadius, circleColor);
            }
            Refresh();
            usleep(delta_t*1000);
        }
        FlushDisplay();
        CloseDisplay();
    }
return 0;
}
