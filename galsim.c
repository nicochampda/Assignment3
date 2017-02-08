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
}*particule;


int main (int argc, char *argv[]){

    if (argc != 6){
        printf("Use %s nbr_of_star filename nsteps delta_t graphics_0/1", argv[0]);
        return -1;
    }
    int N = atoi(argv[1]);
    const char* fileName = argv[2];
//    const char* result.gal;
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
  //Definition of Epsilon0
    double E0 = 0.001;
    particule particules[N];
 
    int i, j, p;
    double A, F;
    double rij, cst_j, mj, sum_dist_square, cord_x, cord_y;
    double distancex, distancey;
    double buf[5*N];
    double output[5*N];
    

    if (read_doubles_from_file(5*N, buf, fileName) != 0){
        printf("Error reading file \n");
        return -1;
    }
    for(i = 0; i<N; i++){
        particules[i] = (struct particule *)malloc(sizeof(struct particule));
    	particules[i]->pos_x = buf[i*5 + 0];
    	particules[i]->pos_y = buf[i*5 + 1];     
    	particules[i]->mass  = buf[i*5 + 2];     
    	particules[i]->vel_x = buf[i*5 + 3];     
    	particules[i]->vel_y = buf[i*5 + 4];     
    }
  printf("INIT\n");
  for (i=0; i<N; i++) {
    printf("---------------------------------\n");
    printf("x:%lf\ny:%lf\nmass:%lf\nvel_x:%lf\nvel_y:%lf\n", particules[i]->pos_x, particules[i]->pos_y, particules[i]->mass, particules[i]->vel_x, particules[i]->vel_y);
  }
  // Here we define sigma(mj/(rpj+Epsilon0)^3),j=0..N-1,j!=p)
    for (p=0; p<nsteps; p++) {
        for (i=0; i<N; i++) {
            sum_dist_square = 0;
            
            for (j=0; j<N; j++) {
                if (j != i) {
                distancex = particules[i]->pos_x - particules[j]->pos_x;
                distancey = particules[i]->pos_y - particules[j]->pos_y;
                rij = pow(pow(distancex,2) + pow(distancey, 2),0.5);
		mj = particules[j]->mass;
                cst_j = mj/(pow(rij+E0,3));
                cord_x = cst_j * (distancex);
                cord_y = cst_j * (distancey);
                sum_dist_square += pow(cord_x,2) + pow(cord_y,2);              

 //utilise la 2e formule pour la stabilite
                }
            }
            F = (-100/N) * (particules[i]->mass) * pow(sum_dist_square,0.5); 
            A = F/particules[i]->mass;
            particules[i]->vel_x += delta_t*A;
            particules[i]->vel_y += delta_t*A;
            particules[i]->pos_x += delta_t*particules[i]->vel_x;
            particules[i]->pos_y += delta_t*particules[i]->vel_y;
        }
  }
printf("\nFINAL\n");
  for (i=0; i<N; i++) {
    printf("---------------------------------\n");
    printf("x:%f\ny:%f\nmass:%f\nvel_x:%f\nvel_y:%f\n", particules[i]->pos_x, particules[i]->pos_y, particules[i]->mass, particules[i]->vel_x, particules[i]->vel_y);
    output[i*5 + 0] = particules[i]->pos_x;
    output[i*5 + 1] = particules[i]->pos_y;
    output[i*5 + 2] = particules[i]->mass;
    output[i*5 + 3] = particules[i]->vel_x;
    output[i*5 + 4] = particules[i]->vel_x;
  }

if (write_doubles_to_file(5*N, output, "result.gal") != 0){
    printf("Error writing file");
    return -1;
}
 // for (i=0;i<N;i++){
    
 //  }

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
                y = 1;
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
