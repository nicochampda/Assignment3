#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "graphics/graphics.h"
#include "file_operations/file_operations.h"

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

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
    const int N = atoi(argv[1]);
    const char* fileName = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const int graphics = atoi(argv[5]);
  //Definition of Epsilon0
    const double E0 = 0.001;
    particule particules[N];
 
    int i, j, p;
    double rij, cst_j, cord_x, cord_y;
    double distancex, distancey;
    double buf[5*N];
    double output[5*N];
    double sum_Fx, sum_Fy, Ax, Ay;

    double time1;
    double positions_x[nsteps][N], positions_y[nsteps][N];

    const float circleRadius = 0.05/N, circleColor = 0;
    const int windowWidth = 800;
    float L=1, W=1;
    double x, y;

    if (graphics == 1){
        InitializeGraphics(argv[0], windowWidth, windowWidth);
        SetCAxes(0,1);
    }

    time1 = get_wall_seconds();
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
  printf("reading files took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    // Here we define sigma(mj/(rpj+Epsilon0)^3),j=0..N-1,j!=p)
  time1 = get_wall_seconds();
    double G = -100/N;
    for (p=0; p<nsteps; p++) {
        if (graphics == 1) ClearScreen();
        for (i=0; i<N; i++) {
            sum_Fx = 0;
            sum_Fy = 0;
            for (j=0; j<N; j++) {
                distancex = particules[i]->pos_x - particules[j]->pos_x;
                distancey = particules[i]->pos_y - particules[j]->pos_y;
                rij = sqrt(distancex*distancex + distancey*distancey);
                cst_j = (particules[j]->mass) * (1.0 / ((rij+E0)*(rij+E0)*(rij+E0)));
                cord_x = cst_j * (distancex);
                cord_y = cst_j * (distancey);
                sum_Fx += cord_x;
                sum_Fy += cord_y;
            }

            Ax = G * sum_Fx;
            Ay = G * sum_Fy;
            particules[i]->vel_x += delta_t*Ax;
            particules[i]->vel_y += delta_t*Ay;
            particules[i]->pos_x += delta_t*particules[i]->vel_x;
            particules[i]->pos_y += delta_t*particules[i]->vel_y;

            if (graphics == 1){
                x = particules[i]->pos_x;
                y = particules[i]->pos_y;
                DrawCircle(x, y, L, W, circleRadius, circleColor);
            }

            if (graphics == 1){
                Refresh();
                usleep(2000);
            }
        }
    }

    if (graphics == 1){
        FlushDisplay();
        CloseDisplay();
    }

  printf("calculations took %7.3f wall seconds.\n", get_wall_seconds()-time1);
  time1 = get_wall_seconds();
    for (i=0; i<N; i++) {
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
  printf("writing took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    
    for (i=0; i<N; i++){
        free(particules[i]);
    }
   /* if (graphics == 1){
        const float circleRadius = 0.05/N, circleColor = 0;
        const int windowWidth = 800;
        float L=1, W=1;
	double x, y;

        InitializeGraphics(argv[0], windowWidth, windowWidth);
        SetCAxes(0,1);
        for (p=0; p<nsteps; p++){
            ClearScreen();
            for (i=0; i<N; i++){
                x = positions_x[p][i];
                y = positions_y[p][i];
                //  printf("p %i x %lf y %lf\n", p, x, y);
                DrawCircle(x, y, L, W, circleRadius, circleColor);
            }
            Refresh();
            usleep(2000);
        }
        FlushDisplay();
        CloseDisplay();
    }*/ 
return 0;
}
