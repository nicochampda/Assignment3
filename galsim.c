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
// Definition of particles by their positions,mass and velocities
typedef struct particule {
    double pos_x;
    double pos_y;
    double mass;
    double vel_x;
    double vel_y;
}*particule;


int main (int argc, char *argv[]){


    if (argc != 6){
        printf("Wrong number of arguments given. Write:%s nbr_of_star filename nsteps delta_t graphics_0/1", argv[0]);
        return -1;
    }

    const int N = atoi(argv[1]);
    const char* fileName = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const int graphics = atoi(argv[5]);

    const double n = atof(argv[1]);
  //Definition of Epsilon0
    const double E0 = 0.001;
    particule particules[N];
  
  //Declaration of variables for the sympletic Euler integration method
 
    int i, j, p;
    double rij, cst_j;
    double distancex, distancey;
    double cord_x, cord_y;

  //Declaration of the inputs and outputs that will be respectively in filename and result.gal
    double input[5*N];
    double output[5*N];
    

  //Declaration of positions for the graphic part 

    const float circleRadius = 0.005, circleColor = 0;
    const int windowWidth = 800;
    const float L=1, W=1;
    double x, y;

    if (graphics == 1){
        InitializeGraphics(argv[0], windowWidth, windowWidth);
        SetCAxes(0,1);
    }

   //for time measures of the program 
    double time1;

    time1 = get_wall_seconds();

    if (read_doubles_from_file(5*N, input, fileName) != 0){
        printf("Error reading file \n");
        return -1;
  
    }

    // Initializing particules data with the input file
    for(i = 0; i<N; i++){
        particules[i] = (struct particule *)malloc(sizeof(struct particule));
    	particules[i]->pos_x = input[i*5 + 0];
    	particules[i]->pos_y = input[i*5 + 1];     
    	particules[i]->mass  = input[i*5 + 2];     
    	particules[i]->vel_x = input[i*5 + 3];     
    	particules[i]->vel_y = input[i*5 + 4];     
    }


    printf("reading files took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    time1 = get_wall_seconds();

    //Euler sympletic integration method
    const double Gdelta_t = (-100.0/n)*delta_t;
    for (p=0; p<nsteps; p++) {
        double sum_Fx[N], sum_Fy[N];
        if (graphics == 1) ClearScreen();
        for (i=0; i<N; i++) {
            for (j=i; j<N; j++) {
                distancex = particules[i]->pos_x - particules[j]->pos_x;
                distancey = particules[i]->pos_y - particules[j]->pos_y;
                rij = sqrt(distancex*distancex + distancey*distancey);
                cst_j = 1.0 /((rij+E0)*(rij+E0)*(rij+E0));
                cord_x = cst_j * distancex;
                cord_y = cst_j * distancey;
                sum_Fx[i] += particules[j]->mass * cord_x; 
                sum_Fy[i] += particules[j]->mass * cord_y;
                sum_Fx[j] -= particules[i]->mass * cord_x;
                sum_Fx[j] -= particules[i]->mass * cord_y;
            } 
        }

        for (i=0; i<N; i++){
            particules[i]->vel_x += Gdelta_t * sum_Fx[i];
            particules[i]->vel_y += Gdelta_t * sum_Fy[i];
            particules[i]->pos_x += delta_t*particules[i]->vel_x;
            particules[i]->pos_y += delta_t*particules[i]->vel_y;

            if (graphics == 1){
                x = particules[i]->pos_x;
                y = particules[i]->pos_y;
                DrawCircle(x, y, L, W, circleRadius, circleColor);
            }
        }
        if (graphics == 1){
            Refresh();
            usleep(2000);
        }
    }

    if (graphics == 1){
        FlushDisplay();
        CloseDisplay();
    }

    printf("calculations took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    time1 = get_wall_seconds();


 //We put in output the result of Euler's method and put output in the file result.gal
    for (i=0; i<N; i++) {
        output[i*5 + 0] = particules[i]->pos_x;
        output[i*5 + 1] = particules[i]->pos_y;
        output[i*5 + 2] = particules[i]->mass;
        output[i*5 + 3] = particules[i]->vel_x;
        output[i*5 + 4] = particules[i]->vel_y;
    }

    if (write_doubles_to_file(5*N, output, "result.gal") != 0){
        printf("Error writing file");
        return -1;
    }


  printf("writing took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    
    for (i=0; i<N; i++){
        free(particules[i]);
    }

return 0;
}
