#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
/*
#include "graphics/graphics.h"
*/
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


//Definition of quadtree
typedef struct node{
    int part_index[];
    double center_x;
    double center_y;
    double center_mass;
    int part_nbr;

    struct node *ul;
    struct node *ur;
    struct node *dl;
    struct node *dr;
}quadtree;
// function to check if a particle p is in a rectangle of a quadtree
int is_in_rectangle(double xmin, double xmax, double ymin, double ymax, double pos_x, double pos_y){
       if(pos_x>xmax || pos_x<xmin || pos_y>ymax || pos_y<ymin){
      return 0 ;
}
       else{
      return 1 ;
}

void insert_in_list
      
    
}
//Function that make recursively the quad tree
quadtree makeQuadtree(quadtree *src, double xmin, double xmax, double ymin, double ymax, int *list, int len){
    if (src->part_nbr > 1){  //general case
        
        //calculate new space
        double xmid = (xmax - xmin)/2;
        double ymid = (ymin - ymax)/2;

        //Separate particules into the 4 quad

        //Allocate memory for the 4 new quadtree
        quadtree *ul = (quadtree *)malloc(sizeof(quadtree));
        quadtree *ur = (quadtree *)malloc(sizeof(quadtree));
        quadtree *dl = (quadtree *)malloc(sizeof(quadtree));
        quadtree *dr = (quadtree *)malloc(sizeof(quadtree));

        //Apply recursively on each subsquare
        quadtree->ul = makeQuadtree(ul ,xmin, xmid, ymin, ymid);
        quadtree->ur = makeQuadtree(ur ,xmid, xmax, ymin, ymid);
        quadtree->dl = makeQuadtree(dl ,xmin, xmid, ymid, ymax);
        quadtree->dr = makeQuadtree(dr ,xmid, xmax, ymid, ymax);

        //merge result
    }else{
        if (src->part_nbr == 1){
            particule current = particules[part_index[0]];
            src->center_x = particule->pos_x;
            src->center_y = particule->pos_y;
            src->center_mass = particule->mass
        }else{
            src->center_x = 0;
            src->center_y = 0;
            src->center_mass = 0;
        }
    }
}

//Compute force recursively
void computeForce(){

}

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
    double sum_Fx[N], sum_Fy[N];

  //Declaration of the inputs and outputs that will be respectively in filename and result.gal
    double input[5*N];
    double output[5*N];
    

  //Declaration of positions for the graphic part 

    /*const float circleRadius = 0.005, circleColor = 0;
    const int windowWidth = 800;
    const float L=1, W=1;
    double x, y;

    if (graphics == 1){
        InitializeGraphics(argv[0], windowWidth, windowWidth);
        SetCAxes(0,1);
    }*/

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

         //make the quadtree
	quadtree repartition = makeQuadtree();

	 //compute forces for each particule recursively
	 computeForce();


	    
	    
    	/*
	for (i=0; i<N; i++) {
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

    if (graphics == 1){
        FlushDisplay();
        CloseDisplay();
    }*/
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
