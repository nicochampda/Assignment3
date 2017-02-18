#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
/*
#include "graphics/graphics.h"
*/
#include "file_operations/file_operations.h"

//Definition of Epsilon0
const double E0 = 0.001;

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
    int *part_index;
    double center_x;
    double center_y;
    double center_mass;
    double box_size;
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
    }else{
        return 1 ;
    }
}



//insert an elment at the end of an array
void insert_in_list(int *list,int new_element, int N){
    int* maliste = list;
    maliste = realloc(maliste,(N+1)*sizeof(int));
    maliste[N] = new_element;
}      
   
//Function that make recursively the quad tree
void makeQuadtree(quadtree *src, double xmin, double xmax, double ymin, double ymax, particule *particules){
    if (src->part_nbr > 1){  //general case
  printf("recurse %i \n",src->part_nbr);
        
        //calculate new space
        double xmid = (xmax - xmin)/2;
        double ymid = (ymax - ymin)/2;

        //Separate particules into the 4 quad
        int i;
        int N = src->part_nbr;
        
        int ul_len = 0;
        int *ul_list = (int*)malloc(sizeof(int));
        int ur_len = 0;
        int *ur_list = (int*)malloc(sizeof(int));
        int dl_len = 0;
        int *dl_list = (int*)malloc(sizeof(int));
        int dr_len = 0;
        int *dr_list = (int*)malloc(sizeof(int));

        for(i = 0; i<N; i++){
            int cur_part = src->part_index[i];
            double pos_x = particules[cur_part]->pos_x;
            double pos_y = particules[cur_part]->pos_y;
            if (is_in_rectangle(xmin, xmid, ymin, ymid, pos_x, pos_y)){
                insert_in_list(ul_list, i, ul_len);
                ul_len++;
            }else{
                if (is_in_rectangle(xmid, xmax, ymin, ymid, pos_x, pos_y)){
                    insert_in_list(ur_list, i, ur_len);
                    ur_len++;
                }else{
                    if (is_in_rectangle(xmin, xmid, ymid, ymax, pos_x, pos_y)){
                        insert_in_list(dl_list, i, dl_len);
                        dl_len++;
                    }else{
                        if (is_in_rectangle(xmid, xmax, ymid, ymax, pos_x, pos_y)){
                            insert_in_list(dr_list, i, dr_len);
                            dr_len++;
                        }
                    }
                }
            }
        }

        //Allocate memory for the 4 new quadtree
        quadtree *ul = (quadtree *)malloc(sizeof(quadtree));
        quadtree *ur = (quadtree *)malloc(sizeof(quadtree));
        quadtree *dl = (quadtree *)malloc(sizeof(quadtree));
        quadtree *dr = (quadtree *)malloc(sizeof(quadtree));

        //set ul list of index and nbr
        ul->part_nbr = ul_len;
        ul->part_index = ul_list;
        ur->part_nbr = ur_len;
        ur->part_index = ur_list;
        dl->part_nbr = dl_len;
        dl->part_index = dl_list;
        dr->part_nbr = dr_len;
        dr->part_index = dr_list;

        //Apply recursively on each subsquare
        makeQuadtree(ul ,xmin, xmid, ymin, ymid, particules);
        makeQuadtree(ur ,xmid, xmax, ymin, ymid, particules);
        makeQuadtree(dl ,xmin, xmid, ymid, ymax, particules);
        makeQuadtree(dr ,xmid, xmax, ymid, ymax, particules);
        //makeQuadtree(dr ,0, 0, 0, 0, particules);

        src->center_mass = ul->center_mass + ur->center_mass + dl->center_mass + dr->center_mass; 
        src->center_x = (ul->center_mass * ul->center_x + 
                         ur->center_mass * ur->center_x +
                         dl->center_mass * dl->center_x +
                         dr->center_mass * dr->center_x) / src->center_mass; 
        src->center_y = (ul->center_mass * ul->center_y + 
                         ur->center_mass * ur->center_y +
                         dl->center_mass * dl->center_y +
                         dr->center_mass * dr->center_y) / src->center_mass; 
        src->box_size = xmax -xmin;

        src->ul = ul;
        src->ur = ur;
        src->dl = dl;
        src->dr = dr;

    }else{
        if (src->part_nbr == 1){
            particule current = particules[src->part_index[0]];
            src->center_x = current->pos_x;
            src->center_y = current->pos_y;
            src->center_mass = current->mass;
            src->box_size = xmax -xmin;
            src->ul = NULL;
            src->ur = NULL;
            src->dl = NULL;
            src->dr = NULL;
        }else{
            src = NULL;
            /*src->center_x = 0;
            src->center_y = 0;
            src->center_mass = 0;
            src->box_size = xmax -xmin;
            src->ul = NULL;
            src->ur = NULL;
            src->dl = NULL;
            src->dr = NULL;*/
        }
    }
}

//Compute force recursively
void computeForce(double x, double y, double mass, double theta, quadtree *src, double Fx, double Fy, particule *particules){
    if (src->part_nbr > 1){
        double distancex = x - src->center_x;
        double distancey = y - src->center_y;
        double dist_box = sqrt(distancex*distancex + distancey*distancey); 
        if(src->box_size / dist_box > theta){//compare theta to the threshold
            //Apply recursion on all sub trees
            computeForce(x, y, mass, theta, src->ul, Fx, Fy, particules);
            computeForce(x, y, mass, theta, src->ur, Fx, Fy, particules);
            computeForce(x, y, mass, theta, src->dl, Fx, Fy, particules);
            computeForce(x, y, mass, theta, src->dr, Fx, Fy, particules);
        }else{
            //or calculate force from G on the particule
            double distancex = x - src->center_x;
            double distancey = y - src->center_y;
            double rij = sqrt(distancex*distancex + distancey*distancey);
            double cst_j = src->center_mass * 1.0 /((rij+E0)*(rij+E0)*(rij+E0));
            double cord_x = cst_j * distancex;
            double cord_y = cst_j * distancey;
            Fx += cord_x;
            Fy += cord_y;
        }

    }else{
        if (src->part_nbr ==1){
            //calculate force apply by the only particule of the quadtree
            double distancex = x - src->center_x;
            double distancey = y - src->center_y;
            double rij = sqrt(distancex*distancex + distancey*distancey);
            double cst_j = src->center_mass * 1.0 /((rij+E0)*(rij+E0)*(rij+E0));
            double cord_x = cst_j * distancex;
            double cord_y = cst_j * distancey;
            Fx += cord_x;
            Fy += cord_y;
        }
    }
}

int main (int argc, char *argv[]){


    if (argc != 7){
        printf("Wrong number of arguments given. Write:%s nbr_of_star filename nsteps delta_t theta graphics_0/1\n", argv[0]);
        return -1;
    }

    const int N = atoi(argv[1]);
    const char* fileName = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const double theta = atof(argv[5]);
    const int graphics = atoi(argv[6]);

    const double n = atof(argv[1]);
    particule particules[N];
  
  //Declaration of variables for the sympletic Euler integration method
 
    int i, p;
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
        printf("particules %i: x,y = %f,%f\n",i,particules[i]->pos_x, particules[i]->pos_y);
    }


    printf("reading files took %7.3f wall seconds.\n", get_wall_seconds()-time1);
    time1 = get_wall_seconds();

    //Euler sympletic integration method
    const double Gdelta_t = (-100.0/n)*delta_t;
    for (p=0; p<nsteps; p++) {
        //make the quadtree
        quadtree *root = malloc(sizeof(quadtree));
        root->part_nbr = N;
        int index[N];
        for (i = 0; i<N; i++)
            index[i] = i;
        root->part_index = index;
        printf("run make quadtree\n");
        makeQuadtree(root, 0, 1, 0, 1, particules);
        printf("quadtree %d %d %d %d %d \n",root->part_nbr, root->ul->part_nbr, root->ur->part_nbr, root->dl->part_nbr, root->dr->part_nbr);
	    //compute forces for each particule recursively 
        for (i = 0; i<N; i++){
            sum_Fx[i] = 0;
            sum_Fy[i] = 0;
            computeForce(particules[i]->pos_x, particules[i]->pos_y, particules[i]->mass, theta, root, sum_Fx[i], sum_Fy[i], particules);
        }
	    
	    //Update of the position and velocity of each particule
        for (i=0; i<N; i++){
            particules[i]->vel_x += Gdelta_t * sum_Fx[i];
            particules[i]->vel_y += Gdelta_t * sum_Fy[i];
            particules[i]->pos_x += delta_t*particules[i]->vel_x;
            particules[i]->pos_y += delta_t*particules[i]->vel_y;
        }

        //Free quad tree


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
    
    /*for (i=0; i<N; i++){
        free(particules[i]);
    }*/

return 0;
}
