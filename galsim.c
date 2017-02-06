#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct particules {
    double pos_x;
    double pos_y;
    double mass;
    double vel_x;
    double vel_y;
}particule;


int main (int argc, char *argv[]){

    if (argc != 6){
        printf("Use %s nbr_of_star filename nsteps delta_t graphics_0/1", argv[0]);
        exit(1);
    }
    int N = atoi(argv[1]);
    int nsteps = atoi(argv[3]);
    int delta_t = atoi(argv[4]);
    int graphics = atoi(argv[5]);
    double d;
    FILE *file = fopen(argv[2],"r");
    fscanf(file, "%lf", &d);
    printf("file open \n");
    printf("%lf\n", d);
    fclose(file);

    int A;
    int sum;
    int F;
    double X[10], Y[10], UX[10], UY[10];

//vu que je comprends pas trop encore comment extraire les données du fichier binaire j'ai posé :
// UX -> liste velocité selon x
//UY -> liste vélocité selon y
//X-> liste position x
//Y -> liste position y

    for (int m=0; m<nsteps; m++){
        for (int i=0; i<N; i++){
        sum = 0;
            for(int j=0; j<N-1; j++){
                if (j != i){
                sum += m/pow((pow(X[i]-X[j],2)-pow(Y[i]-Y[j],2)),2);
                }
            }  
            F = -100/N*m*sum;
            A = F/m;
            X[i] = X[i] + delta_t*UX[i] + pow(delta_t,2)*A;
            Y[i] = Y[i] + delta_t*UY[i] + pow(delta_t,2)*A;
            UX[i] = UX[i] + delta_t*A;
            UY[i] = UY[i] + delta_t*A;
        }
    }
 
return 0;
}
