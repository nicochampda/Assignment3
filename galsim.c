int N;
int nsteps;
int delta_t;

int A;
int sum;
int F;


//vu que je comprends pas trop encore comment extraire les données du fichier binaire j'ai posé :
// UX -> liste velocité selon x
//UY -> liste vélocité selon y
//X-> liste position x
//Y -> liste position y

for (int m=0,m<nsteps,m++){
	for (i=0,i<N,i++){
  sum==0;
  	for(j=0,j<N-1,j++){
    	if (j!=i){
       sum=sum+m[j]/pow((pow(X[i]-X[j],2)-pow(Y[i]-Y[j],2)),2);
      }
      else{
      sum=sum+0
      }
    }  
    F=-100/N*m[i]*sum;
    A=F/m[i];
    X[i]=X[i]+delta_t*UX[i]+pow(delta_t,2)*A;
    Y[i]=Y[i]+delta_t*UY[i]+pow(delta_t,2)*A;
		UX[i]=UX[i]+delta_t*A;
    UY[i]=UY[i]delta_t*A;
	}
}

