// Created by Raffaele Dragone
#include <stdio.h>
#include "utilBody.h"
#include <mpi.h>
void write_file(Body* local_bodies,int size);
int main(const int argc, const char** argv) {

    double start,end;
    int nBodies = atoi(argv[1]);
    const int nIters = atoi(argv[2]);  // simulation iterations
    const int print_bodies = argv[3]!=NULL ? atoi(argv[3]) : 0; //1 -> print | 0 -> no
    const int print_iterations = argv[4]!=NULL ? atoi(argv[4]) : 0; //1 -> print | 0 -> no
    const float dt = 0.01f; // time step

    MPI_Init(NULL,NULL);
    start = MPI_Wtime();
    Body *p = (Body*) malloc(sizeof(Body)*nBodies);
    randomizeBodies(p,nBodies);

    for (int iter = 0; iter < nIters; iter++) {
    	if(print_iterations==1){ //print computation state
            printf("Iteration : %d \n",iter);
            fflush(stdout);
        }
        bodyForce(p, dt, nBodies,iter); // compute interbody forces
        //Change positions
        for (int i = 0 ; i < nBodies; i++) { // integrate position
            p[i].x += p[i].vx*dt;
            p[i].y += p[i].vy*dt;
            p[i].z += p[i].vz*dt;
        }
    }
    end=MPI_Wtime();
    printf(" \n Time in ms = %f \n",(end-start));
    if(print_bodies==1)
        write_file(p,nBodies);
    MPI_Finalize();
    free(p);
}
void write_file(Body* bodies,int size) {
    FILE * fp;
    // open the file for writing
    fp = fopen ("output.csv","w");
    for(int i=0; i<size;++i){
        fprintf (fp, "p[%d]\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\n",i,bodies[i].x,bodies[i].y,bodies[i].z,bodies[i].vx,bodies[i].vy,bodies[i].vz);
    }
    // close the file
    fclose (fp);

}