//
// Created by Raffaele Dragone
//

#ifndef NBODY_UTIL_BODY_H
#define NBODY_UTIL_BODY_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define SOFTENING 1e-9f
typedef struct { int x, y, z; } PBody;//Position bodies.
typedef struct {int vx,vy,vz; } VBody;//Velocity bodies.

void randomizeBodies(PBody *p,VBody *vp,int n){//randomizza posizione e velocità di ogni particella.
    for(int i=0; i<n; ++i){
        p[i].x = rand() % (5000 ) + 1;
        p[i].y = rand() % (5000 ) + 1;
        p[i].z = rand() % (5000 ) + 1;
        vp[i].vx = rand() % (5000 ) + 1;
        vp[i].vy = rand() % (5000 ) + 1;
        vp[i].vz = rand() % (5000 ) + 1;
    }
}

void bodyForce(PBody *plocal_bodies,VBody *vlocal_bodies,int size_local, int dt,PBody *other_bodies, int size_other) { //Body Force . Se Other_bodies=NULL allora la computazione è relativa soltanto alle particelle locali al processo.
    int size_it = other_bodies!=NULL && size_other>0 ? size_other : size_local; // Dimensione della taglia di particelle sulla quale effettuare la computazione
    for(int i=0;i<size_local;++i){ //Per ogni particella locale
        int Fx = 0; int Fy = 0; int Fz = 0; //Init Force
        for(int j=0; j<size_it;++j){
            int dx = other_bodies!=NULL && size_other>0 ? (other_bodies[j].x-plocal_bodies[i].x) : (plocal_bodies[j].x - plocal_bodies[i].x);
            int dy = other_bodies!=NULL && size_other>0 ? (other_bodies[j].y-plocal_bodies[i].y) : (plocal_bodies[j].y - plocal_bodies[i].y);
            int dz = other_bodies!=NULL && size_other>0 ? (other_bodies[j].z-plocal_bodies[i].z) : (plocal_bodies[j].z - plocal_bodies[i].z);
            //int distSqr = dx + dy + dz ;
            //float invDist = (  1.0f /  sqrtf(distSqr));
            //float invDist3 = invDist * invDist * invDist;
            Fx +=  dx; Fy +=  dy ; Fz +=  dz ;
        }
        vlocal_bodies[i].vx += (dt*Fx); vlocal_bodies[i].vy += (dt*Fy); vlocal_bodies[i].vz += (dt*Fz);
    }
}

void split_qta(int nBodies,int tasks,int bodies_gsize[],int myRank){ //Calcola la taglia di particelle da assegnate ad ogni processo.
    //bodies_gsize è un array globale, contenente la taglia di ogni processo, visibile da tutti i processi. Questo riduce comunicazione ridondante durante la fase di computazione perchè non viene scambiata ogni volta la taglia.
    if(myRank==0){
        for( int rank = 0; rank<tasks; ++rank){
            int resto = (nBodies) % (tasks);
            bodies_gsize[rank] = (resto>0 && rank<resto ) ?   nBodies/tasks+1  : nBodies / tasks ;
        }
    }
    //send size from rank0 to all
    MPI_Bcast(bodies_gsize,tasks,MPI_INT,0,MPI_COMM_WORLD);
}

void init_bodies(int tasks, int myRank,int global_size, MPI_Datatype bodyes_datatp,int bodies_gsize[],PBody bodies[],VBody vbodies[]){
    int displs[tasks];
    PBody *b_app=NULL;
    VBody *vb_app=NULL;
    if(myRank==0){ //Master
        b_app = (PBody*) malloc(sizeof(PBody)* (global_size));
        vb_app = (VBody*) malloc(sizeof(VBody)* (global_size));
        //init bodies random.
        randomizeBodies(b_app,vb_app,(global_size));

        displs[0]=0;
        for( int rank = 1; rank<tasks; ++rank){
            displs[rank]=bodies_gsize[rank-1]+displs[rank-1];
        }
    }
    MPI_Scatterv(&b_app[0], bodies_gsize, displs, bodyes_datatp, &bodies[0] , bodies_gsize[myRank] , bodyes_datatp, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&vb_app[0], bodies_gsize, displs, bodyes_datatp, &vbodies[0] , bodies_gsize[myRank] , bodyes_datatp, 0, MPI_COMM_WORLD);
    if(myRank==0) {
        for(int i=0; i<bodies_gsize[myRank];++i){
            //Il rank0 mantiene soltanto le proprie particelle
            bodies[i]=b_app[i];
            vbodies[i]=vb_app[i];
        }
    }
    free(b_app);
    free(vb_app);
}

#endif //NBODY_UTIL_BODY_H
