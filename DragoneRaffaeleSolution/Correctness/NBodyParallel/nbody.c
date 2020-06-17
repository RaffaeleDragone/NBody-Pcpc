// Created by Raffaele Dragone
#include <stdio.h>
#include "util_body.h"
#include <mpi.h>

void run_step(PBody *plocal_bodies,VBody *vlocal_bodies,PBody *other_bodies,int bodies_gsize[],MPI_Datatype basic_bodiesdp,int tasks,int myRank);
void collect_bodies(PBody *plocal_bodies,VBody *vlocal_bodies,int bodies_gsize[],MPI_Datatype bodies_datatp,int tasks,int myRank,int size_global);
void write_bodies_csv(PBody *bodies,VBody *vbodies,int size);
int main(const int argc, const char** argv) {

    const int num_bodies = argv[1]!=NULL ? atoi(argv[1]) : 0; //Numero di particelle
    const int nIters = argv[2]!=NULL ? atoi(argv[2]) : 0;  // Numero di iterazioni
    const int print_bodies = argv[3]!=NULL ? atoi(argv[3]) : 0; //1 -> scrivi su csv posizione e velocità di ogni particella | 0 -> no
    const int print_iterations = argv[4]!=NULL ? atoi(argv[4]) : 0; //1 -> scrivi lo stato delle iterazioni su stdout | 0 -> no

    //init MPI
    int myRank,tasks;
    MPI_Status status;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&tasks);
    int bodies_gsize[tasks];//array with body size of every rank.
    double start,end;



    MPI_Datatype bodyes_datatp, old_types[1];
    int blockcounts[1];
    MPI_Aint offsets[1];
    // setup description of the  MPI_FLOAT fields velocityX,Y,Z | PositionX,Y,Z
    offsets[0] = 0;
    old_types[0] = MPI_INT;
    blockcounts[0] = 3;

    MPI_Type_create_struct(1, blockcounts, offsets, old_types, &bodyes_datatp);
    MPI_Type_commit(&bodyes_datatp);
    //Commit data type utilizzato per la comunicazione. Il tipo di dato derivato conterrà 3 campi float.



    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    split_qta(num_bodies,tasks,bodies_gsize,myRank); //split qta
    PBody *plocal_bodies = (PBody*) malloc(sizeof(PBody)*bodies_gsize[myRank]);
    VBody *vlocal_bodies = (VBody*) malloc(sizeof(VBody)*bodies_gsize[myRank]);

    init_bodies(tasks,myRank,num_bodies,bodyes_datatp,bodies_gsize,plocal_bodies,vlocal_bodies); //init bodies

    int size_other = (num_bodies - bodies_gsize[myRank] ); //size total for other bodies
    PBody *other_bodies = (PBody*) malloc(sizeof(PBody) * size_other);//Other bodies



    for(int it=0; it<nIters; ++it){ //for each iteration
        if(myRank==0 && print_iterations==1){ //print computation state
            printf("Iteration : %d \n",it);
            fflush(stdout);
        }
        run_step(plocal_bodies,vlocal_bodies,other_bodies,bodies_gsize,bodyes_datatp,tasks,myRank); //run step and refresh position of local bodies
    }



    if(myRank==0){
        plocal_bodies = (PBody*) realloc(plocal_bodies,sizeof(PBody) * num_bodies); //size global
        vlocal_bodies = (VBody*) realloc(vlocal_bodies,sizeof(VBody) * num_bodies); //size global
    }
    //Collect all the bodies
    collect_bodies(plocal_bodies,vlocal_bodies,bodies_gsize,bodyes_datatp,tasks,myRank,num_bodies);
    MPI_Type_free(&bodyes_datatp);



    MPI_Barrier(MPI_COMM_WORLD);
    end=MPI_Wtime();
    if(myRank==0){
        printf(" \n Time in s = %f \n",(end-start));
    }
    if(print_bodies==1 && myRank==0)
        write_bodies_csv(plocal_bodies,vlocal_bodies,num_bodies); //write on csv

    free(vlocal_bodies);
    free(plocal_bodies);

    MPI_Finalize();
}


void run_step(PBody *local_bodies,VBody *vlocal_bodies,PBody *other_bodies,int bodies_gsize[],MPI_Datatype basic_bodiesdp,int tasks,int myRank){
    MPI_Request requests[tasks];
    int dt = 1; // time step
    int idx=0;
    for(int root=0; root<tasks;++root){
        if(root==myRank){
            //Send Asyncronus Broadcast
            MPI_Ibcast(local_bodies,bodies_gsize[myRank],basic_bodiesdp,myRank,MPI_COMM_WORLD,&requests[myRank]);
            //Invia le proprie particelle in modo asincrono a tutti gli altri processori tramite una IBcast
        }else{
            int pos = root==0 ? 0 :  (idx*(bodies_gsize[root-1]));
            //Receive Asyncronus Broadcast
            MPI_Ibcast(&other_bodies[pos],(bodies_gsize[root]) ,basic_bodiesdp,root,MPI_COMM_WORLD,&requests[root]);
            //Riceve in modo asincrono le posizioni delle particelle di tutti gli altri processori tramite delle IBcast
            ++idx;
        }
    }
    bodyForce(local_bodies,vlocal_bodies,bodies_gsize[myRank],dt,NULL,0);//Inizia a computare sulle proprie particelle
    int req_left=tasks-1;
    int idxReq;
    MPI_Status stats;
    while(req_left>0){//Finché non sono terminate tutte le ricezioni
        MPI_Waitany(tasks,requests,&idxReq,&stats); //Wait any receive
        if(idxReq!=myRank && idxReq<tasks){
            int nElem = bodies_gsize[idxReq]; //num di particelle ricevute dallo specifico processo
            idxReq+= (idxReq>myRank) ? -1 : 0;
            int start = idxReq*nElem; //posizione iniziale delle particelle ricevute.
            bodyForce(local_bodies,vlocal_bodies,bodies_gsize[myRank],dt,&other_bodies[start],nElem);
            req_left--;
        }
    }
    for (int i = 0 ; i < bodies_gsize[myRank]; i++) { // integrate position
        local_bodies[i].x += vlocal_bodies[i].vx*dt;
        local_bodies[i].y += vlocal_bodies[i].vy*dt;
        local_bodies[i].z += vlocal_bodies[i].vz*dt;

    }
}

void collect_bodies(PBody *local_bodies,VBody *vlocal_bodies,int bodies_gsize[],MPI_Datatype bodyes_datatp,int tasks,int myRank,int size_global){
    int displs_recv[tasks];
    if(myRank==0){
        displs_recv[0]=0;
        for( int rank = 1; rank<tasks; ++rank){
            displs_recv[rank]=bodies_gsize[rank-1]+displs_recv[rank-1];
        }
    }
    MPI_Gatherv(local_bodies,bodies_gsize[myRank],bodyes_datatp,local_bodies,bodies_gsize,displs_recv,bodyes_datatp,0,MPI_COMM_WORLD);
    //Collect delle positions
    MPI_Gatherv(vlocal_bodies,bodies_gsize[myRank],bodyes_datatp,vlocal_bodies,bodies_gsize,displs_recv,bodyes_datatp,0,MPI_COMM_WORLD);
    //Collect delle velocities
}

void write_bodies_csv(PBody *bodies,VBody *vbodies,int size){
    FILE * fp;
    // open the file for writing
    fp = fopen ("output.csv","w");
    for(int i=0; i<size;++i){
        fprintf (fp, "p[%d]\t%d\t%d\t%d\t%d\t%d\t%d\n",i,bodies[i].x,bodies[i].y,bodies[i].z,vbodies[i].vx,vbodies[i].vy,vbodies[i].vz);
    }
    // close the file
    fclose (fp);
}