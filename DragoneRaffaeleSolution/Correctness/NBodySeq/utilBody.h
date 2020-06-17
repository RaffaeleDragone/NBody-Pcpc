//
// Created by Raffaele Dragone
//

#ifndef NBODYSEQ_UTILBODY_H
#define NBODYSEQ_UTILBODY_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define SOFTENING 1e-9f
typedef struct { int x, y, z, vx, vy, vz; } Body;

void randomizeBodies(Body *p,int n){
    for(int i=0; i<n; ++i){
        p[i].x = rand() % (5000 ) + 1;
        p[i].y = rand() % (5000 ) + 1;
        p[i].z = rand() % (5000 ) + 1;
        p[i].vx = rand() % (5000 ) + 1;
        p[i].vy = rand() % (5000 ) + 1;
        p[i].vz = rand() % (5000 ) + 1;
    }
}

void bodyForce(Body *p, int dt, int n,int iter) {
    for (int i =0; i < n; i++) {
        int Fx = 0; int Fy = 0; int Fz = 0;
        for (int j = 0; j < n; j++) {
            int dx = p[j].x - p[i].x;
            int dy = p[j].y - p[i].y;
            int dz = p[j].z - p[i].z;
            //int distSqr =  dx + dy + dz ;
            //int invDist =   (int) sqrt(distSqr);
            //double invDist3 =  invDist * invDist * invDist;
            Fx +=  dx ; Fy +=  dy ; Fz +=  dz ;
        }
        p[i].vx += (dt*Fx); p[i].vy += (dt*Fy); p[i].vz += (dt*Fz);
    }

}

#endif //NBODYSEQ_UTILBODY_H
