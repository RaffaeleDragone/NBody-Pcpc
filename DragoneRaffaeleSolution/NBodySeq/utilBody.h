//
// Created by Raffaele Dragone
//

#ifndef NBODYSEQ_UTILBODY_H
#define NBODYSEQ_UTILBODY_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define SOFTENING 1e-9f
typedef struct { float x, y, z, vx, vy, vz; } Body;

void randomizeBodies(Body *p,int n){
    for(int i=0; i<n; ++i){
        p[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
}

void bodyForce(Body *p, float dt, int n,int iter) {
    for (int i =0; i < n; i++) {
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;
        for (int j = 0; j < n; j++) {
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr =  dx*dx + dy*dy + dz*dz + SOFTENING;
            float invDist =     (  1.0f /  sqrtf(distSqr));
            float invDist3 =  invDist * invDist * invDist;
            Fx +=  dx * invDist3; Fy +=  dy * invDist3; Fz +=  dz * invDist3;

        }
        p[i].vx += (dt*Fx); p[i].vy += (dt*Fy); p[i].vz += (dt*Fz);
    }

}

#endif //NBODYSEQ_UTILBODY_H
