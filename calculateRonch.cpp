#include <math.h>
#include <stdint.h>
#include <iostream>
#include <emscripten.h>

extern "C" {

float* findMinMax(float *image,int numPx) {
    float minmax[2] = {image[0], image[0]};
    float valCur;

    for (int j=0; j<numPx; j++) {
        for (int i=0; i<numPx; i++) {
            valCur = *(image+ i + (j*numPx));
            if (valCur < minmax[0]) {
                minmax[0] = valCur;
            }
            if (valCur > minmax[1]) {
                minmax[1] = valCur;
            }
        }
        std::cout << "\n" <<std::endl;
    }
    auto minmaxPtr = &minmax[0];
    return minmaxPtr;
}

int normalizeImage0255 (float* image, int numPx) {
    float* minmaxPtr = findMinMax(image, numPx);
    for (int j=0; j<numPx; j++) {
        for (int i=0; i<numPx; i++) {
            image[i + (j*numPx)] = image[i + (j*numPx)]-minmaxPtr[0];
            image[i + (j*numPx)] = image[i + (j*numPx)]/minmaxPtr[1]*255;
        }
    }
    return 0;
}

float* createMeshXX (int numPx) {
    float xx[numPx][numPx];
    for (int j=0; j<numPx; j++) {
        for (int i=0; i<numPx; i++) {
            xx[i][j] = i-(numPx/2+1);
        }
    }
    auto imagePtr = &xx[0][0];
    return imagePtr;
}
float* createMeshYY (int numPx) {
    float yy[numPx][numPx];
    for (int j=0; j<numPx; j++) {
        for (int i=0; i<numPx; i++) {
            yy[i][j] = j-(numPx/2+1);
        }
    }
    auto imagePtr = &yy[0][0];
    return imagePtr;
}

float* createMeshRR (int numPx, float al_max) {
    float rr[numPx][numPx];
    int center = numPx/2;
    for (int j=0; j < numPx; j++) {
        for (int i=0; i < numPx; i++) {
            rr[i][j] = sqrt( pow(i-center,2)+pow(j-center,2) )*al_max/center;
            }        
    }
    auto imagePtr = &rr[0][0];
    return imagePtr;
}

float* createMeshPP (int numPx) {
    float pp[numPx][numPx];
    int center = numPx/2;
    for (int j=0; j < numPx; j++) {
        for (int i=0; i < numPx; i++) {
            pp[i][j] = atan2( j-center,i-center );
            }        
    }
    auto imagePtr = &pp[0][0];
    return imagePtr;
}
float* createObjAp (float *imageRR, int numPx, float objApR) {
    float objAp[numPx][numPx];
    for (int j=0; j < numPx; j++) {
        for (int i=0; i < numPx; i++) {
            if (imageRR[ j+i*numPx] < objApR) {
                objAp[i][j] = 1;
            } else {
                objAp[i][j] = 0;
            }
        }
    }
    auto imagePtr = &objAp[0][0];
    return imagePtr;

}


float* calcRonch (int numPx, float al_max, float objApR) {
    float* rrPtr = createMeshRR(numPx, al_max);
    float* imagePtr = createObjAp(rrPtr,numPx, objApR);
    return imagePtr;
}   

}
