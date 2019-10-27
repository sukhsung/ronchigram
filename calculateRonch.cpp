
#include <math.h>
#include <stdint.h>
#include <iostream>
#include <stdlib.h>



extern "C" {

    float* deepCopyFloat(float* base, int dimX, int dimY)
    {
        int sz = dimX*dimY;
        float* deepCopy = new float[sz];
        for(size_t i = 0; i < sz; i++)
            deepCopy[i] = base[i]; // copy the allocated memory 
        return deepCopy;
    }


    int sub2ind(int x, int y, int z, int dimX, int dimY, int dimZ) {

        return dimX*dimY*z+dimX*y+x;
    }


    float getMin(float* base, int dimX, int dimY)
    {
        float minV = base[sub2ind(0,0,0,dimX,dimY,1)];
        for(int j=0; j<dimY; j++)
        {
            for(int i = 0; i < dimX; i++)
            {
                float trial = base[sub2ind(i,j,0,dimX,dimY,1)];
                if( trial < minV)
                {
                    minV = trial;
                }
            }
        }
        return minV;
    }

    float getMax(float* base, int dimX, int dimY)
    {
        float maxV = base[sub2ind(0,0,0,dimX,dimY,1)];
        for(int j=0; j<dimY; j++)
        {
            for(int i = 0; i < dimX; i++)
            {
                float trial = base[sub2ind(i,j,0,dimX,dimY,1)];
                if( trial > maxV)
                {
                    maxV = trial;
                }
            }
        }
        return maxV;
    }

    float* normalize(float* base, float scale, int dimX, int dimY)
    {
        float minV = getMin(base, dimX, dimY);
        for(int j=0; j<dimY; j++)
        {
            for(int i = 0; i < dimX; i++)
            {
                int idx = sub2ind(i,j,0,dimX,dimY,1);
                base[idx] = base[idx] - minV;
            }
        }
        float maxV = getMax(base, dimX, dimY);
        for(int j=0; j<dimY; j++)
        {
            for(int i = 0; i < dimX; i++)
            {
                int idx = sub2ind(i,j,0,dimX,dimY,1);
                base[idx] = base[idx] / maxV * scale;
            }
        }
    return base;

    }

    float* linspace (float start, float stop, int steps){
        float* vals = new float[steps];
        for(int i = 0; i < steps; i++)
        {
            vals[i] = start + (stop-start)/(steps-1)*i;
        }

        return vals;
    }

    float* meshgrid(float* x, float* y, int sizeX, int sizeY)
    {
        float* vals = new float[sizeX*sizeY*2];
        for(int j=0; j<sizeX; j++)
        {
            for(int i=0; i<sizeY; i++)
            {
                vals[sub2ind(i,j,0,sizeX,sizeY,2)] = x[i];
                vals[sub2ind(i,j,1,sizeX,sizeY,2)] = y[j];
            }
        }
        return vals;
    }

    float* cValArray(float val, int dimX, int dimY) {
        float* vals = new float[dimX*dimY];
        for(int j=0; j<dimX; j++)
        {
            for(int i=0; i<dimY; i++)
            {
                vals[sub2ind(i,j,0,dimX,dimY,1)] = val;
            }
        }
        return vals;
    }

    float* inPlaceAperture(float* base, float* alrr, float limit, int dimX, int dimY)
    {
        for(int j=0; j<dimY; j++)
        {
            for(int i = 0; i < dimX; i++)
            {
                //float x = mesh[sub2ind(i,j,0,dimX,dimY,2)];
                //float y = mesh[sub2ind(i,j,1,dimX,dimY,2)];
                if(alrr[sub2ind(i,j,0,dimX,dimY,1)]> limit)
                {
                    base[sub2ind(i,j,0,dimX,dimY,1)] = 0;
                }
            }
        }
        return base;
    }

    float* meshToRadii(float* cart_mesh, int dimX, int dimY)
    {
        float* vals = new float[dimX*dimY];
        for(int j=0; j<dimX; j++)
        {
            for(int i=0; i<dimY; i++)
            {
                int idxX = sub2ind(i,j,0,dimX,dimY,1);
                int idxY = sub2ind(i,j,1,dimX,dimY,1);
                float vX = cart_mesh[idxX];
                float vY = cart_mesh[idxY];
                vals[sub2ind(i,j,0,dimX,dimY,1)] =sqrt( vX*vX+vY*vY);
            }
        }
        return vals;

    }
    float* meshToAngles(float* cart_mesh, int dimX, int dimY)
    {
        float* vals = new float[dimX*dimY];
        for(int j=0; j<dimX; j++)
        {
            for(int i=0; i<dimY; i++)
            {
                int idxX = sub2ind(i,j,0,dimX,dimY,1);
                int idxY = sub2ind(i,j,1,dimX,dimY,1);
                float vX = cart_mesh[idxX];
                float vY = cart_mesh[idxY];
                vals[sub2ind(i,j,0,dimX,dimY,1)] =atan2(vY,vX);
            }
        }
        return vals;    
    }



    float* calcRonch(int numPx,float al_max, float objApR) {
        float obj_ap_r = objApR; //mrad
        float simdim = al_max; //mrad


        //numPx x numPx x 2 meshgrid, index into w/ sub2ind
        float* cart_mesh = meshgrid(linspace(-simdim,simdim,numPx),linspace(-simdim,simdim,numPx),numPx,numPx);
        float* oapp = cValArray(1,numPx,numPx);
        float* alrr = meshToRadii(cart_mesh,numPx,numPx);
        float* alpp = meshToAngles(cart_mesh,numPx,numPx);
        float values[numPx][numPx];

        oapp = inPlaceAperture(oapp,alrr,obj_ap_r,numPx,numPx);

        float* res = alpp;
        oapp = normalize(oapp, 255, numPx, numPx);
        for (int j=0; j<numPx; j++) {
            for (int i=0; i<numPx; i++) {

                values[i][j] = res[sub2ind(i,j,0,numPx,numPx,1)];
            }
        }
        auto arrayPtr = &values[0][0];

        //delete res;
        return arrayPtr;

    }


}
