#include <math.h>
#include <stdint.h>
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <time.h>
#include <emscripten.h>
#include "kiss_fftnd.h"
using namespace std;
extern "C" {
int sub2ind(int x, int y, int z, int dimX, int dimY, int dimZ) {
    return dimX * dimY * z + dimX * y + x;
}

float* getMinMax(float* base, int dimX, int dimY) {
    float* minMax = new float[2];
    minMax[0] = base[0]; //min
    minMax[1] = base[0]; //max

    for (int i = 0; i < dimX * dimY; i++) {
        if (base[i] < minMax[0]) {
            minMax[0] = base[i];
        } else if (base[i] > minMax[1]) {
            minMax[1] = base[i];
        }
    }
    return minMax;
}

float* normalize(float* base, float scale, int dimX, int dimY) {
    float* minMax = getMinMax(base, dimX, dimY);
    for (int i = 0; i < dimX * dimY; i++) {
        base[i] = (base[i] - minMax[0]) / minMax[1] * scale;
    }
    return base;
}
float calculateLambda(float keV) {
    float lambda = 12.3986 / sqrt((2 * 511 + keV) * keV) * 1e-10;
    return lambda;
}

int polarMeshnOapp(float* rr, float* pp, float* oapp, float r_max, float obj_ap_r, int numPx) {
    float center = numPx / 2;
    int idx;
    float xval;
    float yval;
    for (int j = 0; j < numPx; j++) {
        for (int i = 0; i < numPx; i++) {
            idx = i + numPx * j;
            xval = (i - center) / center * r_max;
            yval = (j - center) / center * r_max;
            rr[idx] = sqrt(pow(xval, 2) + pow(yval, 2));
            pp[idx] = atan2(yval, xval);
            if (rr[idx] > obj_ap_r) {
                oapp[idx] = 0;
            } else {
                oapp[idx] = 1;
            }
        }
    }
    return 0;
}

float* noisyGrating(int dimX, int dimY) {
    // srand(time(NULL));
    //srandom(time(NULL));
    float* vals = new float[dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        //vals[i] = rand() / float(RAND_MAX);
        //vals[i] = random();
        vals[i] = emscripten_random();
    }
    return vals;
}

float* generateSample(int dimX, int dimY, int scaleFactor) {
    float* subsample = noisyGrating(dimX / scaleFactor, dimX / scaleFactor);
    float* supersample = new float[dimX * dimY];
    for (int j = 0; j < dimY; j++) {
        for (int i = 0; i < dimX; i++) {
            int idx = sub2ind(i, j, 0, dimX, dimY, 1);
            int idxsub = sub2ind(i / scaleFactor, j / scaleFactor, 0, dimX / scaleFactor, dimY / scaleFactor, 1);
            supersample[idx] = subsample[idxsub];
        }
    }
    return supersample;
}

complex<float>* generateTransmissionFn(float* sample, int dimX, int dimY, float interactionParam) {
    complex<float>* trans = new complex<float> [dimX * dimY];
    complex<float> imag(0.0, 1.0);
    for (int i = 0; i < dimX * dimY; i++) {
        complex<float> real_part(M_PI / 4 * sample[i] * interactionParam, 0.0);
        trans[i] = exp(-imag * real_part);
    }
    return trans;
}

float* complexToReal(complex<float>* orig, int dimX, int dimY) {
    float* real = new float[dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        real[i] = abs(orig[i]);
    }
    return real;
}

complex<float>* realToComplex(float* orig, int dimX, int dimY) {
    complex<float>* comp = new complex<float> [dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        comp[i].real(orig[i]);
        comp[i].imag(0);
    }

    return comp;
}

kiss_fft_cpx* complexToKiss(complex<float>* orig, int dimX, int dimY) {
    kiss_fft_cpx* kisscpx = new kiss_fft_cpx[dimX * dimY];
    kiss_fft_cpx holder;
    for (int i = 0; i < dimX * dimY; i++) {
        holder.r = real(orig[i]);
        holder.i = imag(orig[i]);
        kisscpx[i] = holder;
    }
    return kisscpx;
}

complex<float>* kissToComplex(kiss_fft_cpx* orig, int dimX, int dimY) {
    complex<float>* comp = new complex<float> [dimX * dimY];
    for (int i = 0; i < dimX * dimY; i++) {
        //comp[i] = (orig[i].r,orig[i].i);
        comp[i].real(orig[i].r);
        comp[i].imag(orig[i].i);
    }
    return comp;
}

float* packageOutput(float* base1, float* base2, float* base3, float* scalars, int dimX, int dimY, int nScal) {
    int sz = dimX * dimY;
    float* imageStack = new float[sz * 2 + nScal];
    for (size_t i = 0; i < sz; i++) {
        imageStack[i] = base1[i]; // copy the allocated memory
        imageStack[sz + i] = base2[i];
        imageStack[2 * sz + i] = base3[i];
    }
    for (size_t i = 0; i < nScal; i++) {
        imageStack[3 * sz + i] = scalars[i];
    }
    return imageStack;
}

float* calculateChi0(float* magptr, float* angleptr, float* alrr, float* alpp, int numPx, int numAbs, float keV) {
    float* chi0 = new float[numPx * numPx];
    int n[14] = {1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5};
    int m[14] = {0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, 2, 4, 6};
    float lambda = calculateLambda(keV);

    for (int i = 0; i < numPx*numPx; i++) {
        chi0[i] = 0;
        for (int k = 0; k < numAbs; k++) {
            chi0[i] = chi0[i] + 2 * M_PI / lambda * magptr[k] * pow(alrr[i], n[k] + 1) * cos(m[k] * (alpp[i] - angleptr[k])) / (n[k] + 1);
        }
    }
    return chi0;
}

complex<float>* calculateChi(float* chi0, int numPx) {
    complex<float>* chi = new complex<float> [numPx * numPx];
    complex<float> imag(0.0, 1.0);
    for (int i = 0; i < numPx * numPx; i++) {
        chi[i] = exp(-imag * chi0[i]);
    }
    return chi;
}

float maskChi0(float* chi0, float* alrr, int numPx, float threshold) {
    float rmax = 1e5;
    for (int i = 0; i < numPx * numPx; i++) {
        chi0[i] = 255 * ((abs(chi0[i]) > threshold) ? 0 : 1);
        float cv = (1 - chi0[i]) * alrr[i];
        if (cv > 0 && cv < rmax) {
            rmax = cv;
        }
    }
    return rmax * 1000;
}

float* fftShift(float* original, int numPx) {
    // only for square matrices...
    float* shifted = new float[numPx * numPx];

    for (int j = 0; j < numPx / 2; j++) {
        for (int i = 0; i < numPx / 2; i++) {
            int idx_q1 = sub2ind(i, j, 0, numPx, numPx, 1);
            int idx_q3 = sub2ind(i + numPx / 2, j + numPx / 2, 0, numPx, numPx, 1);
            int idx_q2 = sub2ind(i + numPx / 2, j, 0, numPx, numPx, 1);
            int idx_q4 = sub2ind(i, j + numPx / 2, 0, numPx, numPx, 1);
            shifted[idx_q1] = original[idx_q3];
            shifted[idx_q3] = original[idx_q1];
            shifted[idx_q4] = original[idx_q2];
            shifted[idx_q2] = original[idx_q4];
        }
    }
    return shifted;
}

complex<float>* cmplxFFT(complex<float>* comp, int dimX, int dimY) {
    int isInverseFFT = 0;
    int ndims = 2;
    int dims[2];
    dims[0] = dimX;
    dims[1] = dimY;

    kiss_fft_cpx* cxin;
    kiss_fft_cpx* cxout;

    kiss_fftnd_cfg cfg = kiss_fftnd_alloc(dims, ndims, isInverseFFT, 0, 0);
    // preallocation
    cxin = complexToKiss(comp, dimX, dimY);
    cxout = complexToKiss(comp, dimX, dimY);
    kiss_fftnd(cfg, cxin, cxout);
    complex<float>* fftResult = kissToComplex(cxout, dimX, dimY);

    return fftResult;
}

float* calcDiffract(complex<float>* chi, complex<float>* trans, float* oapp, int numPx) {
    // want: abs(fft(trans*fft(chi)))

    //in place b/c chi won't be reused
    chi = cmplxFFT(chi, numPx, numPx);
    for (int i = 0; i < numPx * numPx; i++) {
        chi[i] = chi[i] * trans[i];
    }
    chi = cmplxFFT(chi, numPx, numPx);
    float* diffInt = new float[numPx * numPx];
    diffInt = complexToReal(chi, numPx, numPx);
    for (int i = 0; i < numPx * numPx; i++) {
        diffInt[i] = diffInt[i] * oapp[i] * diffInt[i];
    }
    return diffInt;
}

float calculateInteractionParam(float keV) {
    float keV_300 = 300;
    float c = 3e8;
    float mass_e = 9.11e-31;
    float charge_e = 1.602e-19;
    float lambda = calculateLambda(keV);
    float lambda_300 = calculateLambda(keV_300);
    float param = 2 * M_PI / (lambda * keV / charge_e * 1000) * (mass_e * c * c + keV * 1000) / (2 * mass_e * c * c + keV * 1000);
    float param_300 = 2 * M_PI / (lambda_300 * keV_300 / charge_e * 1000) * (mass_e * c * c + keV_300 * 1000) / (2 * mass_e * c * c + keV_300 * 1000);
    return param / param_300;
}

float singleStrehl(float rmax, complex<float>* chi, float al_max, int numPx){
    float strr[numPx * numPx];
    float stpp[numPx * numPx];
    float sapp[numPx * numPx];
    float  strehl_radius = rmax/1000;
    polarMeshnOapp(strr, stpp, sapp, al_max, strehl_radius, numPx);
    float* strehl_aperture_real = sapp;
    complex<float>* strehl_aperture = realToComplex(strehl_aperture_real, numPx, numPx);
    complex<float>* strehl_inner = new complex<float> [numPx * numPx];
     for(int i = 0; i < numPx * numPx; i++){
         strehl_inner[i] = chi[i] * strehl_aperture[i];
     };
    complex<float>* strehl_inner_fft = cmplxFFT(strehl_inner, numPx, numPx);
    float* strehl_inner_fft_real = complexToReal(strehl_inner_fft, numPx, numPx);
    float strehl_inner_max = getMinMax(strehl_inner_fft_real, numPx, numPx)[1];
    complex<float>*strehl_bottom = cmplxFFT(strehl_aperture, numPx, numPx);
    float* strehl_bottom_real = complexToReal(strehl_bottom, numPx, numPx);
    float strehl_bottom_max = getMinMax(strehl_bottom_real, numPx, numPx)[1];
    float strehl = strehl_inner_max/strehl_bottom_max;
    strehl = pow(strehl,2);
    return strehl;
}

float* probeGeneration(complex<float>* chi, int numPx, float* oapp){
    complex<float>* obj_aperture = realToComplex(oapp, numPx, numPx);
    complex<float>* probe = new complex<float> [numPx * numPx];
     for(int i = 0; i < numPx * numPx; i++){
         probe[i] = chi[i] * obj_aperture[i];
     };
    complex<float>* probe_fft = cmplxFFT(probe, numPx, numPx);
    float* probe_abs = complexToReal(probe_fft, numPx, numPx);
    float* probe_shift = fftShift(probe_abs, numPx);
    float* probe_out = normalize(probe_shift, 255, numPx, numPx);
    return probe_out;
}

float* calcRonch(float* buffer, int bufSize) {
    int numPx = static_cast < int > (buffer[0]);
    float al_max = buffer[1]; //mrad
    float obj_ap_r = buffer[2]; //mrad
    int scalefactor = buffer[3];
    float keV = buffer[4];
    bool calcStrehl = buffer[5];
    float* outputScalars = new float[5];
    float alrr[numPx * numPx];
    float alpp[numPx * numPx];
    float oapp[numPx * numPx];
    polarMeshnOapp(alrr, alpp, oapp, al_max, obj_ap_r, numPx);

    float* sample = generateSample(numPx, numPx, scalefactor);

    complex<float>* trans = generateTransmissionFn(sample, numPx, numPx, calculateInteractionParam(keV));

    float* chi0 = calculateChi0( & buffer[6], & buffer[20], alrr, alpp, numPx, 14, keV);

    complex<float>* chi = calculateChi(chi0, numPx);
    // Calculate r_max and return,  and turn chi0 into pi/4map Normalized
    outputScalars[0] = maskChi0(chi0, alrr, numPx, M_PI / 4);
    
    float strehl = singleStrehl(outputScalars[0], chi, al_max, numPx);
    outputScalars[1] = strehl;


    float r_strehl;
    float middle;
    float count = 0;
    if(calcStrehl == true){
    float upper_bound = 1.8 * outputScalars[0]/1000;
    float lower_bound = 1 * outputScalars[0]/1000;
    float upper_strehl = 1;
    float lower_strehl = 0.5;
    float tolerance = 0.01;
    float ans;
    while((upper_strehl-lower_strehl)>=tolerance){
    

        //The "middle" variable is set to the average of the bounds with each iteration
        middle = (upper_bound+lower_bound)/2;

        //Find the strehl ratio for the given radius
        ans = singleStrehl(1000*middle, chi, al_max, numPx);
    
        if(ans > 0.8){
            lower_bound = middle;
            upper_strehl = ans;
        }

        if(ans < 0.8){
            upper_bound = middle;
            lower_strehl = ans;
        }

        // if(ans == 0.8){
        //     break;
        // }
        count = count + 1;
        if(count == 10){
            break;
        }
    }

    middle = (upper_bound+lower_bound)/2;
    r_strehl = singleStrehl(1000*middle, chi, al_max, numPx);
    }
    else{
        r_strehl = -1;
        middle = 0;
        count = 0;
    }
    outputScalars[2] = r_strehl; //strehl ratio achieved
    outputScalars[3] = middle*1000; //strehl radius
    outputScalars[4] = count;

    float strr[numPx * numPx];
    float stpp[numPx * numPx];
    float sapp[numPx * numPx];
    float  strehl_radius = middle;
    
    if(calcStrehl == true){
    polarMeshnOapp(strr, stpp, sapp, al_max, strehl_radius, numPx);
    }
    else{
    polarMeshnOapp(strr, stpp, sapp, al_max, outputScalars[0]/1000, numPx);
    
    }
    float* probe = probeGeneration(chi, numPx, sapp);

    // Normalize to 0-255 for output
    float* ronch = normalize(calcDiffract(chi, trans, oapp, numPx), 255, numPx, numPx);
    // Package results and return
    auto arrayPtr = packageOutput(ronch, chi0, probe, outputScalars, numPx, numPx, 5);
    return arrayPtr;
}

}