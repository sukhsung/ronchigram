function fft2_wrap(X) {
    let X_pass = fft2(X);
    let X_arr = math.matrix(X_pass);
    X_arr = math.transpose(X_arr);
    X_arr = X_arr.toArray();
    X_pass = fft2(X_arr);
    return X_pass;
}
/*
    From https://gist.github.com/mrquincle/b11fff96209c9d1396b0
    @mrquincle
*/
function fft2(X) {
    let N = X.length;
    if (!(N > 1)) {
        return X;
    }
    let M = N / 2;
    let even = [];
    let odd = [];
    even.length = M;
    odd.length = M;
    for (let i = 0; i < M; ++i) {
        even[i] = X[i * 2];
        odd[i] = X[i * 2 + 1];
    }
    even = fft2(even);
    odd = fft2(odd);
    let a = -2 * PI;
    for (let k = 0; k < M; ++k) {
        let t = math.exp(math.complex(0, (a * k) / N));
        t = math.multiply(t, odd[k]);
        X[k] = odd[k] = math.add(even[k], t);
        X[k + M] = even[k] = math.subtract(even[k], t);
    }
    return X;
}

function normalizeScale(matharray, scale) {
    let output = math.subtract(matharray, math.min(matharray));
    output = math.dotDivide(output, math.max(output) / scale);
    return output;
}

function drawOverlays(
    ctx1,
    ctx2,
    numPx,
    al_max,
    disp_size_mrad,
    obj_ap_r,
    rmax,
    best_r,
    ctx3,
    obj_ap_r
) {
    let scalar = 256;
    let keV = Number(document.getElementById("beamvolt").value);
    let lambda = (12.3986 / Math.sqrt((2 * 511 + keV) * keV)) * ang;

    ctx1.font = (numPx / scalar) * 14 + "px Arial";
    ctx1.fillStyle = "#e6eaeb";
    ctx1.fillText(
        math.round((disp_size_mrad / 0.07) * 30) + " mrad",
        numPx - (70 / scalar) * numPx,
        numPx - (10 / scalar) * numPx
    );

    ctx1.beginPath();
    ctx1.moveTo(numPx - (70 / scalar) * numPx, numPx - (30 / scalar) * numPx);
    ctx1.lineTo(numPx - (15 / scalar) * numPx, numPx - (30 / scalar) * numPx);
    ctx1.strokeStyle = "#e6eaeb";
    ctx1.lineWidth = (5 * numPx) / scalar;
    ctx1.stroke();
    ctx1.beginPath();
    ctx1.arc(
        numPx / 2,
        numPx / 2,
        ((rmax * numPx) / (2 * al_max)) * mrad,
        0,
        2 * PI
    );
    ctx1.strokeStyle = "#71B6D6"; //from "blue"
    ctx1.lineWidth = (1 * numPx) / scalar;
    ctx1.stroke();

    ctx1.beginPath();
    ctx1.arc(
        numPx/2,
        numPx/2,
        ((best_r * numPx)/(2 * al_max)) * mrad,
        0,
        2 * PI
    );
    ctx1.strokeStyle = "#B6D671"; //for the strehl radius
    ctx1.stroke();

    /// right panel
    ctx2.beginPath();
    ctx2.arc(
        numPx / 2,
        numPx / 2,
        ((rmax * numPx) / (2 * al_max)) * mrad,
        0,
        2 * PI
    );
    ctx2.strokeStyle = "#71B6D6"; //from "blue"
    ctx2.lineWidth = 2;
    ctx2.stroke();
    ctx2.beginPath();
    ctx2.arc(
        numPx / 2,
        numPx / 2,
        (obj_ap_r * numPx) / (2 * al_max),
        0,
        2 * PI
    );
    ctx2.strokeStyle = "#D671B6"; //from "red"
    ctx2.lineWidth = 2;
    ctx2.stroke();

    ctx2.beginPath();
    ctx2.arc(
        numPx/2,
        numPx/2,
        ((best_r * numPx)/(2 * al_max)) * mrad,
        0,
        2 * PI
    );
    ctx2.strokeStyle ="#B6D671"; //for the strehl radius
    ctx2.stroke();

    ctx3.font = (numPx / scalar) * 5 + "px Arial"; // this is where to put in the correct length scale
    ctx3.fillStyle = "#e6eaeb";
    ctx3.fillText(
         math.round((numPx*lambda)/(2*obj_ap_r/1000)*(39/(2*numPx))*10000000,1) + " Å",
         numPx/4 - (55 / scalar) * numPx/4,
         numPx/4 - (10 / scalar) * numPx/4
    );
    ctx3.beginPath();
    ctx3.moveTo(numPx/4 - (100 / scalar) * numPx/4, numPx/4 - (40 / scalar) * numPx/4);
    ctx3.lineTo(numPx/4 - (22 / scalar) * numPx/4, numPx/4 - (40 / scalar) * numPx/4);
    ctx3.strokeStyle = "#e6eaeb";
    ctx3.lineWidth = (1 * numPx) / scalar;
    ctx3.stroke();
}

function drawGrayscaleBitmap2( ctx, input, numPx ) {
    var imData = ctx.getImageData(0,0,numPx,numPx)
    var data = imData.data;

    let red
    for (let i = 0; i < numPx; i++ ){
        for (let j = 0; j <numPx; j++ ){
            red = i * (numPx * 4) + j * 4

            data[red] = input[i][j]
            data[red+1] = input[i][j]
            data[red+2] = input[i][j]
            data[red+3] = 255
        }

    }

    ctx.putImageData(imData,0,0)

}

function getAberrations() {
    let ab_mags = [];
    let ab_angles = [];
    let abers = [];
    let numAbers = aberrations.length;
    for (let it = 0; it < numAbers; it++) {
        let aberration = aberrations[it];
        let mag_val = Number(aberration.mag_el.value) * aberration.mag_unit;
        let arg_val =
            (aberration.arg_el ? Number(aberration.arg_el.value) : 0) * deg;
        ab_mags.push(mag_val);
        ab_angles.push(arg_val);
    }
    abers.push(ab_mags);
    abers.push(ab_angles);
    abers.push(numAbers);
    return abers;
}

function lambdaCalc(keV) {
    return (12.3986 / Math.sqrt((2 * 511 + keV) * keV)) * ang;
}

function energyUI() {
    let keV = Number(document.getElementById("beamvolt").value);

    if (keV < 0) {
        keV = 0;
        document.getElementById("beamvolt").value = 0;
    }

    let lambda = lambdaCalc(keV);
    document.getElementById("wavlen").value = math.round(lambda / pm, 4);

    let alpha = Number(document.getElementById("aperture").value) * mrad;
    //resolution calculation:
    let d = (0.61 * lambda) / pm / alpha;
    document.getElementById("diffres").value = math.round(d, 4);
    return lambda;
}

function randButton() {
    //document.getElementById('loading').innerHTML = "Calculating..."
    if(real.checked){
        setTimeout(function () {
            randomize_real();
        }, 0);
    } else {
        setTimeout(function () {
            randomize_realistic(); //randomize_realistic
        }, 0);
    }
}

function generateSample(scalefactor, numPx) {
    let subsample = math.random([numPx / scalefactor, numPx / scalefactor]);
    let supersample = math.zeros(numPx, numPx);
    //quick nearest neightbours interpolation
    supersample = supersample.map(function (value, index, matrix) {
        return subsample[
            math.floor(index[0] / scalefactor)
        ][math.floor(index[1] / scalefactor)];
    });
    return supersample;
}

//Normalized to 300 keV
function calculateInteractionParam(keV) {
    let keV_300 = 300;
    let c = 3e8;
    let mass_e = 9.11e-31;
    let charge_e = 1.602e-19;
    let lambda = lambdaCalc(keV);
    let lambda_300 = lambdaCalc(keV_300);
    let param =
        (((2 * PI) / (((lambda * keV) / charge_e) * 1000)) *
            (mass_e * c * c + keV * 1000)) /
        (2 * mass_e * c * c + keV * 1000);
    let param_300 =
        (((2 * PI) / (((lambda_300 * keV_300) / charge_e) * 1000)) *
            (mass_e * c * c + keV_300 * 1000)) /
        (2 * mass_e * c * c + keV_300 * 1000);
    return param / param_300;
}

function getObjAperture() {
    let obj_ap_r = Number(document.getElementById("aperture").value) * mrad;

    if (obj_ap_r < 0) {
        obj_ap_r = 0;
        document.getElementById("aperture").value = 0;
    }
    return obj_ap_r;
}

function getDispSizePx() {
    let disp_size_px = Number(document.getElementById("disp_size_px").value);
    if ((disp_size_px & (disp_size_px - 1)) != 0 || disp_size_px < 2) {
        alert(
            "Select a display size in pixels that is a power of 2 greater than 0"
        );
        return;
    } else {
        return disp_size_px;
    }
}

function getDispSizeMrad() {
    let disp_size_mrad =
        (Number(document.getElementById("disp_size_mrad").value) * mrad) / 2;
    if (disp_size_mrad < 0.0000001) {
        disp_size_mrad = 0.0000001;
    }
    return disp_size_mrad;
}

function polarMeshOapp(r_max, obj_ap_r, numPx) {
    const center = numPx / 2;
    let idx;
    let xval;
    let yval;
    let rval;

    let curRR = [];
    let curPP = [];
    let curOb = [];

    let rr = [];
    let pp = [];
    let oapp = [];
    let rrPPObj = [];
    for (let j = 0; j < numPx; j++) {
        for (let i = 0; i < numPx; i++) {
            xval = ((i - center) / center) * r_max;
            yval = ((j - center) / center) * r_max;
            rval = Math.sqrt(math.pow(xval, 2) + math.pow(yval, 2));
            curRR.push(rval);
            curPP.push(math.atan2(yval, xval));
            curOb.push(rval > obj_ap_r ? 0 : 1);
        }
        rr.push(curRR);
        pp.push(curPP);
        oapp.push(curOb);
        curRR = [];
        curPP = [];
        curOb = [];
    }
    rrPPObj.push(rr);
    rrPPObj.push(pp);
    rrPPObj.push(oapp);
    return rrPPObj;
}

function hasWASM() {
    //per @JF-Bastien https://stackoverflow.com/questions/47879864/how-can-i-check-if-a-browser-supports-webassembly
    const supported = (() => {
        try {
            if (
                typeof WebAssembly === "object" &&
                typeof WebAssembly.instantiate === "function"
            ) {
                const module = new WebAssembly.Module(
                    Uint8Array.of(0x0, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00)
                );
                if (module instanceof WebAssembly.Module)
                    return (
                        new WebAssembly.Instance(module) instanceof
                        WebAssembly.Instance
                    );
            }
        } catch (e) {
            console.log(e);
        }
        return false;
    })();
    if (!supported) {
        console.log("WASM not supported, Foricng Javascript");
    } else {
        console.log("WASM supported");
    }
    return supported;
}

function calcButton() {
    calculate();
}

function maskChi0(chi0, alrr, numPx) {
    let rmax = 1e5;
    let maskedChi0 = [];
    let curMasked = [];
    let cv;
    for (let i = 0; i < numPx; i++) {
        for (let j = 0; j < numPx; j++) {
            cv = math.abs(chi0[i][j]) > PI / 4 ? 0 : 1;
            curMasked.push(cv * 255);
            cv = (1 - cv) * alrr[i][j];
            if (cv > 0 && cv < rmax) {
                rmax = cv;
            }
        }
        maskedChi0.push(curMasked);
        curMasked = [];
    }
    let returns = [];
    returns.push(maskedChi0);
    returns.push(rmax * 1000);
    return returns;
}

function calculateChi0(mags, angles, alrr, alpp, numPx, numAbers, keV) {
    let n = [1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5];
    let m = [0, 2, 1, 3, 0, 2, 4, 1, 3, 5, 0, 2, 4, 6];
    let lambda = lambdaCalc(keV);

    let curRow = [];
    let curVal;
    let chi0 = [];

    for (let i = 0; i < numPx; i++) {
        for (let j = 0; j < numPx; j++) {
            curVal = 0;
            for (let k = 0; k < numAbers; k++) {
                curVal =
                    curVal +
                    (((2 * PI) / lambda) *
                        mags[k] *
                        math.pow(alrr[i][j], n[k] + 1) *
                        math.cos(m[k] * (alpp[i][j] - angles[k]))) /
                        (n[k] + 1);
            }
            curRow.push(curVal);
        }
        chi0.push(curRow);
        curRow = [];
    }
    return chi0;
}

function probeZoom(probe, numPx ){
    let len = numPx/4;
    let probe_new = Array(len).fill().map(() => Array(len).fill(0));
    for(let i = 0; i<numPx; i++){
        for(let j = 0; j<numPx; j++){
            if(i>=(3*numPx/8) && i<(5*numPx/8) && j>=(3*numPx/8) && j<(5*numPx/8)){
                probe_new[i-3*numPx/8][j-3*numPx/8] = probe[i][j];
            };
        };
    };
    return probe_new;
}

function drawEverything(
    out_ronch,
    out_phase_map,
    numPx,
    al_max,
    disp_size_mrad,
    obj_ap_r,
    rmax,
    out_probe,
    best_strehl,
    best_r,
    current_strehl,
    count,
    draw_overlay
) {
    // drawing part
    canvas1.width = numPx;
    canvas1.height = numPx;
    canvas2.width = numPx;
    canvas2.height = numPx;
    canvas3.width = numPx/4; //if changing probe zoom, need to change this as well
    canvas3.height = numPx/4;

    drawGrayscaleBitmap2(ctx1, out_ronch, numPx);
    drawGrayscaleBitmap2(ctx2, out_phase_map, numPx);
    // if (draw_overlay) {
    //     drawOverlays(ctx1, ctx2, numPx, al_max, disp_size_mrad, obj_ap_r, rmax, best_r);
    // }
    let probe = probeZoom(out_probe, numPx);

    drawGrayscaleBitmap2(ctx3, probe, probe.length);
    
    if (draw_overlay) {
        drawOverlays(ctx1, ctx2, numPx, al_max, disp_size_mrad, obj_ap_r, rmax, best_r, ctx3, obj_ap_r);
    }
    document.getElementById("alpha_max").value = math.round(rmax, 2);
    //document.getElementById("best_strehl").value = math.round(best_strehl, 3);
    document.getElementById("current_strehl").value = math.round(current_strehl, 3);
    document.getElementById("best_r").value = math.round(best_r, 2);
    //document.getElementById("ratio").value = math.round(best_r/rmax, 4);
    //document.getElementById("count").value = count;
}

function rotate(src, n) {
    var len = src.length
    reverse(src, 0, len)
    reverse(src, 0, n)
    reverse(src, n, len)
    return src
}

function reverse(src, from, to) {
    --from
    while (++from < --to) {
      var tmp = src[from]
      src[from] = src[to]
      src[to] = tmp
    }
}

function ifftshift(src) {
    const len = src.length
    return rotate(src, math.floor((len+1) / 2))
}

function fftshift(src){
    const len = src.length;
    let src_new = Array(len).fill().map(() => Array(len).fill(0));
    for(let i = 0; i<len; i++){
        for(let j = 0; j<len; j++){
            if(i<len/2 && j<len/2){
                src_new[i+len/2][j+len/2] = src[i][j];
            }
            if(i<len/2 && j>=len/2){
                src_new[i+len/2][j-len/2] = src[i][j];
            }
            if(i>=len/2 && j<len/2){
                src_new[i-len/2][j+len/2] = src[i][j];
            }
            if(i>=len/2 && j>=len/2){
                src_new[i-len/2][j-len/2] = src[i][j];
            }
        }
    }
    return src_new;
}

function singleStrehl(rmax, chi, al_max, numPx){
    let  strehl_radius = rmax/1000;
    let rrStrehl = polarMeshOapp(al_max, strehl_radius, numPx);
    let strehl_aperture = rrStrehl[2];
    let strehl_inner = fft2_wrap((math.dotMultiply(chi.toArray(),strehl_aperture)));
    strehl_inner = math.max(math.abs(strehl_inner));
    let strehl_bottom = fft2_wrap(strehl_aperture);
    strehl_bottom = math.max(math.abs(strehl_bottom));
    let strehl = strehl_inner/strehl_bottom;
    strehl = math.pow(strehl,2);
    return strehl;
}

function pointEightStrehl(rmax, chi, al_max, numPx){
    let best_strehl = -1;
    let best_r = -1;
    let upper_bound = 1.8 * rmax/1000;
    let lower_bound = 1 * rmax/1000;
    let upper_strehl = 1;
    let lower_strehl = 0.5;
    let tolerance = 0.01;
    let middle;
    let ans;
    let count = 0;
    while((upper_strehl-lower_strehl)>=tolerance){
        middle = (upper_bound+lower_bound)/2;
        ans = singleStrehl(1000*middle, chi, al_max, numPx);
        if(ans > 0.8){
            lower_bound = middle;
            upper_strehl = ans;
        }
        if(ans < 0.8){
            upper_bound = middle;
            lower_strehl = ans;
        }
        count = count + 1;
        if(count == 20){
            break;
        }
    }
    middle = (upper_bound+lower_bound)/2;
    best_r = middle*1000;
    best_strehl = singleStrehl(best_r, chi, al_max, numPx);
    let return_array = [best_strehl, best_r, count];
    return return_array;
}

function generateProbe(al_max, best_r, numPx, chi) {
    let probe_ap = polarMeshOapp(al_max, best_r/1000, numPx);
    out_probe = fft2_wrap((math.dotMultiply(chi.toArray(),probe_ap[2])));
    out_probe = math.abs(fftshift(out_probe));
    out_probe = math.matrix(out_probe);
    out_probe = math.round(normalizeScale(out_probe, 255));
    out_probe = out_probe.toArray()
    return out_probe
}

function calculateJS() {
    ////////
    //reading in constants from ui:
    ////////
    let numPx = getDispSizePx();
    let disp_size_mrad = getDispSizeMrad();
    let al_max = disp_size_mrad;
    let obj_ap_r = getObjAperture();
    let keV = Number(document.getElementById("beamvolt").value);
    let scalefactor = Number(
        document.getElementById("sample_scale_factor").value
    );

    let draw_overlay = document.getElementById("draw_overlay").checked; 
    let rrPPObj = polarMeshOapp(al_max, obj_ap_r, numPx);
    let alrr = rrPPObj[0];
    let alpp = rrPPObj[1];
    let obj_ap = math.matrix(rrPPObj[2]);

    let sample = generateSample(scalefactor, numPx);

    let trans = math.exp(
        math.multiply(
            math.complex(0, -1),
            PI,
            0.25,
            calculateInteractionParam(keV),
            sample
        )
    );

    let abers = getAberrations();
    let ab_mags = abers[0];
    let ab_angles = abers[1];
    let numAbers = abers[2];

    let chi0 = calculateChi0(
        ab_mags,
        ab_angles,
        alrr,
        alpp,
        numPx,
        numAbers,
        keV
    );

    let chi = math.dotPow(
        math.E,
        math.dotMultiply(math.complex(0, -1), math.matrix(chi0))
    );
    // //To place objective before sample:
    //var chi = math.dotMultiply(chi, obj_ap);

    let out_probe = math.matrix(fft2_wrap(chi.toArray()));
    let out_ronch = math.dotPow(math.abs(math.dotMultiply(math.matrix(fft2_wrap(math.dotMultiply(trans, out_probe).toArray())),obj_ap)),2);

    out_ronch = math.round(normalizeScale(out_ronch, 255));
    out_ronch = out_ronch.toArray();

    let returnVals = maskChi0(chi0, alrr, numPx);
    let out_phase_map = returnVals[0];
    let rmax = returnVals[1];
    let current_strehl = singleStrehl(rmax, chi, al_max, numPx);    
    let point_eight_strehl = pointEightStrehl(rmax, chi, al_max, numPx);
    out_probe = generateProbe(al_max, point_eight_strehl[1], numPx, chi);

    drawEverything(
        out_ronch,
        out_phase_map,
        numPx,
        al_max,
        disp_size_mrad,
        obj_ap_r,
        rmax,
        out_probe,
        point_eight_strehl[0],
        point_eight_strehl[1],
        current_strehl,
        point_eight_strehl[2],
        draw_overlay
    );
}

function calculateWASM(Module) {
    ////////
    //reading in constants from ui:
    ////////
    let numPx = getDispSizePx();
    let disp_size_mrad = getDispSizeMrad();
    let al_max = disp_size_mrad;
    let obj_ap_r = getObjAperture();
    let keV = Number(document.getElementById("beamvolt").value);
    let scalefactor = Number(
        document.getElementById("sample_scale_factor").value
    );
    let draw_overlay = document.getElementById("draw_overlay").checked; //figure out how to read from checkbox

    // getting aberrations into tidy arrays of magnitude, angle. degree, order  are assumed based on order in C++ section, units are baked in!
    let abers = getAberrations();
    let ab_mags = abers[0];
    let ab_angles = abers[1];
    globalTest = 250;
    let params = [numPx, al_max, obj_ap_r, scalefactor, keV];
    const arrayDataToPass = params.concat(ab_mags, ab_angles);
    let buffer;
    let error;
    let result;
    try {
        const typedArray = new Float32Array(arrayDataToPass.length);
        for (let i = 0; i < arrayDataToPass.length; i++) {
            typedArray[i] = arrayDataToPass[i];
        }
        buffer = Module._malloc(
            typedArray.length * typedArray.BYTES_PER_ELEMENT
        );
        Module.HEAPF32.set(typedArray, buffer >> 2);
        result = Module.ccall(
            "calcRonch",
            null,
            ["number", "number"],
            [buffer, arrayDataToPass.length]
        );
    } catch (e) {
        error = e;
    } finally {
        // To avoid memory leaks we need to always clear out the allocated heap data
        // This needs to happen in the finally block, otherwise thrown errors will stop code execution before this happens
        Module._free(buffer);
    }
    if (error) throw error;

    let arrayData1 = [];
    let arrayData2 = [];
    let arrayData3 = [];
    let out_ronch = [];
    let out_phase_map = [];
    let out_probe = [];
    let im2Offset = numPx * numPx;
    for (let j = 0; j < numPx; j++) {
        for (let i = 0; i < numPx; i++) {
            arrayData1.push(
                Module.HEAPF32[
                    result / Float32Array.BYTES_PER_ELEMENT + i + numPx * j
                ]
            );
            arrayData2.push(
                Module.HEAPF32[
                    result / Float32Array.BYTES_PER_ELEMENT +
                        i +
                        numPx * j +
                        im2Offset
                ]
            );
            arrayData3.push(
                Module.HEAPF32[
                    result / Float32Array.BYTES_PER_ELEMENT +
                        i +
                        numPx * j +
                        im2Offset * 2
                ]
            );
        }
        out_ronch.push(arrayData1);
        out_phase_map.push(arrayData2);
        out_probe.push(arrayData3);
        arrayData1 = [];
        arrayData2 = [];
        arrayData3 = [];
    }
    let rmax =
        Module.HEAPF32[
            result / Float32Array.BYTES_PER_ELEMENT + 3 * (numPx * numPx)
        ];
    let current_strehl = 
        Module.HEAPF32[
            result / Float32Array.BYTES_PER_ELEMENT + 3 * (numPx * numPx) + 1
        ];
    let best_strehl = 
        Module.HEAPF32[
            result / Float32Array.BYTES_PER_ELEMENT + 3 * (numPx * numPx) + 2
        ];
    let best_r = 
        Module.HEAPF32[
            result / Float32Array.BYTES_PER_ELEMENT + 3 * (numPx * numPx) + 3
        ];
    let count = 
        Module.HEAPF32[
            result / Float32Array.BYTES_PER_ELEMENT + 3 * (numPx * numPx) + 4
        ];

    drawEverything(
        out_ronch,
        out_phase_map,
        numPx,
        al_max,
        disp_size_mrad,
        obj_ap_r,
        rmax,
        out_probe,
        best_strehl,
        best_r,
        current_strehl,
        count,
        draw_overlay
    );
}

function calculate() {
    //let t0 = performance.now();
    energyUI();
    if (hasWASM && !forceJS.checked) {
        //document.getElementById("loading").innerHTML =
        //    "Calculating with WebAssembly...";
        //console.log("Calculating with WebAssembly...");
        let curInstance = ronchModule().then(function (Module) {
            calculateWASM(Module);
            Module.delete;
        });
    } else {
        document.getElementById("loading").innerHTML =
            "Calculating with Javascript...";
        console.log("Calculating with Javascript...");
        calculateJS();
        document.getElementById("loading").innerHTML = " ";
    }
    //console.log("T = " + (performance.now() - t0) + " ms");
}

function initialize() {
    hasWASM = hasWASM();
    let url = window.location.href
    let urlparts = url.split('?')
    if (urlparts.length == 1) {
        calculate();
    } 
    else if (urlparts.length == 2) {
        let ronchID = parseInt(urlparts[1])
        if (~isNaN(ronchID)) {
            randomize_realistic(terms=-1, aberration_set_index=ronchID)
        }
    }
    else {
        calculate();
    }
}

function randomize() {
    for (let it = 0; it < aberrations.length; it++) {
        let aberration = aberrations[it];
        aberration.mag_el.value = Math.round(Math.random() * 100);
        if (aberration.arg_el) {
            aberration.arg_el.value = Math.round(Math.random() * 180);
        }
    }
    calculate();
}

function randomize_realistic(terms = -1, aberration_set_index=-1) {
    
    if(terms==-1)
    {
        terms = number_aberration_terms;
    }
    let randVal = 99; //or however large the data set is minus one

    if (aberration_set_index==-1) {
        aberration_set_index = Math.round(Math.random() * randVal);
    }
    console.log(aberration_set_index)
    let aberration_coefs = random_reasonable_aberration_coefs[aberration_set_index];
    for (let it =0; it < terms; it++) {
        let aberration = aberrations[it];
        aberration.mag_el.value = aberration_coefs[2*it];
        if (aberration.arg_el) {
            aberration.arg_el.value = aberration_coefs[2*it+1];   }
    
        }
    calculate();
}

function randomize_real(terms = -1) {
    
    if(terms==-1)
    {
        terms = number_aberration_terms;
    }
    let randVal = 8; //or however large the data set is minus one
    aberration_set_index = Math.round(Math.random() * randVal);
    let aberration_coefs = real_aberration_coefs[aberration_set_index];
    for (let it =0; it < terms; it++) {
        let aberration = aberrations[it];
        aberration.mag_el.value = aberration_coefs[2*it];
        if (aberration.arg_el) {
            aberration.arg_el.value = aberration_coefs[2*it+1];   }
    
        }
    calculate();
}

function allZero() {
    for (let it = 0; it < aberrations.length; it++) {
        let aberration = aberrations[it];
        aberration.mag_el.value = 0;
        if (aberration.arg_el) {
            aberration.arg_el.value = 0;
        }
    }
}

function setC(c_in) {
    for (let it = 0; it < aberrations.length; it++) {
        if (c_in == "" + aberrations[it].m + aberrations[it].n) {
            aberrations[it].mag_el.value =
                Number(aberrations[it].mag_el.value) + 50;
        }
    }
}


var pm = math.pow(10, -12);
var ang = math.pow(10, -10);
var nm = math.pow(10, -9);
var um = math.pow(10, -6);
var mm = math.pow(10, -3);
var mrad = math.pow(10, -3);
var PI = math.pi;
var deg = PI / 180;
var correction_factor = 1;

var canvas1 = document.getElementById("canvas1");
var ctx1 = canvas1.getContext("2d");
var canvas2 = document.getElementById("canvas2");
var ctx2 = canvas2.getContext("2d");
var canvas3 = document.getElementById("canvas3");
var ctx3 = canvas3.getContext("2d");

var forceJS = document.getElementById("forceJS"); //figure out how to read from checkbox
var interactiveMode = document.getElementById("interactiveMode");
var real = document.getElementById("realMode");
var aberration_list = [
    "C10",
    "C12",
    "C21",
    "C23",
    "C30",
    "C32",
    "C34",
    "C41",
    "C43",
    "C45",
    "C50",
    "C52",
    "C54",
    "C56",
];

var number_aberration_terms = aberration_list.length;
var aberrations = [];

// var random_dummy_aberration_coefs = [ [1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [0.000000, 0.000000, 2.000000, 150.000000, 0.500000, 10.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [2.000000, 0.000000, 1.500000, 120.000000, 0.000000, 0.000000, 0.300000, 35.000000, 1.100000, 0.000000, 0.500000, 20.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [5.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 3.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [0.000000, 0.000000, 0.200000, 10.000000, 0.100000, 60.000000, 0.300000, 10.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [1.000000, 0.000000, 0.500000, 0.000000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [8.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [0.000000, 0.000000, 1.000000, 0.000000, 2.000000, -70.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] , 
// [0.500000, 0.000000, 0.250000, -25.000000, 0.600000, 25.000000, 0.220000, -15.000000, 0.100000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.100000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000] ];

var real_aberration_coefs = [ [-2.200000, 0.000000, 1.013000, -73.060000, 19.124000, -174.090000, 19.124000, -44.330000, 0.199000, 0.000000, 0.188000, 4.040000, 0.055000, -28.080000, 1.526000, -7.290000, 1.408000, 12.360000, 0.165000, -25.200000, 0.005000, 0.000000, 0.018000, -7.510000, 0.000000, 0.000000, 0.000000, 0.000000] ,
[2.648000, 0.000000, 0.739000, -15.990000, 38.458000, 50.850000, 76.944000, 39.820000, 7.920000, 0.000000, 5.151000, 1.330000, 0.428000, -7.870000, 60.289000, -117.570000, 77.468000, -6.700000, 26.582000, 29.660000, -3.324000, 0.000000, 4.739000, -24.460000, 0.000000, 0.000000, 0.000000, 0.000000] ,
 [0.324000, 0.000000, 3.436000, -82.280000, 65.349000, 60.680000, 81.814000, 21.840000, 8.574000, 0.000000, 7.875000, 3.170000, 0.621000, -5.510000, 76.514000, -99.710000, 24.948000, -27.110000, 25.185000, 27.670000, -2.816000, 0.000000, 4.751000, -24.060000, 0.000000, 0.000000, 0.000000, 0.000000] , 
 [-17.421000, 0.000000, 2.084000, -12.520000, 9.878000, 29.980000, 50.932000, 28.830000, -3.057000, 0.000000, 3.478000, 76.050000, 1.136000, -2.300000, 52.522000, 5.050000, 43.404000, -12.170000, 4.306000, 30.970000, 1.876000, 0.000000, 0.439000, -0.320000, 0.000000, 0.000000, 0.000000, 0.000000] , 
 [-21.172000, 0.000000, 3.241000, -63.690000, 61.031000, -85.350000, 38.950000, 16.040000, 3.248000, 0.000000, 0.379000, 67.230000, 0.217000, 24.770000, 7.211000, 57.740000, 7.318000, -45.860000, 1.015000, 21.810000, -0.202000, 0.000000, 0.071000, -5.050000, 0.000000, 0.000000, 0.000000, 0.000000] ,  
 [0.145000, 0.000000, 3.229000, 77.490000, 63.263000, -129.920000, 20.832000, 57.440000, 9.248000, 0.000000, 7.526000, -3.410000, 1.036000, -3.250000, 45.978000, -43.440000, 65.089000, 11.660000, 27.597000, 33.250000, -3.620000, 0.000000, 5.164000, -24.780000, 0.000000, 0.000000, 0.000000, 0.000000] , 
 [-0.452000, 0.000000, 3.578000, 84.160000, 209.004000, 160.040000, 128.827000, 43.550000, 12.420000, 0.000000, 6.666000, -1.830000, 0.500000, -7.710000, 192.217000, -44.520000, 116.666000, -5.590000, 13.357000, 23.120000, -5.239000, 0.000000, 5.025000, -24.290000, 0.000000, 0.000000, 0.000000, 0.000000] , 
 [0.309000, 0.000000, 5.365000, -81.810000, 1524.835000, -163.820000, 80.842000, -56.660000, 5.275000, 0.000000, 2.925000, 4.430000, 1.426000, -2.970000, 616.565000, 14.170000, 89.681000, 2.090000, 66.762000, -33.050000, -3.171000, 0.000000, 3.061000, -25.290000, 0.000000, 0.000000, 0.000000, 0.000000] ,
 [2.027000, 0.000000, 1.000000, -74.220000, 18.867000, 43.920000, 90.089000, -21.610000, 11.680000, 0.000000, 5.229000, -4.010000, 0.345000, 36.600000, 71.049000, -123.790000, 70.408000, 27.530000, 24.894000, 32.140000, -5.988000, 0.000000, 4.805000, -24.110000, 0.000000, 0.000000, 0.000000, 0.000000] ];

var random_reasonable_aberration_coefs = [ [2.860000, 0.000000, -7.520000, 258.230000, -31.900000, 218.510000, -14.550000, 87.440000, -1.490000, 0.000000, -1.610000, 144.340000, 1.100000, 39.180000, 0.010000, 78.580000, 0.030000, 302.090000, 0.030000, 106.580000, -0.020000, 0.000000, -0.870000, 350.200000, 0.120000, 277.690000, -0.410000, 277.170000] ,
[-4.550000, 0.000000, 6.470000, 126.890000, -28.370000, 330.850000, 1.050000, 341.760000, -1.340000, 0.000000, 1.710000, 238.770000, -0.950000, 125.610000, -0.000000, 7.210000, -0.050000, 22.720000, 0.090000, 349.430000, -0.980000, 0.000000, -0.520000, 194.320000, 0.060000, 273.690000, 0.140000, 240.390000] ,
[-1.300000, 0.000000, 9.050000, 335.470000, -14.220000, 259.550000, 5.620000, 265.880000, 0.830000, 0.000000, -0.260000, 59.760000, 0.730000, 316.830000, -0.010000, 118.930000, -0.030000, 321.620000, 0.100000, 247.210000, 0.670000, 0.000000, -0.250000, 309.130000, -0.420000, 332.630000, 0.780000, 293.320000] ,
[-0.360000, 0.000000, 9.460000, 77.700000, -23.230000, 331.250000, 5.670000, 317.190000, 0.530000, 0.000000, 1.250000, 101.190000, -0.110000, 110.690000, -0.010000, 81.400000, 0.020000, 99.450000, -0.060000, 149.950000, -1.660000, 0.000000, -0.010000, 45.390000, 2.040000, 273.780000, 0.780000, 336.610000] ,
[5.250000, 0.000000, -4.320000, 132.720000, 5.580000, 83.620000, -13.940000, 87.990000, -1.630000, 0.000000, -1.460000, 285.650000, -0.930000, 268.230000, -0.010000, 342.040000, -0.060000, 187.770000, 0.070000, 86.430000, 0.660000, 0.000000, 0.590000, 348.270000, -1.700000, 273.510000, 0.090000, 48.570000] ,
[-9.030000, 0.000000, -0.810000, 73.680000, 4.840000, 295.090000, -17.970000, 272.010000, 2.190000, 0.000000, 1.710000, 73.560000, 1.100000, 45.170000, 0.020000, 19.470000, -0.100000, 26.040000, 0.020000, 332.310000, -1.410000, 0.000000, 1.740000, 141.010000, -0.590000, 295.100000, 0.340000, 198.900000] ,
[3.940000, 0.000000, 0.650000, 35.880000, -12.900000, 272.630000, 3.090000, 357.210000, 1.090000, 0.000000, -2.020000, 226.410000, 0.740000, 269.210000, 0.020000, 333.140000, 0.110000, 299.180000, 0.090000, 267.780000, 0.710000, 0.000000, -1.410000, 179.020000, 1.640000, 298.810000, 0.630000, 27.720000] ,
[2.720000, 0.000000, 4.220000, 82.490000, -11.290000, 114.080000, -17.080000, 83.320000, -1.200000, 0.000000, 0.050000, 234.410000, -0.480000, 349.730000, 0.010000, 196.600000, -0.000000, 40.790000, 0.100000, 213.320000, -0.680000, 0.000000, -0.270000, 305.170000, -0.650000, 1.170000, 1.400000, 215.460000] ,
[3.690000, 0.000000, -0.370000, 173.870000, 13.990000, 109.790000, 4.890000, 65.730000, -0.370000, 0.000000, 0.760000, 250.560000, -0.330000, 229.560000, 0.010000, 66.470000, 0.050000, 225.780000, 0.060000, 118.220000, 1.770000, 0.000000, 0.650000, 246.520000, 0.140000, 92.620000, -1.000000, 31.560000] ,
[3.910000, 0.000000, -8.160000, 33.750000, 5.060000, 130.180000, 6.690000, 213.560000, 1.210000, 0.000000, -0.740000, 103.820000, 1.070000, 68.320000, 0.020000, 1.290000, -0.060000, 119.340000, 0.100000, 157.140000, 1.120000, 0.000000, -1.590000, 251.670000, -0.480000, 246.890000, 1.870000, 278.740000] ,
[7.610000, 0.000000, 6.170000, 73.280000, -13.360000, 197.300000, 16.440000, 325.780000, -0.000000, 0.000000, -1.470000, 207.440000, 0.800000, 98.610000, -0.000000, 177.270000, -0.000000, 305.620000, -0.060000, 104.780000, 0.960000, 0.000000, 0.440000, 50.070000, 1.430000, 177.280000, -1.350000, 260.740000] ,
[-0.030000, 0.000000, -7.590000, 43.660000, -11.550000, 129.760000, 16.390000, 335.490000, 1.480000, 0.000000, -0.010000, 294.530000, -0.230000, 120.590000, 0.010000, 237.180000, -0.070000, 93.210000, -0.080000, 26.120000, -0.570000, 0.000000, -0.700000, 103.780000, -0.300000, 32.820000, 0.360000, 336.420000] ,
[3.340000, 0.000000, -0.260000, 274.250000, 25.270000, 56.620000, 0.720000, 225.240000, 0.260000, 0.000000, 1.450000, 153.440000, -0.560000, 141.980000, 0.010000, 117.370000, 0.110000, 229.920000, 0.090000, 121.770000, -0.370000, 0.000000, 1.200000, 1.950000, -0.860000, 278.780000, 1.540000, 41.290000] ,
[-9.460000, 0.000000, 10.220000, 161.680000, -19.060000, 254.850000, 14.660000, 170.610000, -1.740000, 0.000000, -0.860000, 137.850000, 0.680000, 236.570000, -0.010000, 47.420000, -0.000000, 19.240000, 0.050000, 281.120000, -1.600000, 0.000000, 0.340000, 212.270000, 0.400000, 190.770000, -0.820000, 130.300000] ,
[-0.490000, 0.000000, 2.310000, 61.140000, 7.850000, 189.270000, -10.680000, 214.640000, -1.870000, 0.000000, 1.860000, 35.590000, -0.040000, 61.080000, 0.020000, 81.180000, -0.020000, 104.700000, -0.030000, 316.180000, 0.670000, 0.000000, -1.010000, 13.080000, 0.530000, 280.180000, -0.800000, 301.000000] ,
[-6.340000, 0.000000, -8.190000, 220.480000, 18.620000, 242.860000, -12.020000, 259.010000, -0.800000, 0.000000, -1.170000, 156.250000, 0.080000, 138.870000, 0.010000, 55.710000, -0.020000, 5.250000, -0.040000, 137.590000, -1.000000, 0.000000, 0.220000, 233.880000, 0.790000, 331.060000, 0.840000, 291.530000] ,
[3.070000, 0.000000, 0.700000, 2.170000, 7.820000, 303.810000, 0.750000, 231.370000, -0.600000, 0.000000, 1.320000, 258.800000, -0.760000, 244.020000, -0.010000, 11.850000, -0.060000, 246.860000, 0.050000, 222.830000, -2.090000, 0.000000, -0.810000, 2.060000, 0.660000, 94.170000, -1.330000, 308.720000] ,
[3.520000, 0.000000, 3.230000, 316.450000, 25.400000, 112.770000, -13.750000, 67.060000, 1.440000, 0.000000, 1.770000, 243.240000, -0.230000, 68.810000, 0.020000, 254.190000, 0.060000, 197.070000, -0.050000, 335.700000, 0.220000, 0.000000, -0.020000, 336.040000, 1.850000, 198.930000, 1.320000, 287.880000] ,
[3.310000, 0.000000, 9.160000, 358.310000, 24.690000, 116.840000, 5.530000, 212.100000, 1.210000, 0.000000, -1.030000, 286.170000, -0.060000, 217.580000, 0.020000, 60.110000, 0.080000, 311.440000, -0.010000, 239.190000, 0.410000, 0.000000, 0.170000, 232.420000, 0.340000, 53.410000, 0.850000, 11.870000] ,
[6.990000, 0.000000, -8.130000, 185.420000, 0.690000, 176.340000, 12.620000, 17.460000, 0.610000, 0.000000, -1.550000, 162.770000, -0.550000, 148.710000, -0.020000, 146.440000, 0.020000, 258.340000, 0.020000, 292.670000, -0.090000, 0.000000, -1.850000, 358.260000, 0.600000, 26.740000, -1.160000, 215.030000] ,
[2.740000, 0.000000, 5.000000, 332.470000, -4.050000, 166.630000, 18.170000, 306.220000, 1.750000, 0.000000, -0.720000, 276.130000, -0.610000, 193.230000, 0.020000, 171.920000, 0.090000, 167.830000, -0.060000, 348.220000, 1.180000, 0.000000, 1.710000, 275.930000, 1.100000, 92.740000, -0.700000, 346.870000] ,
[1.280000, 0.000000, 2.570000, 199.610000, -1.450000, 68.770000, 6.180000, 129.640000, -1.260000, 0.000000, 1.600000, 218.360000, -0.270000, 39.520000, 0.010000, 71.650000, 0.040000, 213.370000, -0.090000, 214.690000, 0.270000, 0.000000, -2.020000, 87.350000, -2.060000, 123.790000, 0.430000, 332.530000] ,
[8.140000, 0.000000, -6.860000, 335.980000, -0.800000, 161.280000, 5.580000, 286.290000, -1.500000, 0.000000, 1.670000, 105.450000, 0.880000, 131.780000, -0.000000, 269.150000, 0.100000, 98.280000, 0.080000, 44.040000, 0.920000, 0.000000, -1.320000, 332.840000, -1.390000, 101.630000, 0.530000, 73.080000] ,
[-7.840000, 0.000000, 9.390000, 81.920000, -22.370000, 5.000000, -1.530000, 43.200000, 1.820000, 0.000000, 0.500000, 36.310000, -0.230000, 25.410000, -0.000000, 178.720000, -0.050000, 105.550000, 0.020000, 328.470000, -1.950000, 0.000000, 1.320000, 155.470000, -0.600000, 271.220000, -1.950000, 359.270000] ,
[-6.280000, 0.000000, 4.200000, 238.250000, 25.680000, 117.950000, -6.350000, 232.820000, 1.170000, 0.000000, 0.870000, 289.200000, 0.890000, 245.500000, 0.020000, 112.660000, 0.070000, 107.260000, 0.020000, 68.070000, -1.670000, 0.000000, -1.070000, 56.620000, 0.380000, 49.030000, 1.640000, 20.900000] ,
[-9.320000, 0.000000, -0.650000, 333.080000, 5.740000, 92.510000, 3.400000, 60.790000, 1.380000, 0.000000, 0.120000, 333.390000, 0.510000, 209.620000, -0.000000, 81.090000, 0.010000, 228.100000, 0.100000, 6.000000, -1.230000, 0.000000, -0.410000, 188.150000, 0.760000, 110.590000, -0.230000, 232.250000] ,
[4.270000, 0.000000, -5.990000, 119.850000, -15.950000, 273.320000, -19.410000, 246.090000, 1.540000, 0.000000, -0.770000, 216.280000, 0.060000, 240.470000, -0.010000, 305.290000, -0.080000, 92.250000, 0.090000, 185.180000, 0.140000, 0.000000, 1.930000, 295.680000, -0.650000, 265.230000, -1.910000, 129.580000] ,
[-10.390000, 0.000000, -4.350000, 175.420000, 21.120000, 254.980000, -1.300000, 182.670000, -1.350000, 0.000000, 1.850000, 174.120000, -0.910000, 15.630000, 0.010000, 88.150000, -0.080000, 220.050000, 0.080000, 346.170000, -0.100000, 0.000000, 1.170000, 272.630000, 0.340000, 2.520000, 1.030000, 265.130000] ,
[-8.470000, 0.000000, 9.380000, 283.560000, -14.840000, 36.540000, 12.430000, 86.160000, 1.080000, 0.000000, -1.420000, 99.800000, 0.580000, 337.580000, 0.030000, 34.810000, -0.030000, 304.300000, -0.000000, 249.290000, 0.650000, 0.000000, 0.430000, 116.640000, -1.750000, 46.800000, -1.520000, 136.080000] ,
[-8.090000, 0.000000, 1.740000, 317.050000, 11.110000, 75.920000, -7.520000, 190.400000, 1.180000, 0.000000, -2.010000, 43.960000, -0.220000, 185.380000, -0.000000, 76.170000, -0.040000, 57.660000, -0.100000, 156.160000, -1.570000, 0.000000, -1.740000, 166.310000, 1.210000, 280.900000, 1.730000, 246.490000] ,
[-9.180000, 0.000000, 0.570000, 16.990000, 28.220000, 63.770000, -4.490000, 39.430000, 0.610000, 0.000000, 0.180000, 252.130000, 0.710000, 146.390000, 0.020000, 207.970000, 0.070000, 79.830000, -0.020000, 134.720000, 1.300000, 0.000000, 1.040000, 139.440000, 1.470000, 336.310000, 0.910000, 299.330000] ,
[0.340000, 0.000000, 0.060000, 224.880000, 11.340000, 208.380000, 10.270000, 10.620000, -1.570000, 0.000000, -1.320000, 51.820000, 0.760000, 279.870000, -0.010000, 264.720000, -0.030000, 240.020000, 0.090000, 38.010000, -0.080000, 0.000000, -1.260000, 97.320000, 0.670000, 103.590000, -1.160000, 340.930000] ,
[5.850000, 0.000000, -7.510000, 337.090000, 4.030000, 106.060000, 14.990000, 232.030000, 0.770000, 0.000000, 1.460000, 226.570000, 0.600000, 292.680000, 0.010000, 10.250000, 0.060000, 222.330000, -0.040000, 278.800000, 1.600000, 0.000000, -1.340000, 62.180000, -0.230000, 129.530000, 0.620000, 136.340000] ,
[-3.660000, 0.000000, 2.270000, 313.000000, 20.170000, 37.510000, -4.060000, 269.900000, -0.460000, 0.000000, 0.440000, 97.990000, 0.890000, 24.570000, 0.020000, 155.600000, -0.000000, 242.950000, -0.040000, 169.790000, -1.060000, 0.000000, -0.830000, 207.620000, 1.050000, 45.390000, -0.050000, 172.860000] ,
[7.540000, 0.000000, -4.340000, 48.020000, -32.790000, 66.450000, -8.470000, 324.280000, 0.660000, 0.000000, -1.410000, 247.510000, 1.010000, 16.110000, -0.010000, 279.030000, -0.060000, 310.060000, 0.010000, 357.920000, -0.140000, 0.000000, 1.650000, 298.280000, -0.390000, 256.290000, 1.420000, 69.660000] ,
[3.650000, 0.000000, -3.200000, 116.590000, -2.820000, 191.710000, 8.710000, 230.540000, 0.180000, 0.000000, -1.010000, 2.110000, 0.430000, 68.990000, -0.010000, 159.880000, -0.050000, 235.120000, 0.010000, 125.760000, -0.650000, 0.000000, -1.700000, 18.170000, 1.600000, 291.370000, 0.090000, 355.200000] ,
[-6.710000, 0.000000, -1.280000, 129.330000, 12.080000, 307.370000, -4.560000, 283.030000, -0.320000, 0.000000, 0.440000, 37.230000, -0.850000, 348.370000, -0.010000, 312.880000, -0.030000, 288.940000, -0.000000, 142.450000, 1.030000, 0.000000, 1.930000, 349.220000, -2.010000, 31.880000, 1.680000, 19.110000] ,
[-9.660000, 0.000000, -6.640000, 150.900000, 18.370000, 79.150000, -17.100000, 224.210000, -0.390000, 0.000000, 0.450000, 165.580000, 0.140000, 300.800000, 0.020000, 72.810000, -0.090000, 242.130000, -0.000000, 346.520000, -0.590000, 0.000000, 1.760000, 335.110000, 1.690000, 37.170000, -1.290000, 315.220000] ,
[9.010000, 0.000000, 6.750000, 32.280000, 0.590000, 348.630000, -12.500000, 1.710000, 0.560000, 0.000000, 1.090000, 94.280000, 0.490000, 13.140000, 0.020000, 126.190000, 0.020000, 102.660000, 0.100000, 243.100000, 0.580000, 0.000000, -2.060000, 155.840000, -0.570000, 358.580000, -1.680000, 333.160000] ,
[6.330000, 0.000000, -6.590000, 104.790000, -7.860000, 262.700000, 7.950000, 2.510000, -2.130000, 0.000000, -1.410000, 319.090000, -0.770000, 230.120000, 0.010000, 333.140000, -0.040000, 313.220000, -0.070000, 197.190000, -1.160000, 0.000000, -0.190000, 40.650000, 0.200000, 309.890000, 1.500000, 166.050000] ,
[-2.380000, 0.000000, -6.840000, 344.200000, -16.670000, 67.570000, 16.750000, 129.660000, -0.580000, 0.000000, 0.080000, 221.580000, 0.010000, 136.370000, -0.020000, 249.880000, 0.030000, 234.030000, 0.010000, 249.730000, -1.080000, 0.000000, 0.950000, 243.120000, 1.220000, 167.190000, -0.730000, 63.600000] ,
[2.900000, 0.000000, -5.380000, 65.530000, -23.860000, 8.850000, -11.540000, 150.360000, 1.220000, 0.000000, -2.120000, 46.680000, 0.540000, 52.190000, 0.010000, 191.000000, -0.050000, 88.920000, 0.090000, 195.030000, -1.110000, 0.000000, 1.010000, 150.190000, -1.710000, 171.670000, 1.860000, 166.840000] ,
[2.950000, 0.000000, -5.990000, 95.920000, -3.400000, 100.090000, -7.920000, 164.980000, 1.650000, 0.000000, 2.070000, 185.350000, 0.320000, 235.360000, 0.020000, 354.600000, 0.010000, 140.280000, 0.030000, 116.090000, 1.530000, 0.000000, 0.360000, 192.140000, 1.900000, 35.870000, 0.330000, 116.560000] ,
[6.040000, 0.000000, 8.910000, 186.390000, 18.750000, 153.830000, 8.240000, 49.990000, 0.430000, 0.000000, 1.400000, 2.840000, 0.050000, 170.300000, -0.010000, 284.740000, 0.060000, 60.210000, 0.020000, 333.190000, -1.050000, 0.000000, -1.960000, 69.690000, 0.960000, 40.380000, -1.210000, 221.730000] ,
[9.240000, 0.000000, 4.660000, 283.780000, -12.790000, 261.670000, -13.870000, 232.440000, 1.250000, 0.000000, -0.770000, 309.400000, 0.810000, 137.380000, 0.020000, 51.050000, 0.100000, 140.420000, 0.020000, 353.330000, 1.710000, 0.000000, -2.100000, 41.120000, -0.800000, 55.780000, -0.970000, 11.400000] ,
[3.100000, 0.000000, 4.480000, 172.400000, -15.350000, 211.400000, -7.460000, 156.510000, -1.830000, 0.000000, 0.720000, 168.550000, -0.370000, 311.360000, -0.010000, 354.950000, -0.070000, 70.190000, 0.020000, 36.410000, -1.750000, 0.000000, 1.140000, 51.230000, -1.020000, 356.240000, -1.490000, 276.890000] ,
[6.320000, 0.000000, 2.090000, 148.160000, 25.760000, 161.400000, -18.630000, 112.500000, 1.090000, 0.000000, 1.010000, 125.850000, 0.700000, 297.950000, 0.010000, 141.820000, -0.020000, 52.340000, -0.030000, 338.030000, -1.930000, 0.000000, 0.290000, 333.120000, -1.810000, 44.160000, -1.620000, 121.380000] ,
[-8.270000, 0.000000, 8.160000, 92.410000, 11.100000, 326.540000, -5.800000, 175.550000, 1.320000, 0.000000, 0.950000, 140.720000, 0.530000, 197.050000, -0.000000, 160.880000, 0.080000, 200.960000, 0.050000, 24.290000, -2.070000, 0.000000, 0.960000, 15.810000, 1.980000, 16.750000, -1.380000, 106.700000] ,
[2.810000, 0.000000, -2.920000, 101.300000, 2.510000, 266.330000, 5.900000, 89.740000, 1.040000, 0.000000, 1.160000, 0.800000, -0.230000, 193.340000, -0.010000, 173.210000, 0.090000, 48.080000, 0.050000, 165.380000, 0.120000, 0.000000, 2.060000, 142.320000, -0.250000, 191.980000, 0.340000, 258.410000] ,
[-10.200000, 0.000000, -1.150000, 269.620000, -30.970000, 233.170000, -4.510000, 290.220000, -1.940000, 0.000000, 1.200000, 12.360000, -1.010000, 285.240000, -0.000000, 194.260000, -0.100000, 154.670000, 0.020000, 148.950000, 1.570000, 0.000000, 0.920000, 91.770000, 0.640000, 143.980000, 1.820000, 254.520000] ,
[-1.580000, 0.000000, -0.040000, 184.660000, -15.700000, 158.020000, -16.580000, 339.510000, 1.090000, 0.000000, 1.580000, 170.670000, 0.020000, 289.680000, 0.010000, 228.690000, -0.060000, 196.160000, 0.100000, 333.470000, -0.670000, 0.000000, -1.490000, 193.560000, 1.230000, 76.730000, 0.480000, 310.240000] ,
[1.450000, 0.000000, 0.680000, 196.640000, -29.070000, 357.670000, -2.700000, 10.830000, 2.070000, 0.000000, 0.500000, 186.860000, 1.090000, 260.060000, 0.010000, 170.590000, -0.080000, 180.540000, -0.090000, 44.630000, -1.880000, 0.000000, -0.080000, 7.520000, 0.250000, 140.510000, 0.440000, 224.540000] ,
[-0.300000, 0.000000, 6.180000, 77.930000, 17.330000, 19.520000, -18.140000, 117.400000, -0.620000, 0.000000, -1.560000, 172.870000, 1.100000, 118.560000, -0.010000, 272.220000, -0.050000, 357.250000, -0.000000, 208.340000, 0.860000, 0.000000, 0.520000, 194.480000, 0.260000, 271.060000, -1.490000, 36.880000] ,
[4.760000, 0.000000, -0.610000, 268.790000, 30.550000, 76.180000, -9.390000, 95.240000, 0.400000, 0.000000, 0.210000, 112.890000, 0.650000, 340.860000, -0.000000, 248.590000, -0.010000, 285.480000, 0.060000, 118.010000, 1.710000, 0.000000, -2.030000, 100.510000, 1.350000, 219.310000, -0.570000, 229.930000] ,
[5.070000, 0.000000, 0.010000, 181.930000, -16.350000, 168.220000, 7.180000, 349.450000, -1.150000, 0.000000, 1.930000, 340.100000, -0.530000, 297.330000, 0.020000, 175.850000, -0.060000, 18.730000, 0.080000, 17.580000, -0.510000, 0.000000, -0.830000, 226.930000, 1.730000, 101.930000, 1.650000, 114.460000] ,
[9.600000, 0.000000, -8.200000, 91.600000, 17.140000, 251.190000, -13.430000, 219.370000, 0.350000, 0.000000, -1.190000, 303.810000, 0.640000, 293.510000, 0.010000, 60.360000, 0.040000, 16.870000, -0.080000, 149.010000, -0.150000, 0.000000, -0.160000, 339.800000, 0.320000, 302.170000, -1.580000, 274.640000] ,
[-8.800000, 0.000000, -6.030000, 5.340000, -15.210000, 13.470000, 9.440000, 115.920000, 1.890000, 0.000000, -2.000000, 324.060000, 0.540000, 184.100000, -0.000000, 96.330000, -0.010000, 191.760000, 0.020000, 348.440000, -0.490000, 0.000000, 1.800000, 244.210000, 0.520000, 304.020000, 0.780000, 139.330000] ,
[-7.350000, 0.000000, 3.790000, 290.280000, 21.180000, 83.050000, -9.300000, 321.420000, -1.240000, 0.000000, -1.660000, 1.170000, 0.890000, 52.610000, 0.010000, 57.650000, -0.010000, 233.090000, -0.060000, 354.320000, 1.900000, 0.000000, 1.910000, 273.460000, -0.260000, 120.180000, -1.350000, 317.020000] ,
[-4.880000, 0.000000, -5.110000, 310.380000, 11.120000, 150.300000, 5.950000, 337.240000, 1.440000, 0.000000, 0.350000, 329.860000, 0.920000, 15.980000, -0.010000, 251.770000, 0.080000, 296.240000, -0.000000, 305.290000, 1.400000, 0.000000, -1.460000, 234.160000, -0.370000, 282.350000, -1.600000, 37.560000] ,
[10.190000, 0.000000, 8.580000, 106.720000, -22.510000, 90.920000, 6.830000, 59.290000, -1.740000, 0.000000, -0.800000, 114.550000, -0.580000, 335.170000, 0.020000, 322.790000, 0.030000, 55.420000, 0.110000, 59.720000, -2.040000, 0.000000, 1.020000, 160.970000, -1.690000, 304.610000, 0.550000, 321.380000] ,
[10.450000, 0.000000, -3.030000, 269.480000, -30.000000, 353.520000, -10.650000, 94.910000, -0.310000, 0.000000, -2.050000, 325.640000, -0.600000, 190.720000, -0.000000, 125.030000, -0.050000, 164.590000, 0.030000, 157.020000, 0.430000, 0.000000, -0.850000, 211.650000, -1.950000, 221.740000, -0.150000, 13.660000] ,
[9.050000, 0.000000, 2.950000, 33.300000, 22.500000, 333.980000, 16.440000, 358.430000, -1.050000, 0.000000, -1.790000, 51.920000, -0.980000, 247.520000, -0.010000, 43.290000, -0.000000, 97.730000, 0.050000, 244.110000, -1.620000, 0.000000, 0.220000, 122.310000, 1.800000, 257.850000, -0.700000, 304.790000] ,
[8.080000, 0.000000, -8.500000, 285.930000, 25.510000, 293.480000, 9.750000, 264.720000, 1.010000, 0.000000, 0.500000, 234.940000, 0.660000, 356.460000, 0.010000, 253.230000, 0.070000, 304.390000, -0.080000, 322.210000, 0.310000, 0.000000, -2.080000, 304.430000, -1.320000, 247.430000, -0.400000, 199.420000] ,
[-2.080000, 0.000000, -8.910000, 264.870000, 1.980000, 102.630000, -16.030000, 61.500000, 1.780000, 0.000000, -1.550000, 303.550000, -0.000000, 229.340000, 0.010000, 272.990000, -0.000000, 105.320000, -0.070000, 109.780000, -1.260000, 0.000000, 0.330000, 91.830000, -1.670000, 278.290000, -0.160000, 208.510000] ,
[-7.700000, 0.000000, -2.120000, 310.990000, 6.570000, 192.110000, 12.790000, 170.890000, 0.640000, 0.000000, 0.210000, 330.840000, 0.580000, 198.590000, -0.010000, 248.670000, -0.010000, 92.860000, 0.020000, 293.730000, -1.590000, 0.000000, -1.270000, 247.890000, 1.870000, 115.480000, -1.700000, 278.160000] ,
[-3.110000, 0.000000, -3.560000, 199.480000, -5.410000, 174.580000, -19.310000, 262.620000, -1.290000, 0.000000, -1.540000, 301.420000, 0.570000, 305.890000, -0.000000, 249.600000, -0.020000, 186.640000, 0.060000, 212.410000, -0.390000, 0.000000, 1.030000, 117.680000, -2.030000, 189.640000, 1.460000, 338.360000] ,
[4.470000, 0.000000, -10.420000, 338.880000, -0.300000, 23.170000, 1.980000, 120.210000, -1.900000, 0.000000, 1.890000, 204.420000, -0.200000, 96.330000, -0.010000, 28.060000, 0.060000, 58.090000, 0.040000, 186.700000, -1.640000, 0.000000, 2.060000, 164.220000, 1.190000, 185.230000, -1.920000, 263.570000] ,
[3.620000, 0.000000, 5.420000, 319.970000, 21.780000, 138.440000, -2.050000, 273.980000, -1.760000, 0.000000, -1.140000, 359.340000, 0.150000, 57.000000, -0.000000, 61.110000, -0.050000, 349.580000, -0.070000, 195.540000, -1.840000, 0.000000, -1.360000, 292.390000, -0.970000, 61.340000, 0.780000, 58.530000] ,
[8.300000, 0.000000, 7.350000, 261.430000, 7.590000, 204.030000, 18.110000, 218.680000, -0.570000, 0.000000, 1.190000, 144.060000, -0.220000, 166.820000, 0.010000, 3.410000, -0.060000, 347.130000, 0.100000, 17.140000, 0.390000, 0.000000, -0.220000, 24.380000, 1.600000, 234.910000, 1.360000, 223.230000] ,
[5.940000, 0.000000, -6.300000, 182.800000, -17.540000, 243.850000, -10.260000, 17.820000, 0.690000, 0.000000, -0.900000, 70.330000, 0.210000, 9.130000, -0.000000, 23.260000, -0.080000, 357.760000, 0.100000, 39.760000, -1.430000, 0.000000, 0.760000, 321.600000, 1.940000, 269.560000, -0.410000, 45.370000] ,
[5.170000, 0.000000, 10.280000, 80.700000, 15.300000, 190.340000, -11.790000, 67.960000, -1.140000, 0.000000, 1.380000, 162.050000, 0.320000, 322.250000, 0.020000, 39.170000, -0.060000, 86.740000, -0.050000, 147.320000, 0.090000, 0.000000, 0.010000, 69.390000, -1.730000, 243.960000, -0.300000, 65.120000] ,
[-4.050000, 0.000000, 8.150000, 297.980000, -29.550000, 19.080000, 3.280000, 315.760000, -1.890000, 0.000000, 1.100000, 283.590000, -0.570000, 109.470000, -0.010000, 203.220000, 0.040000, 317.990000, -0.070000, 43.090000, 0.130000, 0.000000, -1.700000, 129.230000, 0.220000, 17.480000, 0.960000, 214.130000] ,
[5.800000, 0.000000, 3.890000, 54.940000, 31.680000, 29.730000, -14.620000, 206.200000, -0.210000, 0.000000, 1.390000, 258.210000, 0.020000, 185.660000, 0.020000, 203.650000, -0.070000, 267.410000, 0.030000, 6.840000, -1.510000, 0.000000, 0.640000, 270.130000, -1.640000, 240.330000, 0.910000, 267.900000] ,
[-3.380000, 0.000000, -4.220000, 159.230000, 6.260000, 41.260000, -12.140000, 354.060000, 0.250000, 0.000000, 0.470000, 119.780000, -0.980000, 2.210000, 0.010000, 152.420000, 0.110000, 199.770000, -0.040000, 122.990000, -0.580000, 0.000000, 1.780000, 332.570000, 0.740000, 179.740000, -1.390000, 206.190000] ,
[-3.840000, 0.000000, -1.160000, 317.860000, -11.520000, 236.100000, 10.850000, 267.510000, 0.990000, 0.000000, -0.020000, 348.330000, -0.430000, 101.830000, -0.000000, 36.970000, 0.010000, 126.080000, 0.090000, 291.890000, -1.460000, 0.000000, 0.570000, 73.010000, -1.560000, 184.690000, -2.040000, 200.530000] ,
[-7.770000, 0.000000, 1.770000, 64.490000, -28.880000, 160.790000, -15.940000, 10.270000, 1.190000, 0.000000, -0.840000, 314.230000, 0.750000, 179.340000, 0.020000, 288.960000, 0.090000, 50.650000, -0.000000, 25.180000, 1.600000, 0.000000, 1.880000, 142.560000, 0.420000, 11.940000, 0.200000, 209.120000] ,
[-7.600000, 0.000000, 0.550000, 201.060000, 15.680000, 96.350000, 8.600000, 39.670000, -1.570000, 0.000000, -1.360000, 294.010000, 0.420000, 327.700000, -0.000000, 177.960000, -0.040000, 190.780000, 0.040000, 68.190000, -1.010000, 0.000000, 0.570000, 58.140000, 0.250000, 50.430000, -1.370000, 208.900000] ,
[-7.110000, 0.000000, 8.020000, 258.990000, -14.650000, 107.650000, 2.680000, 146.770000, 0.080000, 0.000000, 1.380000, 101.110000, -0.880000, 236.360000, 0.010000, 0.770000, -0.000000, 283.330000, -0.030000, 351.080000, -1.400000, 0.000000, 2.090000, 355.520000, -1.470000, 152.630000, 0.270000, 226.610000] ,
[-6.090000, 0.000000, 7.630000, 265.420000, -2.850000, 133.520000, 6.820000, 267.990000, -1.520000, 0.000000, 1.100000, 85.140000, 0.040000, 241.180000, -0.000000, 203.470000, 0.020000, 221.350000, -0.020000, 321.580000, 1.550000, 0.000000, 1.510000, 258.780000, 0.870000, 11.920000, -0.180000, 307.040000] ,
[-0.370000, 0.000000, 1.310000, 7.600000, 17.050000, 249.150000, -17.530000, 152.900000, -1.440000, 0.000000, 1.970000, 263.820000, 1.030000, 115.630000, 0.030000, 251.860000, -0.040000, 175.200000, -0.080000, 114.320000, -1.210000, 0.000000, 1.670000, 181.230000, -1.160000, 28.550000, -1.620000, 271.070000] ,
[2.230000, 0.000000, 5.970000, 205.000000, 5.980000, 31.760000, 11.250000, 124.550000, 1.190000, 0.000000, 0.280000, 303.570000, 0.170000, 338.330000, -0.010000, 184.860000, -0.060000, 196.330000, -0.090000, 90.740000, 0.300000, 0.000000, -0.790000, 48.070000, 0.580000, 280.950000, 1.520000, 74.100000] ,
[9.640000, 0.000000, 6.050000, 349.950000, 2.740000, 203.270000, 6.140000, 205.840000, 1.830000, 0.000000, -0.080000, 155.930000, -0.490000, 247.830000, 0.010000, 16.130000, 0.020000, 306.870000, -0.020000, 320.100000, -1.500000, 0.000000, -1.050000, 154.720000, 0.270000, 255.900000, -0.030000, 201.330000] ,
[5.570000, 0.000000, -9.250000, 126.330000, -11.670000, 260.590000, -8.480000, 304.940000, -0.350000, 0.000000, 2.170000, 338.870000, -0.270000, 121.150000, -0.010000, 171.840000, 0.010000, 234.280000, 0.010000, 235.120000, 0.180000, 0.000000, 1.480000, 324.330000, -0.710000, 246.050000, -0.440000, 36.380000] ,
[3.370000, 0.000000, -8.520000, 318.780000, -17.380000, 153.130000, -14.070000, 55.220000, -2.180000, 0.000000, -0.330000, 157.570000, -0.260000, 150.910000, 0.020000, 151.280000, 0.090000, 271.220000, -0.020000, 281.020000, -0.270000, 0.000000, 0.010000, 326.630000, 1.300000, 351.340000, 0.670000, 126.670000] ,
[-0.980000, 0.000000, -9.810000, 19.320000, 28.970000, 40.760000, 1.550000, 45.660000, -0.960000, 0.000000, 1.030000, 247.600000, -0.820000, 101.020000, 0.010000, 57.100000, 0.020000, 181.590000, -0.050000, 163.980000, -1.100000, 0.000000, 1.670000, 238.290000, 0.300000, 241.120000, 0.620000, 142.600000] ,
[5.210000, 0.000000, 0.350000, 109.830000, -9.310000, 72.310000, -7.450000, 197.620000, 0.770000, 0.000000, -0.110000, 306.950000, 0.090000, 344.140000, -0.010000, 74.950000, -0.020000, 239.700000, -0.060000, 337.080000, -0.450000, 0.000000, -0.320000, 154.080000, -1.720000, 226.040000, -0.850000, 114.490000] ,
[4.500000, 0.000000, -9.730000, 339.600000, -16.590000, 166.340000, 6.540000, 198.970000, -1.740000, 0.000000, 1.770000, 351.700000, 0.730000, 279.860000, 0.010000, 136.200000, 0.020000, 133.190000, -0.090000, 349.150000, 1.710000, 0.000000, 1.680000, 74.160000, 0.570000, 107.270000, 1.790000, 115.510000] ,
[2.920000, 0.000000, -2.450000, 80.440000, -21.970000, 127.070000, 12.640000, 151.030000, 2.070000, 0.000000, -2.080000, 172.040000, -0.090000, 26.580000, 0.020000, 333.810000, -0.040000, 102.370000, 0.040000, 210.430000, 0.230000, 0.000000, 1.460000, 248.430000, 1.670000, 165.800000, 1.090000, 267.920000] ,
[-9.130000, 0.000000, -5.160000, 246.360000, 26.900000, 253.420000, -17.770000, 229.310000, -1.930000, 0.000000, -0.040000, 313.870000, -0.320000, 10.320000, -0.010000, 184.300000, -0.010000, 289.630000, -0.080000, 343.330000, 0.170000, 0.000000, 0.580000, 267.540000, -0.280000, 33.330000, 0.230000, 173.560000] ,
[-5.740000, 0.000000, 7.550000, 112.180000, 32.710000, 104.900000, 19.880000, 40.780000, 1.100000, 0.000000, 2.090000, 18.790000, 0.930000, 290.410000, 0.020000, 168.420000, 0.110000, 98.550000, -0.080000, 337.620000, 1.650000, 0.000000, -1.410000, 60.160000, 1.730000, 154.700000, 0.100000, 288.160000] ,
[-4.560000, 0.000000, 3.630000, 27.150000, -7.450000, 183.200000, 12.050000, 191.560000, 0.640000, 0.000000, -0.930000, 286.400000, 0.770000, 142.430000, -0.010000, 95.310000, 0.050000, 277.830000, -0.090000, 13.560000, 2.000000, 0.000000, -1.540000, 44.540000, -0.870000, 312.120000, -0.850000, 106.170000] ,
[-8.490000, 0.000000, 0.940000, 294.840000, 26.020000, 136.040000, 15.520000, 78.250000, 1.090000, 0.000000, 0.720000, 247.320000, -0.030000, 142.150000, 0.020000, 342.330000, -0.050000, 118.890000, -0.080000, 357.400000, -1.970000, 0.000000, 0.680000, 53.070000, -0.960000, 2.130000, 1.570000, 285.040000] ,
[-8.680000, 0.000000, 5.530000, 61.060000, 32.230000, 65.530000, 2.410000, 110.720000, 0.110000, 0.000000, -1.460000, 160.800000, -0.800000, 98.770000, 0.010000, 298.340000, 0.060000, 221.820000, 0.100000, 165.470000, -0.190000, 0.000000, 1.870000, 17.410000, -1.270000, 173.820000, -0.040000, 217.440000] ,
[-2.180000, 0.000000, -9.890000, 130.200000, 2.070000, 173.680000, 6.600000, 285.030000, -1.750000, 0.000000, -0.420000, 80.770000, 0.380000, 224.690000, 0.020000, 204.220000, -0.020000, 289.670000, -0.010000, 204.360000, -1.980000, 0.000000, -1.380000, 132.190000, -0.350000, 172.690000, 1.950000, 41.160000] ,
[4.590000, 0.000000, -8.860000, 162.460000, 20.640000, 41.060000, 3.880000, 38.490000, 1.750000, 0.000000, 1.600000, 93.940000, -0.030000, 359.240000, 0.010000, 96.080000, 0.080000, 56.250000, 0.100000, 187.380000, -1.730000, 0.000000, -0.390000, 137.400000, 1.050000, 41.550000, 0.400000, 207.580000] ,
[4.020000, 0.000000, 5.790000, 200.750000, -32.390000, 147.320000, 0.940000, 307.360000, 2.020000, 0.000000, -1.450000, 43.120000, -0.260000, 309.350000, 0.020000, 10.790000, -0.040000, 311.740000, -0.100000, 172.510000, 2.100000, 0.000000, -0.580000, 353.900000, -0.560000, 31.610000, -1.240000, 276.080000] ,
[-3.280000, 0.000000, -4.660000, 287.390000, -9.500000, 12.380000, -2.260000, 19.280000, 1.780000, 0.000000, -1.360000, 348.080000, -0.900000, 312.480000, 0.010000, 126.910000, -0.020000, 168.230000, -0.030000, 210.400000, -1.830000, 0.000000, 1.310000, 127.940000, -0.670000, 22.770000, -1.460000, 61.700000] ,
[9.210000, 0.000000, 6.660000, 127.520000, -22.910000, 101.420000, -7.370000, 34.040000, -1.530000, 0.000000, 1.640000, 273.490000, 0.500000, 19.940000, 0.020000, 21.530000, -0.010000, 292.750000, -0.040000, 234.230000, -1.110000, 0.000000, -0.100000, 220.930000, 1.980000, 107.510000, -0.580000, 102.360000] ,
[-0.930000, 0.000000, -10.270000, 184.650000, 24.220000, 141.830000, 7.050000, 297.820000, 1.850000, 0.000000, 0.130000, 357.040000, 0.560000, 25.160000, -0.010000, 144.480000, -0.000000, 172.080000, -0.040000, 17.520000, 0.840000, 0.000000, 0.710000, 205.800000, -1.320000, 336.360000, -1.120000, 51.210000] ,
[-1.540000, 0.000000, -0.800000, 214.140000, -25.870000, 317.610000, -8.970000, 170.110000, -0.870000, 0.000000, 2.170000, 342.670000, 0.980000, 79.780000, -0.000000, 188.140000, 0.090000, 250.220000, -0.060000, 150.040000, -0.140000, 0.000000, -1.010000, 167.120000, -1.940000, 59.280000, -1.680000, 333.580000] ];

var hasWASM;

for (var it = 0; it < aberration_list.length; it++) {
    var ab_name = aberration_list[it];
    var ab_obj = {};
    ab_obj.m = Number(ab_name[1]);
    ab_obj.n = Number(ab_name[2]);
    ab_obj.mag_el = document.getElementById(ab_name);
    if (ab_obj.m == 1 && ab_obj.n == 0) {
        ab_obj.mag_unit = ang;
    } else if (ab_obj.m < 4) {
        ab_obj.mag_unit = nm;
    } else if (ab_obj.m < 5) {
        ab_obj.mag_unit = um;
    } else {
        ab_obj.mag_unit = mm;
    }
    if (ab_obj.n != 0) {
        ab_obj.arg_el = document.getElementById("P" + ab_obj.m + ab_obj.n);
    }
    ab_obj.mag_unit = ab_obj.mag_unit * correction_factor;
    aberrations.push(ab_obj);
}

window.addEventListener(
    "keydown",
    function (event) {
        if (event.defaultPrevented) {
            return; // Do nothing if the event was already processed
        }

        switch (event.key) {
            case "ArrowDown":
                if (interactiveMode.checked) {
                    aberrations[0].mag_el.value =
                        Number(aberrations[0].mag_el.value) - 10;
                    calculate();
                    break;
                }
            case "ArrowUp":
                if (interactiveMode.checked) {
                    aberrations[0].mag_el.value =
                        Number(aberrations[0].mag_el.value) + 10;
                    calculate();
                    break;
                }
            case "w":
                if (interactiveMode.checked) {
                    sx =
                        Number(aberrations[1].mag_el.value) *
                        math.cos(Number(aberrations[1].arg_el.value) * deg);
                    sy =
                        Number(aberrations[1].mag_el.value) *
                        math.sin(Number(aberrations[1].arg_el.value) * deg);

                    sy += 1;
                    aberrations[1].mag_el.value = math.sqrt(sx * sx + sy * sy);
                    aberrations[1].arg_el.value = math.atan2(sy, sx) / deg;
                    calculate();
                    break;
                }
            case "s":
                if (interactiveMode.checked) {
                    sx =
                        Number(aberrations[1].mag_el.value) *
                        math.cos(Number(aberrations[1].arg_el.value) * deg);
                    sy =
                        Number(aberrations[1].mag_el.value) *
                        math.sin(Number(aberrations[1].arg_el.value) * deg);

                    sy -= 1;
                    aberrations[1].mag_el.value = math.sqrt(sx * sx + sy * sy);
                    aberrations[1].arg_el.value = math.atan2(sy, sx) / deg;
                    calculate();
                    break;
                }
            case "a":
                if (interactiveMode.checked) {
                    sx =
                        Number(aberrations[1].mag_el.value) *
                        math.cos(Number(aberrations[1].arg_el.value) * deg);
                    sy =
                        Number(aberrations[1].mag_el.value) *
                        math.sin(Number(aberrations[1].arg_el.value) * deg);

                    sx -= 1;
                    aberrations[1].mag_el.value = math.sqrt(sx * sx + sy * sy);
                    aberrations[1].arg_el.value = math.atan2(sy, sx) / deg;
                    calculate();
                    break;
                }
            case "d":
                if (interactiveMode.checked) {
                    sx =
                        Number(aberrations[1].mag_el.value) *
                        math.cos(Number(aberrations[1].arg_el.value) * deg);
                    sy =
                        Number(aberrations[1].mag_el.value) *
                        math.sin(Number(aberrations[1].arg_el.value) * deg);

                    sx += 1;
                    aberrations[1].mag_el.value = math.sqrt(sx * sx + sy * sy);
                    aberrations[1].arg_el.value = math.atan2(sy, sx) / deg;
                    calculate();
                    break;
                }
            case "i":
                if (interactiveMode.checked) {
                    cx =
                        Number(aberrations[2].mag_el.value) *
                        math.cos(Number(aberrations[2].arg_el.value) * deg);
                    cy =
                        Number(aberrations[2].mag_el.value) *
                        math.sin(Number(aberrations[2].arg_el.value) * deg);

                    cy += 100;
                    aberrations[2].mag_el.value = math.sqrt(cx * cx + cy * cy);
                    aberrations[2].arg_el.value = math.atan2(cy, cx) / deg;
                    calculate();
                    break;
                }
            case "k":
                if (interactiveMode.checked) {
                    cx =
                        Number(aberrations[2].mag_el.value) *
                        math.cos(Number(aberrations[2].arg_el.value) * deg);
                    cy =
                        Number(aberrations[2].mag_el.value) *
                        math.sin(Number(aberrations[2].arg_el.value) * deg);

                    cy -= 100;
                    aberrations[2].mag_el.value = math.sqrt(cx * cx + cy * cy);
                    aberrations[2].arg_el.value = math.atan2(cy, cx) / deg;
                    calculate();
                    break;
                }
            case "j":
                if (interactiveMode.checked) {
                    cx =
                        Number(aberrations[2].mag_el.value) *
                        math.cos(Number(aberrations[2].arg_el.value) * deg);
                    cy =
                        Number(aberrations[2].mag_el.value) *
                        math.sin(Number(aberrations[2].arg_el.value) * deg);

                    cx += 100;
                    aberrations[2].mag_el.value = math.sqrt(cx * cx + cy * cy);
                    aberrations[2].arg_el.value = math.atan2(cy, cx) / deg;
                    calculate();
                    break;
                }
            case "l":
                if (interactiveMode.checked) {
                    cx =
                        Number(aberrations[2].mag_el.value) *
                        math.cos(Number(aberrations[2].arg_el.value) * deg);
                    cy =
                        Number(aberrations[2].mag_el.value) *
                        math.sin(Number(aberrations[2].arg_el.value) * deg);

                    cx -= 100;
                    aberrations[2].mag_el.value = math.sqrt(cx * cx + cy * cy);
                    aberrations[2].arg_el.value = math.atan2(cy, cx) / deg;
                    calculate();
                    break;
                }
            default:
                return; // Quit when this doesn't handle the key event.
        }

        // Cancel the default action to avoid it being handled twice
        event.preventDefault();
    },
    true
);

initialize();
