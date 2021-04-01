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
    rmax
) {
    let scalar = 256;
    ctx1.font = (numPx / scalar) * 14 + "px Arial";
    ctx1.fillStyle = "white";
    ctx1.fillText(
        math.round((disp_size_mrad / 0.07) * 30) + " mrad",
        numPx - (70 / scalar) * numPx,
        numPx - (10 / scalar) * numPx
    );

    ctx1.beginPath();
    ctx1.moveTo(numPx - (70 / scalar) * numPx, numPx - (30 / scalar) * numPx);
    ctx1.lineTo(numPx - (15 / scalar) * numPx, numPx - (30 / scalar) * numPx);
    ctx1.strokeStyle = "white";
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
    ctx1.strokeStyle = "blue";
    ctx1.lineWidth = (1 * numPx) / scalar;
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
    ctx2.strokeStyle = "blue";
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
    ctx2.strokeStyle = "red";
    ctx2.lineWidth = 2;
    ctx2.stroke();
}

function drawGrayscaleBitmap(ctx, bitmap, numPx) {
    for (let it = 0; it < numPx; it++) {
        for (let jt = 0; jt < numPx; jt++) {
            let value = bitmap[it][jt];
            if (value < 0) {
                value = 0;
            } else if (value >= 256) {
                value = 255;
            }
            let part = Number(parseInt(value, 10)).toString(16);
            if (part.length < 2) {
                part = "0" + part;
            }
            let color = "#" + part + part + part;
            ctx.fillStyle = color;
            ctx.fillRect(it, jt, 1, 1);
        }
    }
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
    setTimeout(function () {
        randomize();
    }, 0);
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

function drawEverything(
    out_ronch,
    out_phase_map,
    numPx,
    al_max,
    disp_size_mrad,
    obj_ap_r,
    rmax
) {
    // drawing part
    canvas1.width = numPx;
    canvas1.height = numPx;
    canvas2.width = numPx;
    canvas2.height = numPx;
    drawGrayscaleBitmap(ctx1, out_ronch, numPx);
    drawGrayscaleBitmap(ctx2, out_phase_map, numPx);
    if (draw_overlay) {
        drawOverlays(ctx1, ctx2, numPx, al_max, disp_size_mrad, obj_ap_r, rmax);
    }

    document.getElementById("alpha_max").value = math.round(rmax, 2);
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

    let out_ronch = math.dotPow(
        math.abs(
            math.dotMultiply(
                math.matrix(
                    fft2_wrap(
                        math
                            .dotMultiply(
                                trans,
                                math.matrix(fft2_wrap(chi.toArray()))
                            )
                            .toArray()
                    )
                ),
                obj_ap
            )
        ),
        2
    );

    out_ronch = math.round(normalizeScale(out_ronch, 255));
    out_ronch = out_ronch.toArray();

    let returnVals = maskChi0(chi0, alrr, numPx);
    let out_phase_map = returnVals[0];
    let rmax = returnVals[1];

    drawEverything(
        out_ronch,
        out_phase_map,
        numPx,
        al_max,
        disp_size_mrad,
        obj_ap_r,
        rmax
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
    let out_ronch = [];
    let out_phase_map = [];
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
        }
        out_ronch.push(arrayData1);
        out_phase_map.push(arrayData2);
        arrayData1 = [];
        arrayData2 = [];
    }
    let rmax =
        Module.HEAPF32[
            result / Float32Array.BYTES_PER_ELEMENT + 2 * (numPx * numPx)
        ];

    drawEverything(
        out_ronch,
        out_phase_map,
        numPx,
        al_max,
        disp_size_mrad,
        obj_ap_r,
        rmax
    );
}

function calculate() {
    let t0 = performance.now();
    energyUI();
    if (hasWASM && !forceJS.checked) {
        document.getElementById("loading").innerHTML =
            "Calculating with WebAssembly...";
        console.log("Calculating with WebAssembly...");
        let curInstance = ronchModule().then(function (Module) {
            calculateWASM(Module);
            Module.delete;
        });
    } else {
        document.getElementById("loading").innerHTML =
            "Calculating with Javascript...";
        console.log("Calculating with Javascript...");
        calculateJS();
    }
    console.log("T = " + (performance.now() - t0) + " ms");
    document.getElementById("loading").innerHTML = " ";
}

function initialize() {
    hasWASM = hasWASM();
    calculate();
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

var forceJS = document.getElementById("forceJS"); //figure out how to read from checkbox
var interactiveMode = document.getElementById("interactiveMode");
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
var aberrations = [];

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

                    cy += 1;
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

                    cy -= 1;
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

                    cx += 1;
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

                    cx += 1;
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
