function fft2_wrap(X) {
    X_pass = fft2(X);
    X_arr = math.matrix(X_pass);
    X_arr = math.transpose(X_arr);
    X_arr = X_arr.toArray();
    X_pass = fft2(X_arr);
    return (X_pass);
}
/*
    From https://gist.github.com/mrquincle/b11fff96209c9d1396b0
    @mrquincle
*/
function fft2(X) {
  var N = X.length;
  if (!(N > 1)) {
    return X;
  }
  var M = N/2;
  var even = [];
  var odd = [];
  even.length = M;
  odd.length = M;
  for (var i = 0; i < M; ++i) {
    even[i] = X[i*2];
    odd[i] = X[i*2+1];
  }
  even = fft2(even);
  odd = fft2(odd);
  var a = -2*PI;
  for (var k = 0; k < M; ++k) {
    var t = math.exp(math.complex(0, a*k/N));
    t = math.multiply(t, odd[k]);
    X[k] = odd[k] = math.add(even[k], t);
    X[k+M] = even[k] = math.subtract(even[k], t);
  }
  return X;
}

function normalizeScale(matharray,scale) {
    let output = math.subtract(matharray,math.min(matharray))
    output = math.dotDivide(output,math.max(output)/scale)
    return output;
}

function drawOverlays(ctx1, ctx2, numPx,al_max, disp_size_mrad, obj_ap_r, rmax) {
    let scalar = 256;
    ctx1.font = numPx/scalar*14+"px Arial";
    ctx1.fillStyle = "white";
    ctx1.fillText(math.round(disp_size_mrad/.07*30)+" mrad",numPx-70/scalar*numPx,numPx-10/scalar*numPx);

    ctx1.beginPath()
    ctx1.moveTo(numPx-70/scalar*numPx,numPx-30/scalar*numPx);
    ctx1.lineTo(numPx-15/scalar*numPx,numPx-30/scalar*numPx);
    ctx1.strokeStyle = "white";
    ctx1.lineWidth = 5*numPx/scalar;
    ctx1.stroke();
    ctx1.beginPath()
    ctx1.arc(numPx/2,numPx/2,rmax*numPx/(2*al_max)*mrad,0,2*PI);
    ctx1.strokeStyle = "blue";
    ctx1.lineWidth = 1*numPx/scalar;
    ctx1.stroke();

    /// right panel
    ctx2.beginPath();
    ctx2.arc(numPx/2,numPx/2,rmax*numPx/(2*al_max)*mrad,0,2*PI);
    ctx2.strokeStyle = "blue";
    ctx2.lineWidth = 2;
    ctx2.stroke();
    ctx2.beginPath();
    ctx2.arc(numPx/2,numPx/2,obj_ap_r*numPx/(2*al_max),0,2*PI);
    ctx2.strokeStyle = "red";
    ctx2.lineWidth = 2;
    ctx2.stroke();

}

function getAberrations(){
    var ab_list = [];

    for(var it = 0; it < aberrations.length; it++)
    {
        
        var aberration = aberrations[it];
        var mag_val = Number(aberration.mag_el.value)*aberration.mag_unit;
        var arg_val = (aberration.arg_el ? Number(aberration.arg_el.value) : 0)*deg;
        ab_list.push([aberration.m, aberration.n, mag_val, arg_val]);
    }

    return math.matrix(ab_list);

}

function drawGrayscaleBitmap(ctx,bitmap,numPx) {
    for(let it = 0; it < numPx; it++)
    {
        for(let jt = 0; jt < numPx; jt++)
        {
            let value = bitmap[it][jt];
            if (value < 0) {
                value = 0;
            } else if (value >= 256) {
                value = 255;
            }
            let part = Number(parseInt( value , 10)).toString(16);
            if(part.length <2 )
            {
                part = "0"+part;
            }
            let color = '#' + part+part+part;
            ctx.fillStyle=color;
            ctx.fillRect(it,jt,1,1);
        } 
    }
}

function energyCalc(){
    var keV = Number(document.getElementById("beamvolt").value);

    if(keV<0)
    {
        keV=0;
        document.getElementById("beamvolt").value = 0;
    }

    var lambda = 12.3986/Math.sqrt((2*511+keV)*keV) *ang;
    document.getElementById("wavlen").value = math.round(lambda/pm,4);

    var alpha = Number(document.getElementById("aperture").value)* mrad;
    //resolution calculation:
    var d = .61*lambda/pm/alpha;
    document.getElementById("diffres").value = math.round(d,4);
    return lambda;
}

function randButton(){
     //document.getElementById('loading').innerHTML = "Calculating..."
     setTimeout(function(){
        randomize();
    },0);
}

function loadSample(scalefactor,numPx){
    var subsample = math.random([numPx/scalefactor,numPx/scalefactor]);
    var supersample = math.zeros(numPx,numPx);
    //quick nearest neightbours interpolation
    supersample = supersample.map(function(value,index,matrix){
        return subsample[math.floor(index[0]/scalefactor)][math.floor(index[1]/scalefactor)];
    });
    return supersample;
}

//Normalized to 300 keV
function interactionParam(){
    var c = 3e8;
    var mass_e = 9.11e-31;
    var charge_e = 1.602e-19;
    var lambda = energyCalc();
    var keV = Number(document.getElementById("beamvolt").value);

    var param = 2*PI/(lambda*keV/charge_e*1000)*(mass_e*c*c+keV*1000)/(2*mass_e*c*c+keV*1000);
    var param_300 = 2*PI/(lambda*300/charge_e*1000)*(mass_e*c*c+300*1000)/(2*mass_e*c*c+300*1000);
    return param/param_300;
}

function getObjAperture() {
    let obj_ap_r = Number(document.getElementById("aperture").value)* mrad;
    
    if(obj_ap_r<0) {
        obj_ap_r = 0;
        document.getElementById("aperture").value = 0;
    }
    return obj_ap_r
}

function getDispSizePx() {
    let disp_size_px = Number(document.getElementById("disp_size_px").value);
    if((disp_size_px & (disp_size_px - 1)) != 0  || disp_size_px < 2)
    {
        alert("Select a display size in pixels that is a power of 2 greater than 0");
        return;
    }
    else
    {
        return disp_size_px
    }
}

function getDispSizeMrad() {
    let disp_size_mrad = Number(document.getElementById("disp_size_mrad").value)*mrad/2;
    if(disp_size_mrad<.0000001  )
    {
        disp_size_mrad = .0000001
    }
    return disp_size_mrad
}

function hasWASM()
{
    //per @JF-Bastien https://stackoverflow.com/questions/47879864/how-can-i-check-if-a-browser-supports-webassembly

    const supported = (() => {
        try {
            if (typeof WebAssembly === "object"
                && typeof WebAssembly.instantiate === "function") {
                const module = new WebAssembly.Module(Uint8Array.of(0x0, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00));
                if (module instanceof WebAssembly.Module)
                    return new WebAssembly.Instance(module) instanceof WebAssembly.Instance;
            }
        } catch (e) {
        }
        return false;
    })();

    return supported;
}

function calcButton(){
    let t0 = performance.now();

    if(hasWASM())
    {
        document.getElementById('loading').innerHTML = "Calculating with WebAssembly..."
        console.log("wasm supported");
        let curInstance = ronchModule().then(function(Module){ calculateWASM(Module); Module.delete });
    }
    else
    {
        document.getElementById('loading').innerHTML = "Calculating with Javascript..."
        console.log("wasm not supported");
        calculateJS();
    }
    console.log("dT="+(performance.now()-t0)+" ms");

}


function calculateJS(){
    ////////
    //reading in constants from ui:
    ////////
    let lambda = energyCalc()
    let obj_ap_r = getObjAperture()
    let numPx = getDispSizePx()
    let disp_size_mrad = getDispSizeMrad()
    let scalefactor = Number(document.getElementById("sample_scale_factor").value);
    let draw_overlay = document.getElementById("draw_overlay").checked; //figure out how to read from checkbox
    let al_max = disp_size_mrad;

    let al_vec = math.matrix(math.range(-al_max,al_max,(2*al_max)/(numPx)));
    al_vec.resize([numPx,1])

    let alxx = math.multiply(math.ones(numPx,1),math.transpose(al_vec));
    let alyy = math.transpose(alxx);

    let alrr = math.sqrt(math.add(math.dotPow(alxx,2),math.dotPow(alyy,2)));
    let alpp = math.atan2(alyy,alxx);

    let obj_ap = alrr.map(function (value, index, matrix) {
        if(value < obj_ap_r)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    });

    let out_ronch = math.zeros(numPx,numPx);
    let sample = loadSample(scalefactor,numPx);
    let trans = math.exp(  math.multiply(math.complex(0,-1),PI,.25,interactionParam(), sample)  );
    
    let aber = getAberrations();
    let numAber = aber.size()[0];

    let chi = math.zeros(numPx,numPx);

    for(let it = 0; it < numAber; it++)
    {
        chi = math.add(chi, math.dotMultiply(math.dotMultiply(math.cos(math.dotMultiply(aber.subset(math.index(it,1)),math.subtract(alpp,aber.subset(math.index(it,3))))),math.dotPow(alrr,aber.subset(math.index(it,0))+1)), aber.subset(math.index(it,2))/(aber.subset(math.index(it,0))+1) ));
    }
    let chi0 = math.dotMultiply(2*PI/lambda, chi);
    //To place objective before sample:
    //var expchi0 = math.dotMultiply(math.dotPow(math.E, math.dotMultiply(math.complex(0,-1),chi0) ), obj_ap);
    let expchi0 = math.dotPow(math.E, math.dotMultiply(math.complex(0,-1),chi0) );
    out_ronch = math.add(out_ronch,  math.dotPow(math.abs(math.dotMultiply(math.matrix(fft2_wrap(math.dotMultiply(trans,math.matrix(fft2_wrap(expchi0.toArray()))).toArray())),obj_ap)),2));

    out_ronch = math.round(normalizeScale(out_ronch,255));
    out_ronch = out_ronch.toArray();
    let out_phase_map = chi0.map(function (value, index, matrix) {
        if(value < PI/4 && value > -PI/4)
        {
            return 1;            
        }
        else
        {
            return 0;
        }
    });

    let rmax = math.dotDivide(1,math.dotMultiply(alrr,math.subtract(out_phase_map,1)));
    rmax = math.min(rmax);
    rmax = -1/(rmax*mrad); //mrads

    out_phase_map = math.round(normalizeScale(math.abs(out_phase_map),255));
    out_phase_map = out_phase_map.toArray();
    
    canvas1.width = numPx;
    canvas1.height = numPx;
    canvas2.width = numPx;
    canvas2.height = numPx;
    drawGrayscaleBitmap(ctx1,out_ronch,numPx);
    drawGrayscaleBitmap(ctx2,out_phase_map,numPx);    
    
    if(draw_overlay) {
        drawOverlays(ctx1, ctx2, numPx,al_max,disp_size_mrad, obj_ap_r, rmax)
    }

    document.getElementById('loading').innerHTML = " "  
    document.getElementById("alpha_max").value = math.round(rmax,2);
}

function calculateWASM(Module){
    ////////
    //reading in constants from ui:
    ////////
    let lambda = energyCalc()
    let obj_ap_r = getObjAperture()
    let numPx = getDispSizePx()
    let disp_size_mrad = getDispSizeMrad()
    let scalefactor = Number(document.getElementById("sample_scale_factor").value);
    let draw_overlay = document.getElementById("draw_overlay").checked; //figure out how to read from checkbox
    let al_max = disp_size_mrad;

    // getting aberrations into tidy arrays of magnitude, angle. degree, order  are assumed based on order in C++ section, units are baked in!
    let ab_mags = [];
    let ab_angles = [];

    for(let it = 0; it < aberrations.length; it++)
    {
        let aberration = aberrations[it];
        let mag_val = Number(aberration.mag_el.value)*aberration.mag_unit;
        let arg_val = (aberration.arg_el ? Number(aberration.arg_el.value) : 0)*deg;
        ab_mags.push(mag_val);
        ab_angles.push(arg_val);
    }

    let params = [numPx,al_max,obj_ap_r];
    const arrayDataToPass = params.concat(ab_mags,ab_angles);
    let buffer
    let error
    let result
    try {
        const typedArray = new Float32Array(arrayDataToPass.length)
        for (let i=0; i<arrayDataToPass.length; i++) {
                typedArray[i] = arrayDataToPass[i]
            }
        buffer = Module._malloc(typedArray.length * typedArray.BYTES_PER_ELEMENT)
        Module.HEAPF32.set(typedArray, buffer >> 2)
        result = Module.ccall("calcRonch", null, ["number", "number"], [buffer, arrayDataToPass.length])
    } catch (e) {
        error = e
    } finally {
        // To avoid memory leaks we need to always clear out the allocated heap data
        // This needs to happen in the finally block, otherwise thrown errors will stop code execution before this happens
        Module._free(buffer)
    }
    if (error) throw error

    let arrayData1 =[]
    let arrayData2 = []
    let imData1 = []
    let imData2 = []
    let im2Offset = numPx*numPx;
    for (let j=0; j<numPx;j++) {
        for (let i=0; i<numPx; i++) {
            arrayData1.push(Module.HEAPF32[result/Float32Array.BYTES_PER_ELEMENT+ i+numPx*j])
            arrayData2.push(Module.HEAPF32[result/Float32Array.BYTES_PER_ELEMENT+ i+numPx*j + im2Offset])
        }
        imData1.push(arrayData1)
        imData2.push(arrayData2)
        arrayData1 = []
        arrayData2 = []
    }
    let rmax = Module.HEAPF32[result/Float32Array.BYTES_PER_ELEMENT+ 2*(numPx*numPx)];
    
    canvas1.width = numPx;
    canvas1.height = numPx;
    canvas2.width = numPx;
    canvas2.height = numPx;
    drawGrayscaleBitmap(ctx1,imData1,numPx);
    drawGrayscaleBitmap(ctx2,imData2,numPx);

    if(draw_overlay) {
        drawOverlays(ctx1, ctx2, numPx,al_max,disp_size_mrad, obj_ap_r, rmax)
    }

    document.getElementById('loading').innerHTML = " "
    document.getElementById("alpha_max").value = math.round(rmax,2);
}

function randomize(){

    for(var it = 0; it < aberrations.length; it++)
    {
        var aberration = aberrations[it];
        aberration.mag_el.value = Math.round(Math.random()*100);
        if(aberration.arg_el)
        {
            aberration.arg_el.value = Math.round(Math.random()*180);
        }
    }
    calcButton();
}

function allZero(){

    for(var it = 0; it < aberrations.length; it++)
    {
        var aberration = aberrations[it];
        aberration.mag_el.value = 0;
        if(aberration.arg_el)
        {
            aberration.arg_el.value = 0;
        }
    }

}

function setC(c_in){
    for(var it = 0; it < aberrations.length; it++)
    {
        if(c_in == ""+aberrations[it].m+aberrations[it].n)
        {
            aberrations[it].mag_el.value = Number(aberrations[it].mag_el.value) + 50;
        }
    }
}

var pm = math.pow(10,-12);
var ang = math.pow(10,-10);
var nm = math.pow(10,-9);
var um = math.pow(10,-6);
var mm = math.pow(10,-3);
var mrad = math.pow(10,-3);
var PI = math.pi;
var deg = PI/180;
var correction_factor = 1;

var canvas1 = document.getElementById("canvas1");
var ctx1 = canvas1.getContext("2d");
var canvas2 = document.getElementById("canvas2");
var ctx2 = canvas2.getContext("2d");

var aberration_list = ["C10","C12","C21","C23","C30","C32","C34","C41","C43","C45","C50","C52","C54","C56"];
var aberrations = [];

for(var it = 0; it < aberration_list.length; it++)
{
    var ab_name = aberration_list[it];
    var ab_obj = {};
    ab_obj.m = Number(ab_name[1]);
    ab_obj.n = Number(ab_name[2]);
    ab_obj.mag_el = document.getElementById(ab_name);
    if(ab_obj.m==1 && ab_obj.n ==0)
    {
        ab_obj.mag_unit = ang;
    }
    else if(ab_obj.m < 4)
    {
        ab_obj.mag_unit = nm;
    }
    else if(ab_obj.m < 5)
    {
        ab_obj.mag_unit = um;
    }
    else
    {
        ab_obj.mag_unit = mm;
    }
    if(ab_obj.n != 0)
    {
        ab_obj.arg_el = document.getElementById("P"+ab_obj.m+ab_obj.n)
    }
    ab_obj.mag_unit = ab_obj.mag_unit * correction_factor;
    aberrations.push(ab_obj);
}

calcButton();