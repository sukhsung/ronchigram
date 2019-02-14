//from https://web.eecs.umich.edu/~fessler/course/451/l/pdf/c6.pdf
/*function radix_fft(x)
{

    n = math.matrix(math.range(0,numPx)); //0...255
    W_re = math.cos(math.multiply(2*math.pi/numPx,n));
    W_im = math.sin(math.multiply(2*math.pi/numPx,n));// - i * math.sin(2*math.pi/numPx*n);

    W = math.zeros(numPx);
    W = W.map(function(value, index, matrix){
        return math.complex(W_re.subset(math.index(index)), W_im.subset(math.index(index)));
    });
    X = math.zeros(numPx);
    X = X.map(function(value, index, matrix){

        k = index;
        //console.log(k);
        sum = 0;
        for(var n = 0; n < numPx -1; n++)
        {
            m = math.mod(k*n,numPx);
            sum = math.add(sum, math.multiply(x.subset(math.index(n)), W.subset(math.index(m))  ));
        }
        //m = k*n;
        //sum = math.dotMultiply(x,W)
        return sum
    });
    return X;
    //W = math.complex(W_re, W_im);
}

function radix_fft2(x){

}*/


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


function fft2_wrap(X) {
    //X_arr = X.toArray();
    X_pass = fft2(X);
    X_arr = math.matrix(X_pass);
    X_arr = math.transpose(X_arr);
    X_arr = X_arr.toArray();
    X_pass = fft2(X_arr);
    return (X_pass);
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

function drawGrayscaleBitmap(ctx,bitmap)
{
    for(var it = 0; it < numPx; it++)
    {
        for(var jt = 0; jt < numPx; jt++)
        {
            var value = bitmap[it][jt];
            var part = Number(parseInt( value , 10)).toString(16);
            if(part.length <2 )
            {
                part = "0"+part;
            }
            var color = '#' + part+part+part;
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
    document.getElementById("wavlen").value = (lambda/pm);

    var alpha = Number(document.getElementById("aperture").value)* mrad;
    //resolution calculation:
    var d = .61*lambda/pm/alpha;
    document.getElementById("diffres").value = d;
    return lambda;
}

function calcButton(){
     document.getElementById('loading').innerHTML = "Calculating..."
     setTimeout(function(){
        calculate();
    },0);
}

function randButton(){
     document.getElementById('loading').innerHTML = "Calculating..."
     setTimeout(function(){
        randomize();
    },0);
}


function loadSample(){
    var scalefactor = 8;
         //subsample = math.dotMultiply(math.subtract(math.random([numPx/scalefactor,numPx/scalefactor]),0.5),2);

    var subsample = math.random([numPx/scalefactor,numPx/scalefactor]);

    var supersample = math.zeros(numPx,numPx);
    //quick nearest neightbours interpolation
    supersample = supersample.map(function(value,index,matrix){
        return subsample[math.floor(index[0]/scalefactor)][math.floor(index[1]/scalefactor)];
    });

    return supersample;



}


function calculate(){
    lambda = energyCalc();
    ////////
    //reading in constants from ui:
    ////////
    var obj_ap_r = Number(document.getElementById("aperture").value)* mrad;

    if(obj_ap_r<0)
    {
        obj_ap_r = 0;
        document.getElementById("aperture").value = 0;
    }
    else if (obj_ap_r>65*mrad)
    {
        obj_ap_r= 65*mrad;
        document.getElementById("aperture").value = 65;
    }

    var disp_size_px = Number(document.getElementById("disp_size_px").value);
    if((disp_size_px & (disp_size_px - 1)) != 0  || disp_size_px < 2)
    {
        alert("Select a display size in pixels that is a power of 2 greater than 0");
        return;
    }
    else
    {
        numPx = disp_size_px;
    }

    var disp_size_mrad = Number(document.getElementById("disp_size_mrad").value)*mrad/2;
    if(disp_size_mrad<.0000001  )
    {
        disp_size_mrad = .0000001
    }


    //var ill_angle = Number(document.getElementById("illumination").value)*mrad;
    var draw_overlay = document.getElementById("draw_overlay").checked; //figure out how to read from checkbox

    var al_max = disp_size_mrad;//70*mrad; //= ill_angle;
    var al_vec = math.matrix(math.range(-al_max,al_max,(2*al_max)/(numPx)));
    al_vec.resize([numPx,1])

    var alxx = math.multiply(math.ones(numPx,1),math.transpose(al_vec));
    var alyy = math.transpose(alxx);

    var alrr = math.sqrt(math.add(math.dotPow(alxx,2),math.dotPow(alyy,2)));
    var alpp = math.atan2(alyy,alxx);



    var obj_ap = alrr.map(function (value, index, matrix) {
        if(value < obj_ap_r)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    });

    var sample = loadSample();
    var trans = math.exp(  math.multiply(math.complex(0,-1),PI,.25, sample)  );
    
    var aber = getAberrations();
    var numAber = aber.size()[0];

    var chi = math.zeros(numPx,numPx);

    for(var it = 0; it < numAber; it++)
    {
        chi = math.add(chi, math.dotMultiply(math.dotMultiply(math.cos(math.dotMultiply(aber.subset(math.index(it,1)),math.subtract(alpp,aber.subset(math.index(it,3))))),math.dotPow(alrr,aber.subset(math.index(it,0))+1)), aber.subset(math.index(it,2))/(aber.subset(math.index(it,0))+1) ));
    }
    var chi0 = math.dotMultiply(2*PI/lambda, chi);
    var expchi0 = math.dotMultiply(math.dotPow(math.E, math.dotMultiply(math.complex(0,-1),chi0) ), obj_ap);
    var out_ronch = math.dotPow(math.abs(math.dotMultiply(math.matrix(fft2_wrap(math.dotMultiply(trans,math.matrix(fft2_wrap(expchi0.toArray()))).toArray())),obj_ap)),2);

    out_ronch = math.subtract(out_ronch, math.min(out_ronch));
    out_ronch = math.dotDivide(out_ronch,math.max(out_ronch)/255);
    out_ronch = math.round(out_ronch);
    out_ronch = out_ronch.toArray();
    var out_phase_map = chi0.map(function (value, index, matrix) {
        if(value < PI/4 && value > -PI/4)
        {
            return 1;            
        }
        else
        {
            return 0;
        }
    });

    var rmax = math.dotDivide(1,math.dotMultiply(alrr,math.subtract(out_phase_map,1)));
    rmax = math.min(rmax);
    rmax = -1/(rmax*mrad); //mrads

    document.getElementById("alpha_max").value = math.round(rmax,2);

    out_phase_map = math.abs(out_phase_map);
    out_phase_map = math.subtract(out_phase_map, math.min(out_phase_map));
    out_phase_map = math.dotDivide(out_phase_map,math.max(out_phase_map)/255);
    out_phase_map = math.round(out_phase_map);
    out_phase_map = out_phase_map.toArray();


    var canvas = document.getElementById("canvas1");
    var ctx = canvas.getContext("2d");
    canvas.width = numPx;
    canvas.height = numPx;
    drawGrayscaleBitmap(ctx,out_ronch);
    if(draw_overlay)
    {
        var scalar = 256;
        ctx.font = numPx/scalar*14+"px Arial";
        ctx.fillStyle = "white";
        ctx.fillText(math.round(disp_size_mrad/.07*30)+" mrad",numPx-70/scalar*numPx,numPx-10/scalar*numPx);

        ctx.beginPath()

        ctx.moveTo(numPx-70/scalar*numPx,numPx-30/scalar*numPx);
        ctx.lineTo(numPx-15/scalar*numPx,numPx-30/scalar*numPx);
        ctx.strokeStyle = "white";
        ctx.lineWidth = 5*numPx/scalar;
        ctx.stroke();
        ctx.beginPath()
        ctx.arc(numPx/2,numPx/2,rmax*numPx/(2*al_max)*mrad,0,2*PI);
        ctx.strokeStyle = "blue";
        ctx.lineWidth = 1*numPx/scalar;
        ctx.stroke();
    }

    canvas = document.getElementById("canvas2");
    ctx = canvas.getContext("2d");        
    canvas.width = numPx;
    canvas.height = numPx;
    drawGrayscaleBitmap(ctx,out_phase_map);
    if(draw_overlay)
    {
        ctx.beginPath();
        ctx.arc(numPx/2,numPx/2,rmax*numPx/(2*al_max)*mrad,0,2*PI);
        ctx.strokeStyle = "blue";
        ctx.lineWidth = 2;
        ctx.stroke();
        ctx.beginPath();
        ctx.arc(numPx/2,numPx/2,obj_ap_r*numPx/(2*al_max),0,2*PI);
        ctx.strokeStyle = "red";
        ctx.lineWidth = 2;
        ctx.stroke();
    }


    document.getElementById('loading').innerHTML = " "
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
    calculate();
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
//TODO: parses URL, returning object storing aberration identifiers and values
function parseURL(){

    return {};
}
//TODO: loads parsed aberrations into UI or sets defaults
function loadAberrations(){
    parsed_abs = parseURL();
    //logic to go thru parsed abs and set in UI
}

//TODO: encodes current aberration values into a string, outputs to console (for now)
function generateURL(){
    var str = "?";
     ab_snapshot = [];
    for(var it = 0; it < aberrations.length; it++)
    {
        var aberration = aberrations[it];
        var mag_val = aberration.mag_el.value;
        var arg_val = (aberration.arg_el ? aberration.arg_el.value : "");
        //more?
        ab_snapshot.push([aberration.m, aberration.n, mag_val, arg_val]);

    }
    //ab_snap -> string -> base64 or some encoding?

    return str;
}



var numPx = 256;
var pm = math.pow(10,-12);
var ang = math.pow(10,-10);
var nm = math.pow(10,-9);
var um = math.pow(10,-6);
var mm = math.pow(10,-3);
var mrad = mm;
var PI = math.pi;
var deg = PI/180;
var correction_factor = 1;

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


