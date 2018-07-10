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


function calculate(){
    var lambda = 12.3986/math.sqrt((2*511+300)*300)*ang; // 300 keV
    var al_max = 70*math.pow(10,-3);
    var al_vec = math.matrix(math.range(-al_max,al_max,(2*al_max)/(numPx)));
    al_vec.resize([256,1])

    var alxx = math.multiply(math.ones(256,1),math.transpose(al_vec));
    var alyy = math.transpose(alxx);


    var alrr = math.sqrt(math.add(math.dotPow(alxx,2),math.dotPow(alyy,2)));
    var alpp = math.atan2(alyy,alxx);

    var obj_ap_r = 50 * math.pow(10,-3);
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

    var trans = math.exp(  math.multiply(math.complex(0,-1),PI,.25, math.matrix(math.random([numPx,numPx])))  );
    
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
        if(value < PI/4)
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
    rmax = -1/rmax*1000; //mrads

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

    canvas = document.getElementById("canvas2");
    ctx = canvas.getContext("2d");        
    canvas.width = numPx;
    canvas.height = numPx;
    drawGrayscaleBitmap(ctx,out_phase_map);    

}

function randomize(){

    for(var it = 0; it < aberrations.length; it++)
    {
        var aberration = aberrations[it];
        aberration.mag_el.value = Math.random()*100;
        if(aberration.arg_el)
        {
            aberration.arg_el.value = Math.random()*100;
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
var numPx = 256;

var ang = math.pow(10,-10);
var nm = math.pow(10,-9);
var um = math.pow(10,-6);
var mm = math.pow(10,-3);
var PI = math.pi;
var deg = PI/180;


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
    aberrations.push(ab_obj);
}

window.onload = calculate();

