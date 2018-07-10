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
//C10, C12, C21, C23, C30
    //C10

    //Number(parseInt( value , 10))

    //C10,C21, P21,...
    var C10_in1 = Number(document.getElementById("C10").value);
    var C12_in1 = Number(document.getElementById("C12").value);
    var C12_in2 = Number(document.getElementById("P12").value);
    var C21_in1 = Number(document.getElementById("C21").value);
    var C21_in2 = Number(document.getElementById("P21").value);
    var C23_in1 = Number(document.getElementById("C23").value);
    var C23_in2 = Number(document.getElementById("P23").value);
    var C30_in1 = Number(document.getElementById("C30").value);
    var C32_in1 = Number(document.getElementById("C32").value);
    var C32_in2 = Number(document.getElementById("P32").value);
    var C34_in1 = Number(document.getElementById("C34").value);
    var C34_in2 = Number(document.getElementById("P34").value);
    var C41_in1 = Number(document.getElementById("C41").value);
    var C41_in2 = Number(document.getElementById("P41").value);
    var C43_in1 = Number(document.getElementById("C43").value);
    var C43_in2 = Number(document.getElementById("P43").value);
    var C45_in1 = Number(document.getElementById("C45").value);
    var C45_in2 = Number(document.getElementById("P45").value);
    var C50_in1 = Number(document.getElementById("C50").value);
    var C52_in1 = Number(document.getElementById("C52").value);
    var C52_in2 = Number(document.getElementById("P52").value);
    var C54_in1 = Number(document.getElementById("C54").value);
    var C54_in2 = Number(document.getElementById("P54").value);
    var C56_in1 = Number(document.getElementById("C56").value);
    var C56_in2 = Number(document.getElementById("P56").value);
    var C10 = [1,0,C10_in1*ang, 0];
    var C12 = [1,2,C12_in1*nm,C12_in2*deg];
    var C21 = [2,1,C21_in1*nm,C21_in2*deg];
    var C23 = [2,3,C23_in1*nm,C23_in2*deg];
    var C30 = [3,0,C30_in1*nm,0];
    var C32 = [3,2,C32_in1*nm,C32_in2*deg]
    var C34 = [3,4,C34_in1*nm,C34_in2*deg]
    var C41 = [4,1,C41_in1*um,C41_in2*deg]
    var C43 = [4,3,C43_in1*um,C43_in2*deg]
    var C45 = [4,5,C45_in1*um,C45_in2*deg]
    var C50 = [5,0,C50_in1*mm,0]
    var C52 = [5,2,C52_in1*mm,C52_in2*deg]
    var C54 = [5,4,C54_in1*mm,C54_in2*deg]
    var C56 = [5,6,C56_in1*mm,C56_in2*deg]

    return math.matrix([C10,C12,C21,C23,C30,C32,C34,C41,C43,C45,C50,C52,C54,C56]);
}


function calculate(){

    ang = math.pow(10,-10);
    nm = math.pow(10,-9);
    um = math.pow(10,-6);
    mm = math.pow(10,-3);
    PI = math.pi;
    deg = PI/180;

    numPx = 256;
    lambda = 12.3986/math.sqrt((2*511+300)*300)*ang; // 300 keV


    //alpha meshgrid
    al_max = 70*math.pow(10,-3);
    al_vec = math.matrix(math.range(-al_max,al_max,(2*al_max)/(numPx)));
    al_vec.resize([256,1])

    alxx = math.multiply(math.ones(256,1),math.transpose(al_vec));
    alyy = math.transpose(alxx);


    alrr = math.sqrt(math.add(math.dotPow(alxx,2),math.dotPow(alyy,2)));
    alpp = math.atan2(alyy,alxx);

    obj_ap_r = 50 * math.pow(10,-3);
    obj_ap = alrr.map(function (value, index, matrix) {
        if(value < obj_ap_r)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    });

    trans = math.exp(  math.multiply(math.complex(0,-1),PI,.25, math.matrix(math.random([numPx,numPx])))  );
    
    aber = getAberrations();
    numAber = aber.size()[0];

    chi = math.zeros(numPx,numPx);

    for(var it = 0; it < numAber; it++)
    {
        //n = aber.subset(math.index(it,0));
        //m = aber.subset(math.index(it,1));
       // Cnm = aber.subset(math.index(it,2));
       // phinm = aber.subset(math.index(it,3));        
        chi = math.add(chi, math.dotMultiply(math.dotMultiply(math.cos(math.dotMultiply(aber.subset(math.index(it,1)),math.subtract(alpp,aber.subset(math.index(it,3))))),math.dotPow(alrr,aber.subset(math.index(it,0))+1)), aber.subset(math.index(it,2))/(aber.subset(math.index(it,0))+1) ));
        //chi = math.add(chi, math.dotMultiply(math.dotMultiply(math.cos(math.dotMultiply(m,math.subtract(alpp,phinm))),math.dotPow(alrr,n+1)), Cnm/(n+1) ));
    }
    chi0 = math.dotMultiply(2*PI/lambda, chi);
    expchi0 = math.dotMultiply(math.dotPow(math.E, math.dotMultiply(math.complex(0,-1),chi0) ), obj_ap);


  //   psi_p = math.matrix(fft2_wrap(expchi0.toArray()));
  //  psi_t = math.dotMultiply(trans,psi_p);
  //    ronch = math.dotMultiply(math.matrix(fft2_wrap(psi_t.toArray())),obj_ap); // multiply with obj apperture
  //  ronch = math.dotPow(math.abs(ronch),2);

  //  psi_t = math.dotMultiply(trans,math.matrix(fft2_wrap(expchi0.toArray())));    

   // ronch = math.dotMultiply(math.matrix(fft2_wrap(math.dotMultiply(trans,math.matrix(fft2_wrap(expchi0.toArray()))).toArray())),obj_ap); // multiply with obj apperture
    out_ronch = math.dotPow(math.abs(math.dotMultiply(math.matrix(fft2_wrap(math.dotMultiply(trans,math.matrix(fft2_wrap(expchi0.toArray()))).toArray())),obj_ap)),2);

    out_ronch = math.subtract(out_ronch, math.min(out_ronch));
    out_ronch = math.dotDivide(out_ronch,math.max(out_ronch)/255);
    out_ronch = math.round(out_ronch);
    out_ronch = out_ronch.toArray();


    out_phase_map = chi0.map(function (value, index, matrix) {
        if(value < PI/4)
        {
            return 1;            
        }
        else
        {
            return 0;
        }
    });

    rmax = math.dotDivide(1,math.dotMultiply(alrr,math.subtract(out_phase_map,1)));
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
    for(var it = 0; it < numPx; it++)
    {
        for(var jt = 0; jt < numPx; jt++)
        {
            var value = out_ronch[it][jt];
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


    canvas = document.getElementById("canvas2");
    ctx = canvas.getContext("2d");        
    canvas.width = numPx;
    canvas.height = numPx;
    for(var it = 0; it < numPx; it++)
    {
        for(var jt = 0; jt < numPx; jt++)
        {
            var value = out_phase_map[it][jt];
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

function randomize(){
    document.getElementById("C10").value = Math.random()*100;
    document.getElementById("C12").value = Math.random()*100;
    document.getElementById("P12").value = Math.random()*180;
    document.getElementById("C21").value = Math.random()*100;
    document.getElementById("P21").value = Math.random()*180;
    document.getElementById("C23").value = Math.random()*100;
    document.getElementById("P23").value = Math.random()*180;
    document.getElementById("C30").value = Math.random()*100;
    document.getElementById("C32").value = Math.random()*100;
    document.getElementById("P32").value = Math.random()*100;
    document.getElementById("C34").value = Math.random()*100;
    document.getElementById("P34").value = Math.random()*100;
    document.getElementById("C41").value = Math.random()*100;
    document.getElementById("P41").value = Math.random()*100;
    document.getElementById("C43").value = Math.random()*100;
    document.getElementById("P43").value = Math.random()*100;
    document.getElementById("C45").value = Math.random()*100;
    document.getElementById("P45").value = Math.random()*100;
    document.getElementById("C50").value = Math.random()*100;
    document.getElementById("C52").value = Math.random()*100;
    document.getElementById("P52").value = Math.random()*100;
    document.getElementById("C54").value = Math.random()*100;
    document.getElementById("P54").value = Math.random()*100;
    document.getElementById("C56").value = Math.random()*100;
    document.getElementById("P56").value = Math.random()*100;
    calculate();
}

function allZero(){
    document.getElementById("C10").value = 0;
    document.getElementById("C12").value = 0;
    document.getElementById("P12").value = 0;
    document.getElementById("C21").value = 0;
    document.getElementById("P21").value = 0;
    document.getElementById("C23").value = 0;
    document.getElementById("P23").value = 0;
    document.getElementById("C30").value = 0;
    document.getElementById("C32").value = 0;
    document.getElementById("P32").value = 0;
    document.getElementById("C34").value = 0;
    document.getElementById("P34").value = 0;
    document.getElementById("C41").value = 0;
    document.getElementById("P41").value = 0;
    document.getElementById("C43").value = 0;
    document.getElementById("P43").value = 0;
    document.getElementById("C45").value = 0;
    document.getElementById("P45").value = 0;
    document.getElementById("C50").value = 0;
    document.getElementById("C52").value = 0;
    document.getElementById("P52").value = 0;
    document.getElementById("C54").value = 0;
    document.getElementById("P54").value = 0;
    document.getElementById("C56").value = 0;
    document.getElementById("P56").value = 0;
}

function setC10(){
    document.getElementById("C10").value = 50;
}
function setC12(){
    document.getElementById("C12").value = 50;
}
function setC21(){
    document.getElementById("C21").value = 50;
}
function setC23(){
    document.getElementById("C23").value = 50;
}
function setC30(){
    document.getElementById("C30").value = 50;
}
function setC32(){
    document.getElementById("C32").value = 50;
}
function setC34(){
    document.getElementById("C34").value = 50;
}
function setC41(){
    document.getElementById("C41").value = 50;
}
function setC43(){
    document.getElementById("C43").value = 50;
}
function setC45(){
    document.getElementById("C45").value = 50;
}
function setC50(){
    document.getElementById("C50").value = 50;
}
function setC52(){
    document.getElementById("C52").value = 50;
}
function setC54(){
    document.getElementById("C54").value = 50;
}
function setC56(){
    document.getElementById("56").value = 50;
}

window.onload = calculate();

