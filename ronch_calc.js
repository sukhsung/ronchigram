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
          var a = -2*math.pi;
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

            var C10 = [1,0,C10_in1*ang, 0];
            var C12 = [1,2,C12_in1*nm,C12_in2*deg];
            var C21 = [2,1,C21_in1*nm,C21_in2*deg];
            var C23 = [2,3,C23_in1*nm,C23_in2*deg];
            var C30 = [3,0,C30_in1*nm,0];
            return math.matrix([C10,C12,C21,C23,C30]);

            //return math.matrix([ [1,0,C10,0], [1,2,C12,Ph12] ]);
        }


        function calculate(){
            numPx = 256;
            kev = 300;
            lambda = 12.3986/math.sqrt((2*511+kev)*kev)*math.pow(10,-10);

            ang = math.pow(10,-10);
            nm = math.pow(10,-9);
            deg = math.pi/180;

            //alpha meshgrid
            al_max = 70*math.pow(10,-3);
            al_step = (2*al_max)/(numPx);
            al_vec = math.matrix(math.range(-al_max,al_max,al_step));
            al_vec.resize([256,1])
            ones = math.ones(256,1);

            alxx = math.multiply(ones,math.transpose(al_vec));
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

            C10 = 50*ang;
            C12 = 15*nm;//30*math.random()*nm;
            Ph12 = 0*deg;//180*(2*math.random()-1)/2;


            trans = math.exp(  math.multiply(math.complex(0,-1),math.pi,.25, math.matrix(math.random([numPx,numPx])))  );
            aber = getAberrations();
            //aber = math.matrix([ [1,0,C10,0], [1,2,C12,Ph12] ]);

            chi = math.zeros(numPx,numPx);

            for(var it = 0; it < aber.size()[0]; it++)
            {
                n = aber.subset(math.index(it,0));
                m = aber.subset(math.index(it,1));
                Cnm = aber.subset(math.index(it,2));
                phinm = aber.subset(math.index(it,3));
                chi = math.add(chi, math.dotMultiply(math.dotMultiply(math.cos(math.dotMultiply(m,math.subtract(alpp,phinm))),math.dotPow(alrr,n+1)), Cnm/(n+1) ));
            }
            chi0 = math.dotMultiply(2*math.pi/lambda, chi);
            expchi0 = math.dotMultiply(math.dotPow(math.E, math.dotMultiply(math.complex(0,-1),chi0) ), obj_ap);

            psi_p = math.matrix(fft2_wrap(expchi0.toArray()));
            psi_t = math.dotMultiply(trans,psi_p);

            ronch = math.dotMultiply(math.matrix(fft2_wrap(psi_t.toArray())),obj_ap); // multiply with obj apperture
            ronch = math.dotPow(math.abs(ronch),2);


            out_ronch = math.abs(ronch);
            out_ronch = math.subtract(out_ronch, math.min(out_ronch));
            out_ronch = math.dotDivide(out_ronch,math.max(out_ronch)/255);
            out_ronch = math.round(out_ronch);
            out_ronch = out_ronch.toArray();


            out_phase_map = chi0.map(function (value, index, matrix) {
                if(value < math.pi/4)
                {
                    return 1;            
                }
                else
                {
                    return 0;
                }
            });

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

            calculate();
        }

calculate();

