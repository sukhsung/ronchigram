
var ronchModule = (function() {
  var _scriptDir = typeof document !== 'undefined' && document.currentScript ? document.currentScript.src : undefined;
  return (
function(ronchModule) {
  ronchModule = ronchModule || {};

var b;b||(b=typeof ronchModule !== 'undefined' ? ronchModule : {});var k={},l;for(l in b)b.hasOwnProperty(l)&&(k[l]=b[l]);var p=!1,q=!1,r=!1,aa=!1,t=!1;p="object"===typeof window;q="function"===typeof importScripts;r=(aa="object"===typeof process&&"object"===typeof process.versions&&"string"===typeof process.versions.node)&&!p&&!q;t=!p&&!r&&!q;var u="",v,w;
if(r){u=__dirname+"/";var x,y;v=function(a,c){x||(x=require("fs"));y||(y=require("path"));a=y.normalize(a);a=x.readFileSync(a);return c?a:a.toString()};w=function(a){a=v(a,!0);a.buffer||(a=new Uint8Array(a));assert(a.buffer);return a};1<process.argv.length&&process.argv[1].replace(/\\/g,"/");process.argv.slice(2);process.on("uncaughtException",function(a){throw a;});process.on("unhandledRejection",A);b.inspect=function(){return"[Emscripten Module object]"}}else if(t)"undefined"!=typeof read&&(v=function(a){return read(a)}),
w=function(a){if("function"===typeof readbuffer)return new Uint8Array(readbuffer(a));a=read(a,"binary");assert("object"===typeof a);return a},"undefined"!==typeof print&&("undefined"===typeof console&&(console={}),console.log=print,console.warn=console.error="undefined"!==typeof printErr?printErr:print);else if(p||q)q?u=self.location.href:document.currentScript&&(u=document.currentScript.src),_scriptDir&&(u=_scriptDir),0!==u.indexOf("blob:")?u=u.substr(0,u.lastIndexOf("/")+1):u="",v=function(a){var c=
new XMLHttpRequest;c.open("GET",a,!1);c.send(null);return c.responseText},q&&(w=function(a){var c=new XMLHttpRequest;c.open("GET",a,!1);c.responseType="arraybuffer";c.send(null);return new Uint8Array(c.response)});var B=b.print||console.log.bind(console),C=b.printErr||console.warn.bind(console);for(l in k)k.hasOwnProperty(l)&&(b[l]=k[l]);k=null;var ba={"f64-rem":function(a,c){return a%c},"debugger":function(){}},D;b.wasmBinary&&(D=b.wasmBinary);"object"!==typeof WebAssembly&&C("no native wasm support detected");
var E,F=!1;function assert(a,c){a||A("Assertion failed: "+c)}function ca(a){var c=b["_"+a];assert(c,"Cannot call unknown function "+a+", make sure it is exported");return c}var H="undefined"!==typeof TextDecoder?new TextDecoder("utf8"):void 0;
function I(a,c,d){var f=c+d;for(d=c;a[d]&&!(d>=f);)++d;if(16<d-c&&a.subarray&&H)return H.decode(a.subarray(c,d));for(f="";c<d;){var e=a[c++];if(e&128){var m=a[c++]&63;if(192==(e&224))f+=String.fromCharCode((e&31)<<6|m);else{var n=a[c++]&63;e=224==(e&240)?(e&15)<<12|m<<6|n:(e&7)<<18|m<<12|n<<6|a[c++]&63;65536>e?f+=String.fromCharCode(e):(e-=65536,f+=String.fromCharCode(55296|e>>10,56320|e&1023))}}else f+=String.fromCharCode(e)}return f}"undefined"!==typeof TextDecoder&&new TextDecoder("utf-16le");
function J(a){0<a%65536&&(a+=65536-a%65536);return a}var buffer,K,L,M;function N(a){buffer=a;b.HEAP8=K=new Int8Array(a);b.HEAP16=new Int16Array(a);b.HEAP32=M=new Int32Array(a);b.HEAPU8=L=new Uint8Array(a);b.HEAPU16=new Uint16Array(a);b.HEAPU32=new Uint32Array(a);b.HEAPF32=new Float32Array(a);b.HEAPF64=new Float64Array(a)}var O=b.TOTAL_MEMORY||16777216;b.wasmMemory?E=b.wasmMemory:E=new WebAssembly.Memory({initial:O/65536});E&&(buffer=E.buffer);O=buffer.byteLength;N(buffer);M[3780]=5258032;
function P(a){for(;0<a.length;){var c=a.shift();if("function"==typeof c)c();else{var d=c.F;"number"===typeof d?void 0===c.A?b.dynCall_v(d):b.dynCall_vi(d,c.A):d(void 0===c.A?null:c.A)}}}var R=[],da=[],ea=[],fa=[];function ha(){var a=b.preRun.shift();R.unshift(a)}var S=0,T=null,U=null;b.preloadedImages={};b.preloadedAudios={};function ia(){var a=V;return String.prototype.startsWith?a.startsWith("data:application/octet-stream;base64,"):0===a.indexOf("data:application/octet-stream;base64,")}var V="index.wasm";
if(!ia()){var ja=V;V=b.locateFile?b.locateFile(ja,u):u+ja}function ka(){try{if(D)return new Uint8Array(D);if(w)return w(V);throw"both async and sync fetching of the wasm failed";}catch(a){A(a)}}function la(){return D||!p&&!q||"function"!==typeof fetch?new Promise(function(a){a(ka())}):fetch(V,{credentials:"same-origin"}).then(function(a){if(!a.ok)throw"failed to load wasm binary file at '"+V+"'";return a.arrayBuffer()}).catch(function(){return ka()})}
function ma(a){function c(a){b.asm=a.exports;S--;b.monitorRunDependencies&&b.monitorRunDependencies(S);0==S&&(null!==T&&(clearInterval(T),T=null),U&&(a=U,U=null,a()))}function d(a){c(a.instance)}function f(a){return la().then(function(a){return WebAssembly.instantiate(a,e)}).then(a,function(a){C("failed to asynchronously prepare wasm: "+a);A(a)})}var e={env:a,global:{NaN:NaN,Infinity:Infinity},"global.Math":Math,asm2wasm:ba};S++;b.monitorRunDependencies&&b.monitorRunDependencies(S);if(b.instantiateWasm)try{return b.instantiateWasm(e,
c)}catch(m){return C("Module.instantiateWasm callback failed with error: "+m),!1}(function(){if(D||"function"!==typeof WebAssembly.instantiateStreaming||ia()||"function"!==typeof fetch)return f(d);fetch(V,{credentials:"same-origin"}).then(function(a){return WebAssembly.instantiateStreaming(a,e).then(d,function(a){C("wasm streaming compile failed: "+a);C("falling back to ArrayBuffer instantiation");f(d)})})})();return{}}
b.asm=function(a,c){c.memory=E;c.table=new WebAssembly.Table({initial:314,maximum:314,element:"anyfunc"});c.__memory_base=1024;c.__table_base=0;return ma(c)};var na=[null,[],[]],W=0;function X(){W+=4;return M[W-4>>2]}var oa={};function pa(){return K.length}
var qa=b.asm({},{b:A,n:function(){F=!0;throw"Pure virtual function called!";},d:function(a){b.___errno_location&&(M[b.___errno_location()>>2]=a);return a},i:function(a,c){W=c;try{return oa.D(),X(),X(),X(),X(),0}catch(d){return"undefined"!==typeof FS&&d instanceof FS.B||A(d),-d.C}},c:function(a,c){W=c;try{var d=X(),f=X(),e=X();for(c=a=0;c<e;c++){for(var m=M[f+8*c>>2],n=M[f+(8*c+4)>>2],g=0;g<n;g++){var z=L[m+g],G=na[d];0===z||10===z?((1===d?B:C)(I(G,0)),G.length=0):G.push(z)}a+=n}return a}catch(Q){return"undefined"!==
typeof FS&&Q instanceof FS.B||A(Q),-Q.C}},h:function(a,c){W=c;try{return oa.D(),0}catch(d){return"undefined"!==typeof FS&&d instanceof FS.B||A(d),-d.C}},g:function(){b.abort()},f:pa,e:function(a,c,d){L.set(L.subarray(c,c+d),a)},m:function(a){if(2147418112<a)return!1;for(var c=Math.max(pa(),16777216);c<a;)536870912>=c?c=J(2*c):c=Math.min(J((3*c+2147483648)/4),2147418112);a:{try{E.grow(c-buffer.byteLength+65535>>16);N(E.buffer);var d=1;break a}catch(f){}d=void 0}return d?!0:!1},l:function(){A("trap!")},
k:function(a){var c=Date.now()/1E3|0;a&&(M[a>>2]=c);return c},j:function(){A("OOM")},a:15120},buffer);b.asm=qa;b.___errno_location=function(){return b.asm.o.apply(null,arguments)};b._calcRonch=function(){return b.asm.p.apply(null,arguments)};b._free=function(){return b.asm.q.apply(null,arguments)};b._malloc=function(){return b.asm.r.apply(null,arguments)};
var ra=b.stackAlloc=function(){return b.asm.u.apply(null,arguments)},sa=b.stackRestore=function(){return b.asm.v.apply(null,arguments)},ta=b.stackSave=function(){return b.asm.w.apply(null,arguments)};b.dynCall_v=function(){return b.asm.s.apply(null,arguments)};b.dynCall_vi=function(){return b.asm.t.apply(null,arguments)};b.asm=qa;
b.ccall=function(a,c,d,f){var e={string:function(a){var c=0;if(null!==a&&void 0!==a&&0!==a){var d=(a.length<<2)+1;c=ra(d);var e=c,f=L;if(0<d){d=e+d-1;for(var g=0;g<a.length;++g){var h=a.charCodeAt(g);if(55296<=h&&57343>=h){var G=a.charCodeAt(++g);h=65536+((h&1023)<<10)|G&1023}if(127>=h){if(e>=d)break;f[e++]=h}else{if(2047>=h){if(e+1>=d)break;f[e++]=192|h>>6}else{if(65535>=h){if(e+2>=d)break;f[e++]=224|h>>12}else{if(e+3>=d)break;f[e++]=240|h>>18;f[e++]=128|h>>12&63}f[e++]=128|h>>6&63}f[e++]=128|h&
63}}f[e]=0}}return c},array:function(a){var c=ra(a.length);K.set(a,c);return c}},m=ca(a),n=[];a=0;if(f)for(var g=0;g<f.length;g++){var z=e[d[g]];z?(0===a&&(a=ta()),n[g]=z(f[g])):n[g]=f[g]}d=m.apply(null,n);d=function(a){return"string"===c?a?I(L,a,void 0):"":"boolean"===c?!!a:a}(d);0!==a&&sa(a);return d};var Y;b.then=function(a){if(Y)a(b);else{var c=b.onRuntimeInitialized;b.onRuntimeInitialized=function(){c&&c();a(b)}}return b};U=function ua(){Y||Z();Y||(U=ua)};
function Z(){function a(){if(!Y&&(Y=!0,!F)){P(da);P(ea);if(b.onRuntimeInitialized)b.onRuntimeInitialized();if(b.postRun)for("function"==typeof b.postRun&&(b.postRun=[b.postRun]);b.postRun.length;){var a=b.postRun.shift();fa.unshift(a)}P(fa)}}if(!(0<S)){if(b.preRun)for("function"==typeof b.preRun&&(b.preRun=[b.preRun]);b.preRun.length;)ha();P(R);0<S||(b.setStatus?(b.setStatus("Running..."),setTimeout(function(){setTimeout(function(){b.setStatus("")},1);a()},1)):a())}}b.run=Z;
function A(a){if(b.onAbort)b.onAbort(a);B(a);C(a);F=!0;throw"abort("+a+"). Build with -s ASSERTIONS=1 for more info.";}b.abort=A;if(b.preInit)for("function"==typeof b.preInit&&(b.preInit=[b.preInit]);0<b.preInit.length;)b.preInit.pop()();Z();


  return ronchModule
}
);
})();
if (typeof exports === 'object' && typeof module === 'object')
      module.exports = ronchModule;
    else if (typeof define === 'function' && define['amd'])
      define([], function() { return ronchModule; });
    else if (typeof exports === 'object')
      exports["ronchModule"] = ronchModule;
    