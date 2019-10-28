
var ronchModule = (function() {
  var _scriptDir = typeof document !== 'undefined' && document.currentScript ? document.currentScript.src : undefined;
  return (
function(ronchModule) {
  ronchModule = ronchModule || {};

var b;b||(b=typeof ronchModule !== 'undefined' ? ronchModule : {});var h={},l;for(l in b)b.hasOwnProperty(l)&&(h[l]=b[l]);b.arguments=[];b.thisProgram="./this.program";b.quit=function(a,c){throw c;};b.preRun=[];b.postRun=[];var m=!1,q=!1,r=!1,t=!1;m="object"===typeof window;q="function"===typeof importScripts;r="object"===typeof process&&"function"===typeof require&&!m&&!q;t=!m&&!r&&!q;var u="";
if(r){u=__dirname+"/";var v,w;b.read=function(a,c){v||(v=require("fs"));w||(w=require("path"));a=w.normalize(a);a=v.readFileSync(a);return c?a:a.toString()};b.readBinary=function(a){a=b.read(a,!0);a.buffer||(a=new Uint8Array(a));assert(a.buffer);return a};1<process.argv.length&&(b.thisProgram=process.argv[1].replace(/\\/g,"/"));b.arguments=process.argv.slice(2);process.on("uncaughtException",function(a){if(!(a instanceof x))throw a;});process.on("unhandledRejection",z);b.quit=function(a){process.exit(a)};
b.inspect=function(){return"[Emscripten Module object]"}}else if(t)"undefined"!=typeof read&&(b.read=function(a){return read(a)}),b.readBinary=function(a){if("function"===typeof readbuffer)return new Uint8Array(readbuffer(a));a=read(a,"binary");assert("object"===typeof a);return a},"undefined"!=typeof scriptArgs?b.arguments=scriptArgs:"undefined"!=typeof arguments&&(b.arguments=arguments),"function"===typeof quit&&(b.quit=function(a){quit(a)});else if(m||q)q?u=self.location.href:document.currentScript&&
(u=document.currentScript.src),_scriptDir&&(u=_scriptDir),0!==u.indexOf("blob:")?u=u.substr(0,u.lastIndexOf("/")+1):u="",b.read=function(a){var c=new XMLHttpRequest;c.open("GET",a,!1);c.send(null);return c.responseText},q&&(b.readBinary=function(a){var c=new XMLHttpRequest;c.open("GET",a,!1);c.responseType="arraybuffer";c.send(null);return new Uint8Array(c.response)}),b.readAsync=function(a,c,d){var e=new XMLHttpRequest;e.open("GET",a,!0);e.responseType="arraybuffer";e.onload=function(){200==e.status||
0==e.status&&e.response?c(e.response):d()};e.onerror=d;e.send(null)},b.setWindowTitle=function(a){document.title=a};var A=b.print||("undefined"!==typeof console?console.log.bind(console):"undefined"!==typeof print?print:null),B=b.printErr||("undefined"!==typeof printErr?printErr:"undefined"!==typeof console&&console.warn.bind(console)||A);for(l in h)h.hasOwnProperty(l)&&(b[l]=h[l]);h=void 0;var aa={"f64-rem":function(a,c){return a%c},"debugger":function(){debugger}};
"object"!==typeof WebAssembly&&B("no native wasm support detected");var C,D=!1;function assert(a,c){a||z("Assertion failed: "+c)}function ba(a){var c=b["_"+a];assert(c,"Cannot call unknown function "+a+", make sure it is exported");return c}var E="undefined"!==typeof TextDecoder?new TextDecoder("utf8"):void 0;
function F(a,c,d){var e=c+d;for(d=c;a[d]&&!(d>=e);)++d;if(16<d-c&&a.subarray&&E)return E.decode(a.subarray(c,d));for(e="";c<d;){var f=a[c++];if(f&128){var n=a[c++]&63;if(192==(f&224))e+=String.fromCharCode((f&31)<<6|n);else{var p=a[c++]&63;f=224==(f&240)?(f&15)<<12|n<<6|p:(f&7)<<18|n<<12|p<<6|a[c++]&63;65536>f?e+=String.fromCharCode(f):(f-=65536,e+=String.fromCharCode(55296|f>>10,56320|f&1023))}}else e+=String.fromCharCode(f)}return e}"undefined"!==typeof TextDecoder&&new TextDecoder("utf-16le");
var buffer,G,H,I,J=b.TOTAL_MEMORY||16777216;5242880>J&&B("TOTAL_MEMORY should be larger than TOTAL_STACK, was "+J+"! (TOTAL_STACK=5242880)");b.buffer?buffer=b.buffer:"object"===typeof WebAssembly&&"function"===typeof WebAssembly.Memory?(C=new WebAssembly.Memory({initial:J/65536,maximum:J/65536}),buffer=C.buffer):buffer=new ArrayBuffer(J);b.HEAP8=G=new Int8Array(buffer);b.HEAP16=new Int16Array(buffer);b.HEAP32=I=new Int32Array(buffer);b.HEAPU8=H=new Uint8Array(buffer);b.HEAPU16=new Uint16Array(buffer);
b.HEAPU32=new Uint32Array(buffer);b.HEAPF32=new Float32Array(buffer);b.HEAPF64=new Float64Array(buffer);I[1808]=5250144;function K(a){for(;0<a.length;){var c=a.shift();if("function"==typeof c)c();else{var d=c.C;"number"===typeof d?void 0===c.v?b.dynCall_v(d):b.dynCall_vi(d,c.v):d(void 0===c.v?null:c.v)}}}var L=[],ca=[],da=[],M=[],P=!1;function ea(){var a=b.preRun.shift();L.unshift(a)}var Q=0,R=null,S=null;b.preloadedImages={};b.preloadedAudios={};
function T(){var a=U;return String.prototype.startsWith?a.startsWith("data:application/octet-stream;base64,"):0===a.indexOf("data:application/octet-stream;base64,")}var U="index.wasm";if(!T()){var V=U;U=b.locateFile?b.locateFile(V,u):u+V}function W(){try{if(b.wasmBinary)return new Uint8Array(b.wasmBinary);if(b.readBinary)return b.readBinary(U);throw"both async and sync fetching of the wasm failed";}catch(a){z(a)}}
function fa(){return b.wasmBinary||!m&&!q||"function"!==typeof fetch?new Promise(function(a){a(W())}):fetch(U,{credentials:"same-origin"}).then(function(a){if(!a.ok)throw"failed to load wasm binary file at '"+U+"'";return a.arrayBuffer()}).catch(function(){return W()})}
function ha(a){function c(a){b.asm=a.exports;Q--;b.monitorRunDependencies&&b.monitorRunDependencies(Q);0==Q&&(null!==R&&(clearInterval(R),R=null),S&&(a=S,S=null,a()))}function d(a){c(a.instance)}function e(a){fa().then(function(a){return WebAssembly.instantiate(a,f)}).then(a,function(a){B("failed to asynchronously prepare wasm: "+a);z(a)})}var f={env:a,global:{NaN:NaN,Infinity:Infinity},"global.Math":Math,asm2wasm:aa};Q++;b.monitorRunDependencies&&b.monitorRunDependencies(Q);if(b.instantiateWasm)try{return b.instantiateWasm(f,
c)}catch(n){return B("Module.instantiateWasm callback failed with error: "+n),!1}b.wasmBinary||"function"!==typeof WebAssembly.instantiateStreaming||T()||"function"!==typeof fetch?e(d):WebAssembly.instantiateStreaming(fetch(U,{credentials:"same-origin"}),f).then(d,function(a){B("wasm streaming compile failed: "+a);B("falling back to ArrayBuffer instantiation");e(d)});return{}}
b.asm=function(a,c){c.memory=C;c.table=new WebAssembly.Table({initial:36,maximum:36,element:"anyfunc"});c.__memory_base=1024;c.__table_base=0;return ha(c)};var ia=[null,[],[]],X=0;function Y(){X+=4;return I[X-4>>2]}var ja={};function ka(){z("OOM")}
var la=b.asm({},{b:z,d:function(a){b.___errno_location&&(I[b.___errno_location()>>2]=a);return a},k:function(a,c){X=c;try{return ja.B(),Y(),Y(),Y(),Y(),0}catch(d){return"undefined"!==typeof FS&&d instanceof FS.w||z(d),-d.A}},c:function(a,c){X=c;try{var d=Y(),e=Y(),f=Y();for(c=a=0;c<f;c++){for(var n=I[e+8*c>>2],p=I[e+(8*c+4)>>2],k=0;k<p;k++){var y=H[n+k],N=ia[d];0===y||10===y?((1===d?A:B)(F(N,0)),N.length=0):N.push(y)}a+=p}return a}catch(O){return"undefined"!==typeof FS&&O instanceof FS.w||z(O),-O.A}},
j:function(a,c){X=c;try{return ja.B(),0}catch(d){return"undefined"!==typeof FS&&d instanceof FS.w||z(d),-d.A}},i:function(){b.abort()},h:function(){return G.length},g:function(a,c,d){H.set(H.subarray(c,c+d),a)},f:function(a){ka(a)},e:function(a){var c=Date.now()/1E3|0;a&&(I[a>>2]=c);return c},l:ka,a:7232},buffer);b.asm=la;b.___errno_location=function(){return b.asm.m.apply(null,arguments)};b._calcRonch=function(){return b.asm.n.apply(null,arguments)};b._free=function(){return b.asm.o.apply(null,arguments)};
b._malloc=function(){return b.asm.p.apply(null,arguments)};var ma=b.stackAlloc=function(){return b.asm.s.apply(null,arguments)},na=b.stackRestore=function(){return b.asm.t.apply(null,arguments)},oa=b.stackSave=function(){return b.asm.u.apply(null,arguments)};b.dynCall_v=function(){return b.asm.q.apply(null,arguments)};b.dynCall_vi=function(){return b.asm.r.apply(null,arguments)};b.asm=la;
b.ccall=function(a,c,d,e){var f={string:function(a){var c=0;if(null!==a&&void 0!==a&&0!==a){var d=(a.length<<2)+1;c=ma(d);var e=c;if(0<d){d=e+d-1;for(var f=0;f<a.length;++f){var g=a.charCodeAt(f);if(55296<=g&&57343>=g){var k=a.charCodeAt(++f);g=65536+((g&1023)<<10)|k&1023}if(127>=g){if(e>=d)break;H[e++]=g}else{if(2047>=g){if(e+1>=d)break;H[e++]=192|g>>6}else{if(65535>=g){if(e+2>=d)break;H[e++]=224|g>>12}else{if(e+3>=d)break;H[e++]=240|g>>18;H[e++]=128|g>>12&63}H[e++]=128|g>>6&63}H[e++]=128|g&63}}H[e]=
0}}return c},array:function(a){var c=ma(a.length);G.set(a,c);return c}},n=ba(a),p=[];a=0;if(e)for(var k=0;k<e.length;k++){var y=f[d[k]];y?(0===a&&(a=oa()),p[k]=y(e[k])):p[k]=e[k]}d=n.apply(null,p);d=function(a){return"string"===c?a?F(H,a,void 0):"":"boolean"===c?!!a:a}(d);0!==a&&na(a);return d};b.then=function(a){if(b.calledRun)a(b);else{var c=b.onRuntimeInitialized;b.onRuntimeInitialized=function(){c&&c();a(b)}}return b};
function x(a){this.name="ExitStatus";this.message="Program terminated with exit("+a+")";this.status=a}x.prototype=Error();x.prototype.constructor=x;S=function pa(){b.calledRun||Z();b.calledRun||(S=pa)};
function Z(){function a(){if(!b.calledRun&&(b.calledRun=!0,!D)){P||(P=!0,K(ca));K(da);if(b.onRuntimeInitialized)b.onRuntimeInitialized();if(b.postRun)for("function"==typeof b.postRun&&(b.postRun=[b.postRun]);b.postRun.length;){var a=b.postRun.shift();M.unshift(a)}K(M)}}if(!(0<Q)){if(b.preRun)for("function"==typeof b.preRun&&(b.preRun=[b.preRun]);b.preRun.length;)ea();K(L);0<Q||b.calledRun||(b.setStatus?(b.setStatus("Running..."),setTimeout(function(){setTimeout(function(){b.setStatus("")},1);a()},
1)):a())}}b.run=Z;function z(a){if(b.onAbort)b.onAbort(a);void 0!==a?(A(a),B(a),a=JSON.stringify(a)):a="";D=!0;throw"abort("+a+"). Build with -s ASSERTIONS=1 for more info.";}b.abort=z;if(b.preInit)for("function"==typeof b.preInit&&(b.preInit=[b.preInit]);0<b.preInit.length;)b.preInit.pop()();b.noExitRuntime=!0;Z();


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
    