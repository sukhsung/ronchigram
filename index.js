
var ronchModule = (() => {
  var _scriptDir = typeof document !== 'undefined' && document.currentScript ? document.currentScript.src : undefined;
  if (typeof __filename !== 'undefined') _scriptDir = _scriptDir || __filename;
  return (
function(ronchModule) {
  ronchModule = ronchModule || {};


var b;b||(b=typeof ronchModule !== 'undefined' ? ronchModule : {});var m,n;b.ready=new Promise(function(a,c){m=a;n=c});var v=Object.assign({},b),w="object"==typeof window,x="function"==typeof importScripts,z="",A,B,C,fs,D,E;
if("object"==typeof process&&"object"==typeof process.versions&&"string"==typeof process.versions.node)z=x?require("path").dirname(z)+"/":__dirname+"/",E=()=>{D||(fs=require("fs"),D=require("path"))},A=function(a,c){E();a=D.normalize(a);return fs.readFileSync(a,c?void 0:"utf8")},C=a=>{a=A(a,!0);a.buffer||(a=new Uint8Array(a));return a},B=(a,c,f)=>{E();a=D.normalize(a);fs.readFile(a,function(d,e){d?f(d):c(e.buffer)})},1<process.argv.length&&process.argv[1].replace(/\\/g,"/"),process.argv.slice(2),
process.on("uncaughtException",function(a){throw a;}),process.on("unhandledRejection",function(a){throw a;}),b.inspect=function(){return"[Emscripten Module object]"};else if(w||x)x?z=self.location.href:"undefined"!=typeof document&&document.currentScript&&(z=document.currentScript.src),_scriptDir&&(z=_scriptDir),0!==z.indexOf("blob:")?z=z.substr(0,z.replace(/[?#].*/,"").lastIndexOf("/")+1):z="",A=a=>{var c=new XMLHttpRequest;c.open("GET",a,!1);c.send(null);return c.responseText},x&&(C=a=>{var c=new XMLHttpRequest;
c.open("GET",a,!1);c.responseType="arraybuffer";c.send(null);return new Uint8Array(c.response)}),B=(a,c,f)=>{var d=new XMLHttpRequest;d.open("GET",a,!0);d.responseType="arraybuffer";d.onload=()=>{200==d.status||0==d.status&&d.response?c(d.response):f()};d.onerror=f;d.send(null)};var aa=b.print||console.log.bind(console),F=b.printErr||console.warn.bind(console);Object.assign(b,v);v=null;var G;b.wasmBinary&&(G=b.wasmBinary);var noExitRuntime=b.noExitRuntime||!0;"object"!=typeof WebAssembly&&H("no native wasm support detected");
var J,K=!1,L="undefined"!=typeof TextDecoder?new TextDecoder("utf8"):void 0;
function M(a,c){for(var f=c+NaN,d=c;a[d]&&!(d>=f);)++d;if(16<d-c&&a.buffer&&L)return L.decode(a.subarray(c,d));for(f="";c<d;){var e=a[c++];if(e&128){var g=a[c++]&63;if(192==(e&224))f+=String.fromCharCode((e&31)<<6|g);else{var p=a[c++]&63;e=224==(e&240)?(e&15)<<12|g<<6|p:(e&7)<<18|g<<12|p<<6|a[c++]&63;65536>e?f+=String.fromCharCode(e):(e-=65536,f+=String.fromCharCode(55296|e>>10,56320|e&1023))}}else f+=String.fromCharCode(e)}return f}var N,O,P,Q;
function ba(){var a=J.buffer;N=a;b.HEAP8=O=new Int8Array(a);b.HEAP16=new Int16Array(a);b.HEAP32=new Int32Array(a);b.HEAPU8=P=new Uint8Array(a);b.HEAPU16=new Uint16Array(a);b.HEAPU32=Q=new Uint32Array(a);b.HEAPF32=new Float32Array(a);b.HEAPF64=new Float64Array(a)}var ca,da=[],ea=[],fa=[];function ha(){var a=b.preRun.shift();da.unshift(a)}var R=0,S=null,T=null;
function H(a){if(b.onAbort)b.onAbort(a);a="Aborted("+a+")";F(a);K=!0;a=new WebAssembly.RuntimeError(a+". Build with -sASSERTIONS for more info.");n(a);throw a;}function ia(){return U.startsWith("data:application/octet-stream;base64,")}var U;U="index.wasm";if(!ia()){var ja=U;U=b.locateFile?b.locateFile(ja,z):z+ja}function ka(){var a=U;try{if(a==U&&G)return new Uint8Array(G);if(C)return C(a);throw"both async and sync fetching of the wasm failed";}catch(c){H(c)}}
function la(){if(!G&&(w||x)){if("function"==typeof fetch&&!U.startsWith("file://"))return fetch(U,{credentials:"same-origin"}).then(function(a){if(!a.ok)throw"failed to load wasm binary file at '"+U+"'";return a.arrayBuffer()}).catch(function(){return ka()});if(B)return new Promise(function(a,c){B(U,function(f){a(new Uint8Array(f))},c)})}return Promise.resolve().then(function(){return ka()})}
function V(a){for(;0<a.length;){var c=a.shift();if("function"==typeof c)c(b);else{var f=c.u;"number"==typeof f?void 0===c.s?ma(f)():ma(f)(c.s):f(void 0===c.s?null:c.s)}}}var W=[];function ma(a){var c=W[a];c||(a>=W.length&&(W.length=a+1),W[a]=c=ca.get(a));return c}
var na=[null,[],[]],oa={f:function(){return Date.now()},a:function(){H("")},g:function(a,c,f){P.copyWithin(a,c,c+f)},d:function(a){var c=P.length;a>>>=0;if(2147483648<a)return!1;for(var f=1;4>=f;f*=2){var d=c*(1+.2/f);d=Math.min(d,a+100663296);var e=Math;d=Math.max(a,d);e=e.min.call(e,2147483648,d+(65536-d%65536)%65536);a:{try{J.grow(e-N.byteLength+65535>>>16);ba();var g=1;break a}catch(p){}g=void 0}if(g)return!0}return!1},e:function(){return 52},c:function(){return 70},b:function(a,c,f,d){for(var e=
0,g=0;g<f;g++){var p=Q[c>>2],q=Q[c+4>>2];c+=8;for(var y=0;y<q;y++){var h=P[p+y],r=na[a];0===h||10===h?((1===a?aa:F)(M(r,0)),r.length=0):r.push(h)}e+=q}Q[d>>2]=e;return 0}};
(function(){function a(e){b.asm=e.exports;J=b.asm.h;ba();ca=b.asm.m;ea.unshift(b.asm.i);R--;b.monitorRunDependencies&&b.monitorRunDependencies(R);0==R&&(null!==S&&(clearInterval(S),S=null),T&&(e=T,T=null,e()))}function c(e){a(e.instance)}function f(e){return la().then(function(g){return WebAssembly.instantiate(g,d)}).then(function(g){return g}).then(e,function(g){F("failed to asynchronously prepare wasm: "+g);H(g)})}var d={a:oa};R++;b.monitorRunDependencies&&b.monitorRunDependencies(R);if(b.instantiateWasm)try{return b.instantiateWasm(d,
a)}catch(e){return F("Module.instantiateWasm callback failed with error: "+e),!1}(function(){return G||"function"!=typeof WebAssembly.instantiateStreaming||ia()||U.startsWith("file://")||"function"!=typeof fetch?f(c):fetch(U,{credentials:"same-origin"}).then(function(e){return WebAssembly.instantiateStreaming(e,d).then(c,function(g){F("wasm streaming compile failed: "+g);F("falling back to ArrayBuffer instantiation");return f(c)})})})().catch(n);return{}})();
b.___wasm_call_ctors=function(){return(b.___wasm_call_ctors=b.asm.i).apply(null,arguments)};b._calcRonch=function(){return(b._calcRonch=b.asm.j).apply(null,arguments)};b._malloc=function(){return(b._malloc=b.asm.k).apply(null,arguments)};b._free=function(){return(b._free=b.asm.l).apply(null,arguments)};
var pa=b.stackSave=function(){return(pa=b.stackSave=b.asm.n).apply(null,arguments)},qa=b.stackRestore=function(){return(qa=b.stackRestore=b.asm.o).apply(null,arguments)},X=b.stackAlloc=function(){return(X=b.stackAlloc=b.asm.p).apply(null,arguments)};
b.ccall=function(a,c,f,d){var e={string:function(h){var r=0;if(null!==h&&void 0!==h&&0!==h){var t=(h.length<<2)+1;r=X(t);var l=r,u=P;if(0<t){t=l+t-1;for(var I=0;I<h.length;++I){var k=h.charCodeAt(I);if(55296<=k&&57343>=k){var ra=h.charCodeAt(++I);k=65536+((k&1023)<<10)|ra&1023}if(127>=k){if(l>=t)break;u[l++]=k}else{if(2047>=k){if(l+1>=t)break;u[l++]=192|k>>6}else{if(65535>=k){if(l+2>=t)break;u[l++]=224|k>>12}else{if(l+3>=t)break;u[l++]=240|k>>18;u[l++]=128|k>>12&63}u[l++]=128|k>>6&63}u[l++]=128|k&
63}}u[l]=0}}return r},array:function(h){var r=X(h.length);O.set(h,r);return r}};a=b["_"+a];var g=[],p=0;if(d)for(var q=0;q<d.length;q++){var y=e[f[q]];y?(0===p&&(p=pa()),g[q]=y(d[q])):g[q]=d[q]}f=a.apply(null,g);return f=function(h){0!==p&&qa(p);return"string"===c?h?M(P,h):"":"boolean"===c?!!h:h}(f)};var Y;T=function sa(){Y||Z();Y||(T=sa)};
function Z(){function a(){if(!Y&&(Y=!0,b.calledRun=!0,!K)){V(ea);m(b);if(b.onRuntimeInitialized)b.onRuntimeInitialized();if(b.postRun)for("function"==typeof b.postRun&&(b.postRun=[b.postRun]);b.postRun.length;){var c=b.postRun.shift();fa.unshift(c)}V(fa)}}if(!(0<R)){if(b.preRun)for("function"==typeof b.preRun&&(b.preRun=[b.preRun]);b.preRun.length;)ha();V(da);0<R||(b.setStatus?(b.setStatus("Running..."),setTimeout(function(){setTimeout(function(){b.setStatus("")},1);a()},1)):a())}}b.run=Z;
if(b.preInit)for("function"==typeof b.preInit&&(b.preInit=[b.preInit]);0<b.preInit.length;)b.preInit.pop()();Z();


  return ronchModule.ready
}
);
})();
if (typeof exports === 'object' && typeof module === 'object')
  module.exports = ronchModule;
else if (typeof define === 'function' && define['amd'])
  define([], function() { return ronchModule; });
else if (typeof exports === 'object')
  exports["ronchModule"] = ronchModule;
