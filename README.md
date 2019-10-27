emcc calculateRonch.cpp -o index.js -s WASM=1 -std=c++11 -s "EXPORTED_FUNCTIONS=['_calcRonch']" -s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall']"
