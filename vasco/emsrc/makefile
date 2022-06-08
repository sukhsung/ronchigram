all : calculateRonch.cpp kiss_fftnd.c kiss_fft.c
	emcc -O3 --closure 1 calculateRonch.cpp kiss_fftnd.c kiss_fft.c -o ../index.js \
		-s WASM=1 \
		-s EXPORTED_FUNCTIONS='["_calcRonch","_free", "_malloc"]' \
		-s DEMANGLE_SUPPORT=1 \
		-s MODULARIZE=1 \
		-s EXPORT_NAME="ronchModule" \
		-s ALLOW_MEMORY_GROWTH=1 \
		-s EXPORTED_RUNTIME_METHODS='["ccall"]'

test : calculateRonch.cpp kiss_fftnd.c kiss_fft.c
	emcc -O3 --closure 1 calculateRonch.cpp kiss_fftnd.c kiss_fft.c -o ../index.js \
		-s WASM=1 \
		-s EXPORTED_FUNCTIONS='["_calcRonch","_free", "_malloc"]' \
		-s DEMANGLE_SUPPORT=1 \
		-s MODULARIZE=1 \
		-s EXPORT_NAME="ronchModule" \
		-s EXPORTED_RUNTIME_METHODS='["ccall"]' \
		-s ASSERTIONS=1 \
		-s ALLOW_MEMORY_GROWTH=1
		
debug : calculateRonch.cpp kiss_fftnd.c kiss_fft.c
	emcc -O3 calculateRonch.cpp kiss_fftnd.c kiss_fft.c -o ../index.js \
		-s WASM=1 \
		-s EXPORTED_FUNCTIONS='["_calcRonch","_free", "_malloc"]' \
		-s DEMANGLE_SUPPORT=1 \
		-s MODULARIZE=1 \
		-s EXPORT_NAME="ronchModule" \
		-s EXPORTED_RUNTIME_METHODS='["ccall"]' \
		-s SAFE_HEAP=1 \
		-s ASSERTIONS=1 \
		-s NO_FILESYSTEM=1 \