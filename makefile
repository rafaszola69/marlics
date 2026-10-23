########### Gnu:
 COMPILER :=g++
 FLAGS    :=-Ofast -fipa-pta -fopenmp -static -march=native
 LIB      :=-lm -lgsl -lgslcblas -lgomp -pthread

############ Intel:
# COMPILER := icpc
# FLAGS    := -static -ipo -fast -no-prec-div -xHost -simd -qopenmp -fp-model fast=2
# LIB      := -mkl -lgsl 


CPPS   := $(wildcard src/*.cpp)
HEADER := $(wildcard src/*.h)
OBJS   := $(patsubst src/%.cpp,build/%.o,${CPPS})

marlics: ${OBJS}
	@${COMPILER} ${FLAGS} ${OBJS} ${LIB} -o marlics
######### Build all cpp's into objects
${OBJS}: build/%.o: src/%.cpp | build
	@${COMPILER} ${FLAGS} -c $< -o $@
######### 	
${OBJS}: ${HEADER}
build:
	@mkdir build
	
######### Compile a debug version	
debug:	$(patsubst %.o,%dbg.o,${OBJS})
	@${COMPILER}  -O0 -g $(filter-out -O% -fast -static, ${FLAGS})  $(patsubst %.o,%dbg.o,${OBJS}) ${LIB} -o marlics_debug
	@rm build/*dbg.o
$(patsubst %.o,%dbg.o,${OBJS}): build/%dbg.o: src/%.cpp | build
	@${COMPILER}  -O0 -g $(filter-out -O% -fast, ${FLAGS}) -c $< -o $@
$(patsubst %.o,%dbg.o,${OBJS}): ${HEADER}

######### Consistency tests (make test). Each test prints its results; the full output goes to build/<test>.log
TESTS    := $(patsubst tests/%.cpp,build/%,$(wildcard tests/*.cpp))
LIB_OBJS := $(filter-out build/main.o,${OBJS})

test: ${TESTS}
	@status=0; for t in ${TESTS}; do ./$$t 2>&1 >$$t.log || status=1; done; \
	if [ $$status -eq 0 ]; then echo "All tests passed."; else echo "SOME TESTS FAILED."; fi; exit $$status
${TESTS}: build/%: tests/%.cpp ${LIB_OBJS} ${HEADER} | build
	@${COMPILER} ${FLAGS} -Isrc $< ${LIB_OBJS} ${LIB} -o $@

######### GPU build (make gpu): marlics_gpu runs on the CPU, or on the GPU with "device gpu" in the input.
# CUDA_HOST: host compiler used by nvcc (old nvcc versions need an old g++, e.g. CUDA 9.1: CUDA_HOST=g++-6).
# CUDA_ARCH: compute_70 PTX is compiled by the driver for newer cards; use e.g. -arch=sm_75 with CUDA >= 10.
NVCC      := nvcc
CUDA_HOST ?= g++
CUDA_ARCH ?= -gencode arch=compute_70,code=compute_70
NVFLAGS   := -O3 -std=c++14 --expt-extended-lambda -Xcudafe --diag_suppress=177 -ccbin ${CUDA_HOST} ${CUDA_ARCH} -Xcompiler -fopenmp -Isrc
GPU_FLAGS := $(filter-out -static,${FLAGS})
GPU_OBJS  := $(patsubst src/%.cpp,build_gpu/%.o,$(filter-out src/gpu_stub.cpp,${CPPS})) build_gpu/gpu_cuda.o

gpu: marlics_gpu
marlics_gpu: ${GPU_OBJS}
	@${COMPILER} ${GPU_FLAGS} ${GPU_OBJS} ${LIB} -lcudart -o marlics_gpu
build_gpu/%.o: src/%.cpp ${HEADER} | build_gpu
	@${COMPILER} ${GPU_FLAGS} -c $< -o $@
build_gpu/gpu_cuda.o: src/gpu/gpu.cu $(wildcard src/gpu/*.h) ${HEADER} | build_gpu
	@${NVCC} ${NVFLAGS} -c $< -o $@
build_gpu:
	@mkdir build_gpu

######### GPU consistency test (make test-gpu): GPU and CPU forces, energies and relaxations agree
test-gpu: build_gpu/test_gpu
	@./build_gpu/test_gpu 2>&1 >build_gpu/test_gpu.log && echo "GPU test passed." || (echo "GPU TEST FAILED."; exit 1)
GPU_LIB_OBJS := $(filter-out build_gpu/main.o build_gpu/gpu_cuda.o,${GPU_OBJS})
build_gpu/test_gpu: build_gpu/test_gpu_cuda.o ${GPU_LIB_OBJS}
	@${COMPILER} ${GPU_FLAGS} build_gpu/test_gpu_cuda.o ${GPU_LIB_OBJS} ${LIB} -lcudart -o $@
build_gpu/test_gpu_cuda.o: tests/gpu/test_gpu.cu src/gpu/gpu.cu $(wildcard src/gpu/*.h) ${HEADER} | build_gpu
	@${NVCC} ${NVFLAGS} -Isrc/gpu -c $< -o $@

########
clean:
	@rm -f marlics marlics_debug marlics_gpu
	@rm -fr build build_gpu
	@rm -f *.csv
	@rm -f *.o
