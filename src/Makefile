CXX = g++
NVCC = nvcc -ccbin g++

CFLAGS = -std=c++11 -Wall -pedantic
NVCCFLAGS = -m64 -gencode arch=compute_61,code=sm_61

INCLUDES := -I./include/ -I/opt/cuda-10.1/include/
NVCC_INCLUDES = -I/opt/cuda-10.1/samples/common/inc/

LIBS = -lgsl -lgslcblas -lm  -lcudart

MAIN = dnn

all: $(MAIN)
	@echo DNN compiled.

$(MAIN): grid.o dnn.o main.o dnn_gpu.o helper.o
	$(NVCC) -O2 -o $(MAIN) main.o grid.o dnn.o dnn_gpu.o helper.o $(LIBS)

# $(MAIN): grid.o dnn.o main.o 
# 	$(CXX) $(CFLAGS) -O -o $(MAIN) main.o grid.o dnn.o $(LIBS)

grid.o: grid.cpp grid.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c grid.cpp

dnn.o: dnn.cpp dnn.hpp data_types.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c dnn.cpp

main.o: main.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c main.cpp

helper.o: helper.cpp helper.hpp
	$(CXX) $(CFLAGS) $(NVCC_INCLUDES) $(INCLUDES) -c helper.cpp

dnn_gpu.o: dnn_gpu.cu data_types.hpp
	$(NVCC) $(NVCCFLAGS) $(NVCC_INCLUDES) -c dnn_gpu.cu

clean:
	$(RM) $(MAIN) *.o
