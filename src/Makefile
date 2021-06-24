# Set this value to 1 to build on GPU 
GPU = 0

CXX = g++

MAIN = dnn

all: $(MAIN)
	@echo DNN compiled.

ifeq ($(GPU), 1)
NVCC = nvcc -ccbin g++
INCLUDES := -I./include/ -I/opt/cuda-10.1/include/
LIBS = -lgsl -lgslcblas -lm  -lcudart
NVCCFLAGS = -m64 -gencode arch=compute_61,code=sm_61
NVCC_INCLUDES = -I/opt/cuda-10.1/samples/common/inc/
CFLAGS = -std=c++11 -Wall -pedantic

$(MAIN): grid.o dnn.o main.o dnn_gpu.o helper.o
	$(NVCC) -O2 -o $(MAIN) main.o grid.o dnn.o dnn_gpu.o helper.o $(LIBS)

else
INCLUDES := -I./include/
LIBS = -lgsl -lgslcblas -lm
CFLAGS = -std=c++11 -Wall -pedantic -DGPU=0

$(MAIN): grid.o dnn.o main.o helper.o
	$(CXX) -O2 -o $(MAIN) main.o grid.o dnn.o helper.o $(LIBS)

endif

# $(MAIN): grid.o dnn.o main.o 
# 	$(CXX) $(CFLAGS) -O -o $(MAIN) main.o grid.o dnn.o $(LIBS)

grid.o: grid.cpp grid.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c grid.cpp

dnn.o: dnn.cpp dnn.hpp data_types.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c dnn.cpp

main.o: main.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c main.cpp

ifeq ($(GPU), 1)
helper.o: helper.cpp helper.hpp
	$(CXX) $(CFLAGS) $(NVCC_INCLUDES) $(INCLUDES) -c helper.cpp

dnn_gpu.o: dnn_gpu.cu data_types.hpp
	$(NVCC) $(NVCCFLAGS) $(NVCC_INCLUDES) -c dnn_gpu.cu

else
helper.o: helper.cpp helper.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c helper.cpp

endif

clean:
	$(RM) $(MAIN) *.o
