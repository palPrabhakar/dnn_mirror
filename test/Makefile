SRC_DIR = ../src/

CXX = g++

CFLAGS = -std=c++11 -Wall -pedantic -DGPU=0

INCLUDES = -I../src/include/
NVCC_INCLUDES = -I/opt/cuda-10.1/samples/common/inc/ -I/opt/cuda-10.1/include/

LIBS = -lgsl -lgslcblas -lm

all: test_2d_parabolic_needle.o test_ss_solver.o test_j_integral.o test_transient_solver.o test_sharp_needle.o dnn.o grid.o helper.o main.o 
	$(CXX) -o main_test test_2d_parabolic_needle.o test_ss_solver.o test_j_integral.o test_transient_solver.o test_sharp_needle.o dnn.o grid.o main.o helper.o $(LIBS)

main.o: main.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c main.cpp

2d_parabolic: test_2d_parabolic_needle.o dnn.o grid.o helper.o main.o 
	$(CXX) -o main_test test_2d_parabolic_needle.o dnn.o grid.o main.o helper.o $(LIBS)

# thick_needle: test_thick_needle.o dnn.o grid.o helper.o main.o 
# 	$(CXX) -o main_test test_thick_needle.o dnn.o grid.o main.o helper.o $(LIBS)

thin_needle: test_sharp_needle.o dnn.o grid.o helper.o main.o 
	$(CXX) -o main_test test_sharp_needle.o dnn.o grid.o main.o helper.o $(LIBS)

ss_solver: test_ss_solver.o dnn.o grid.o helper.o main.o 
	$(CXX) -o main_test test_ss_solver.o dnn.o grid.o main.o helper.o $(LIBS)

j_integral: test_j_integral.o dnn.o grid.o helper.o main.o 
	$(CXX) -o main_test test_j_integral.o dnn.o grid.o main.o helper.o $(LIBS)

transient_solver: test_transient_solver.o dnn.o grid.o helper.o main.o 
	$(CXX) -o main_test test_transient_solver.o dnn.o grid.o main.o helper.o $(LIBS)

test_j_integral.o: test_j_integral.cpp $(SRC_DIR)dnn.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c test_j_integral.cpp $(SRC_DIR)dnn.cpp

test_ss_solver.o: test_ss_solver.cpp $(SRC_DIR)dnn.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c test_ss_solver.cpp $(SRC_DIR)dnn.cpp

test_transient_solver.o: test_transient_solver.cpp $(SRC_DIR)dnn.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c test_transient_solver.cpp $(SRC_DIR)dnn.cpp

test_thick_needle.o: test_thick_needle.cpp $(SRC_DIR)dnn.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c test_thick_needle.cpp $(SRC_DIR)dnn.cpp

test_sharp_needle.o: test_sharp_needle.cpp $(SRC_DIR)dnn.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c test_sharp_needle.cpp $(SRC_DIR)dnn.cpp

test_2d_parabolic_needle.o: test_2d_parabolic_needle.cpp $(SRC_DIR)dnn.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c test_2d_parabolic_needle.cpp $(SRC_DIR)dnn.cpp

dnn.o: $(SRC_DIR)dnn.cpp $(SRC_DIR)dnn.hpp $(SRC_DIR)data_types.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRC_DIR)dnn.cpp

helper.o: $(SRC_DIR)helper.cpp $(SRC_DIR)helper.hpp
	$(CXX) $(CFLAGS) $(NVCC_INCLUDES) $(INCLUDES) -c $(SRC_DIR)helper.cpp

grid.o: $(SRC_DIR)grid.cpp $(SRC_DIR)grid.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $(SRC_DIR)grid.cpp

clean: 
	rm -rf *.o main_test
