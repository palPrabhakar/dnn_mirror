# Dendritic Needle Network (DNN) by Tourret et. al.

## Current State of Branch
1. Reproduced the thin/thick needle Ivantsov graph
2. CUDA Implementation of 2d thick isolated needle growth 

### Build DNN
Use included Makefile

### Run all test case 
Use included Makefile in the test folder

#### This is not updated
1. cd to test dir
2. g++ -std=c++11 -Icatch.hpp -Wall test_dnn.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm & ./main_test --success

### Run specific test case
1. cd to test dir
2. g++ -std=c++11 -Icatch.hpp -Wall test_dnn.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm & ./main_test "<test_name>"

## Run SS solver test case
1. cd to test dir
2. g++ -std=c++11 -Icatch.hpp -Wall -DNEEDLE=0 test_ss_solver.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm && ./main_test --success

## Run Transient solver test case
1. cd to test dir
2. g++ -std=c++11 -Icatch.hpp -Wall -DNEEDLE=0 -DGROW=0 -O4 test_transient_solver.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm && ./main_test --success

## Run J-Integral test
1. cd to test dir
2. g++ -std=c++11 -Icatch.hpp -Wall -DNEEDLE=1 -DGROW=0 -DNON_DIM=0 -O4 test_j_integral.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm && ./main_test --success

## run thin needle test
1. g++ -std=c++11 -Icatch.hpp -Wall -O4 test_sharp_needle.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm && ./main_test --success

## run thick needle test
1. g++ -std=c++11 -Icatch.hpp -Wall -O4 -DNON_DIM=0 test_thick_needle.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm && ./main_test --success

## run 2d parabolic needle test
1. g++ -std=c++11 -Icatch.hpp -Wall -O4 test_2d_parabolic_needle.cpp ../src/dnn.cpp ../src/grid.cpp main.cpp -o main_test -lgsl -lgslcblas -lm && ./main_test --success


## Restart Simulations 
To restart the simulation, the following files are required:
1. Original input file 
2. Dump output (how will you know what the simulation time or the no. of iterations done at the time of the restart? Good Question: Let's make a list of information required to restart the simulation. For the time being let's write a module to read the dump file. //one of the idea is to use the restart file name as <dump_current_iter> and once the current iter is known then the simulation time can be inferred.) 


## TODO (Probably this is fixed)
1. Fix the dx = dy = dx in FIF calculation (not doing this right now, it won't matter as long the uniform grid spacing is used in the x and y direction)

### Make key and their values case insensitive in input parameters

