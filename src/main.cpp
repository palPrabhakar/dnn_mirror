#include "dnn.hpp"
#include <chrono>
#include <exception>
#include <iostream>

int main(int argc, char **argv) {
  /*
   * TODO
   * 1. Add provision for restart simulation
   */
  std::chrono::steady_clock::time_point tbegin =
      std::chrono::steady_clock::now();
  std::cout << "\nStarting DNN simulation.\n\n";

  if (argc < 2) {
    std::cerr << "\nInvalid number of arguments, input file missing.\n";
    return 1;
  }

  DNN *dnn;

  if (argc == 3) {
    dnn = new DNN(argv[1], argv[2]);
  } else {
    dnn = new DNN(argv[1]);
  }

  /* DNN dnn(argv[1]); */

  try {
    dnn->read_input_file();
    /* dnn.read_input_file(); */

    if (dnn->get_is_transient()) {
      dnn->init();
      dnn->run_dnn_model();
    } else {
      dnn->init_ss();
      dnn->run_dnn_ss();
    }
  } catch (const std::exception &e) {
    std::cerr << "\nSimulation failed.\n";
    std::cerr << e.what() << "\n";
    delete dnn;
    return 1;
  }

  std::chrono::steady_clock::time_point tend = std::chrono::steady_clock::now();

  std::cout << "\nSimulation finished successfully.\n";
  /* std::cout<<"\nTotal iterations run: "<<dnn->get_citer()<<"\n"; */
  std::cout
      << "\nTotal time taken: "
      << std::chrono::duration_cast<std::chrono::seconds>(tend - tbegin).count()
      << "[s]\n";
  std::cout << "\nTotal time taken: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(tend -
                                                                     tbegin)
                   .count()
            << "[ms]\n";

  delete dnn;

  return 0;
}
