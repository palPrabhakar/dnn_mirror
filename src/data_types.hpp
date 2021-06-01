#ifndef DATA_TYPES
#define DATA_TYPES

#include <string>

struct BC {
  std::string type;
  double value;
  void set(std::string _type, double val) {
    type = _type;
    value = val;
  }
};

struct Needle {
  size_t x0, y0; // needle origin co-ordinates
  size_t xf, yf; // needle tip co-ordinates
  double r;
  double vel;
  double rad;
};

struct Point {
  Point(double c0, bool _needle) {
    u = c0;
    needle = _needle;
  }

  double u;
  bool needle;
};

class exit_exception : public std::exception {
  std::string msg;

public:
  exit_exception() {}
  exit_exception(const std::string &msg) : msg(msg) {}

  virtual const char *what() const throw() { return msg.c_str(); }
};

class type_exception : public std::exception {
  const char *msg;

public:
  type_exception(const char *msg) : msg(msg) {}

  virtual const char *what() const throw() { return msg; }
};

#endif /* ifndef DATA_TYPES                                                    \
                                                                               \
 */
