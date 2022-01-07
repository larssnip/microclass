#include <Rcpp.h>
#include "pow64inthead.h"
using namespace Rcpp;

int pow64int(int K) {
  int x = 0;
  switch(K) {
  case 0 : x = 1;
    break;
  case 1 : x = 64;
    break;
  case 2 : x = 4096;
    break;
  case 3 : x = 262144;
    break;
  case 4 : x = 16777216;
    break;
  case 5 : x = 1073741824;
    break;
  }
  return x;
}
