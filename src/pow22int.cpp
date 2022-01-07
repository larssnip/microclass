#include <Rcpp.h>
#include "pow22inthead.h"
using namespace Rcpp;

int pow22int(int K) {
  int x = 0;
  switch(K) {
  case 0 : x = 1;
    break;
  case 1 : x = 22;
    break;
  case 2 : x = 484;
    break;
  case 3 : x = 10648;
    break;
  case 4 : x = 234256;
    break;
  case 5 : x = 5153632;
    break;
  case 6 : x = 113379904;
    break;
  }
  return x;
}
