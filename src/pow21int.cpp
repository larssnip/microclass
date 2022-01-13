#include <Rcpp.h>
#include "pow21inthead.h"
using namespace Rcpp;

int pow21int(int K) {
  int x = 0;
  switch(K) {
  case 0 : x = 1;
    break;
  case 1 : x = 21;
    break;
  case 2 : x = 441;
    break;
  case 3 : x = 9261;
    break;
  case 4 : x = 194481;
    break;
  case 5 : x = 4084101;
    break;
  case 6 : x = 85766121;
    break;
  }
  return x;
}
