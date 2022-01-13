#include <Rcpp.h>
#include "pow2inthead.h"
using namespace Rcpp;

int pow2int(int K) {
  int x = 0;
  switch(K) {
  case 0 : x = 1;
    break;
  case 1 : x = 2;
    break;
  case 2 : x = 4;
    break;
  case 3 : x = 8;
    break;
  case 4 : x = 16;
    break;
  case 5 : x = 32;
    break;
  case 6 : x = 64;
    break;
  case 7 : x = 128;
    break;
  case 8 : x = 256;
    break;
  case 9 : x = 512;
    break;
  case 10 : x = 1024;
    break;
  case 11 : x = 2048;
    break;
  case 12 : x = 4096;
    break;
  case 13 : x = 8192;
    break;
  case 14 : x = 16384;
    break;
  case 15 : x = 32768;
    break;
  case 16 : x = 65536;
    break;
  case 17 : x =  131072;
    break;
  case 18 : x = 262144;
    break;
  case 19 : x = 524288;
    break;
  case 20 : x = 1048576;
    break;
  case 21 : x = 2097152;
    break;
  case 22 : x = 4194304;
    break;
  case 23 : x = 8388608;
    break;
  case 24 : x = 16777216;
    break;
  case 25 : x = 33554432;
    break;
  case 26 : x = 67108864;
    break;
  case 27 : x = 134217728;
    break;
  case 28 : x = 268435456;
    break;
  case 29 : x = 536870912;
    break;
  case 30 : x = 1073741824;
    break;
  }
  return x;
}
