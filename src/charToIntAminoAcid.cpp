#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List charToIntAminoAcid(CharacterVector Seq) {
  std::vector<std::string> seq = as<std::vector<std::string> >(Seq);
  std::vector<std::vector<int> > seq_int;
  std::vector<int> row;
  int acgt;


  // Convert strings to integers
  for(unsigned i=0; i<seq.size(); ++i){
    std::vector<int> row;
    for(std::string::iterator it = seq[i].begin(), end = seq[i].end(); it != end; ++it) {
      acgt = *it;
      switch (acgt) {
      case 65:
        row.push_back(0); // A
        break;
      case 67:
        row.push_back(1); // C
        break;
      case 68:
        row.push_back(2); // D
        break;
      case 69:
        row.push_back(3); // E
        break;
      case 70:
        row.push_back(4); // F
        break;
      case 71:
        row.push_back(5); // G
        break;
      case 72:
        row.push_back(6); // H
        break;
      case 73:
        row.push_back(7); // I
        break;
      case 75:
        row.push_back(8); // K
        break;
      case 76:
        row.push_back(9); // L
        break;
      case 77:
        row.push_back(10); // M
        break;
      case 78:
        row.push_back(11); // N
        break;
      case 80:
        row.push_back(12); // P
        break;
      case 81:
        row.push_back(13); // Q
        break;
      case 82:
        row.push_back(14); // R
        break;
      case 83:
        row.push_back(15); // S
        break;
      case 84:
        row.push_back(16); // T
        break;
      case 86:
        row.push_back(17); // V
        break;
      case 87:
        row.push_back(18); // W
        break;
      case 89:
        row.push_back(19); // Y
        break;
      case 42:
        row.push_back(20); // *
        break;
      default:
        row.push_back(-113379904); // -22^6
        break;
      }
    }
    seq_int.push_back(row);
  }
  return wrap(seq_int);
}
