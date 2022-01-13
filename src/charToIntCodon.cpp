#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List charToIntCodon(CharacterVector Seq) {
  std::vector<std::string> seq = as<std::vector<std::string> >(Seq);
  std::vector<std::vector<int> > seq_int;
  std::vector<int> row;
  int acgt;
  
  
  // Convert strings to integers
  for(unsigned i=0; i<seq.size(); ++i){
    std::vector<int> row;
    for(std::string::iterator it = seq[i].begin(), end = seq[i].end(); it != end; ++it) {
      acgt = *it;
      if(acgt > 64 && acgt < 91)       // A-Z
        row.push_back(acgt-65);
      else if(acgt > 96 && acgt < 123) // a-z
        row.push_back(acgt-71);
      else if(acgt > 47 && acgt < 58)  // 0-9
        row.push_back(acgt+4);
      else if(acgt == 45)              // -
        row.push_back(62);
      else if(acgt == 43)              // +
        row.push_back(63);
      else
        row.push_back(-1073741824); // -64^5
    }
    seq_int.push_back(row);
  }
  return wrap(seq_int);
}
