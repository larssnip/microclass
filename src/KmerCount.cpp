#include <Rcpp.h>
#include "pow4inthead.h"
#include "pow22inthead.h"
#include "pow64inthead.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix Kmer_count( SEXP seqs, int K, bool names, bool codon ) {
  
  Rcpp::List strings(seqs);
  int num_strings = strings.length();
  IntegerMatrix X(num_strings, pow4int(K));
  int where = 0;
  std::vector<int> Where(K);
  
  for(int i=0; i<K; i++){
    Where[i] = pow4int(K-i-1);  // med K=3 blir Where lik 16 og 4 og 1
  }
  
  for( int i=0; i<num_strings; i++ ) {        // looper over sekvenser
    SEXP s1 = strings[i];                     // s1 er sekvens i
    Rcpp::IntegerVector seq1(s1);             // konverterer til IntegerVector, seq1 er sekvens i
    int num_substr = seq1.length()-K+1;       // antall ord av lengde K i sekvens i
    for( int j=0; j<num_substr; j++ ) {       // looper over alle ord
      if(!codon || j%3==0){                   // sjekk om telling gjøres i leserammen
        where = 0;
        for( int k=0; k<K; k++){                // looper over posisjon i ord
          where += seq1[j+k]*Where[k];          // where blir kolonnen til ordet i X, beregnet i 4-talls systemet
        }                                       // dette er alltid et tall fra 0 til (4^K)-1, med mindre
        if(where >= 0){                         // ett av elementene i sekvensen har verdien -4^K, da blir
          ++X(i, where);                        // where negativ, og ordet ignoreres 
        }
      }
    }
  }

  if(names){
    int N = pow4int(K);
    Rcpp::CharacterVector ACGT = Rcpp::CharacterVector::create("A","C","G","T");
    Rcpp::CharacterVector ACGTs(N);
    std::vector< std::vector< std::string > > matr;
    matr.resize( K , std::vector<std::string>( N ) );
    Rcpp::CharacterVector cnms(N);
    for(int i=0; i<K; i++){
      ACGTs = rep(rep_each(ACGT, pow4int(K-i-1)), pow4int(i));
      matr[i] = Rcpp::as< std::vector< std::string > >(ACGTs);
    }
    for(int i=0; i<N; i++){
      std::stringstream ss;
      for(int j=0; j<K; j++){
        ss << matr[j][i];
      }
      cnms[i] = ss.str();      
    }
    Rcpp::List dimnms = Rcpp::List::create(R_NilValue, cnms);
    X.attr("dimnames") = dimnms;
  }
  
  return X;
}

// [[Rcpp::export]]
IntegerMatrix Kmer_count_amino_acid( SEXP seqs, int K, bool names ) {
  
  Rcpp::List strings(seqs);
  int num_strings = strings.length();
  IntegerMatrix X(num_strings, pow22int(K));
  int where = 0;
  std::vector<int> Where(K);
  
  for(int i=0; i<K; i++){
    Where[i] = pow22int(K-i-1);  // med K=3 blir Where lik 484 og 22 og 1
  }
  
  for( int i=0; i<num_strings; i++ ) {        // looper over sekvenser
    SEXP s1 = strings[i];                     // s1 er sekvens i
    Rcpp::IntegerVector seq1(s1);             // konverterer til IntegerVector, seq1 er sekvens i
    int num_substr = seq1.length()-K+1;       // antall ord av lengde K i sekvens i
    for( int j=0; j<num_substr; j++ ) {       // looper over alle ord
      where = 0;
      for( int k=0; k<K; k++){                // looper over posisjon i ord
        where += seq1[j+k]*Where[k];          // where blir kolonnen til ordet i X, beregnet i 4-talls systemet
      }                                       // dette er alltid et tall fra 0 til (4^K)-1, med mindre
      if(where >= 0){                         // ett av elementene i sekvensen har verdien -4^K, da blir
        ++X(i, where);                        // where negativ, og ordet ignoreres 
      }
    }
  }
  
  if(names){
    int N = pow22int(K);
    Rcpp::CharacterVector ACGT = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z"};
    Rcpp::CharacterVector ACGTs(N);
    std::vector< std::vector< std::string > > matr;
    matr.resize( K , std::vector<std::string>( N ) );
    Rcpp::CharacterVector cnms(N);
    for(int i=0; i<K; i++){
      ACGTs = rep(rep_each(ACGT, pow22int(K-i-1)), pow22int(i));
      matr[i] = Rcpp::as< std::vector< std::string > >(ACGTs);
    }
    for(int i=0; i<N; i++){
      std::stringstream ss;
      for(int j=0; j<K; j++){
        ss << matr[j][i];
      }
      cnms[i] = ss.str();      
    }
    Rcpp::List dimnms = Rcpp::List::create(R_NilValue, cnms);
    X.attr("dimnames") = dimnms;
  }
  
  return X;
}

// [[Rcpp::export]]
IntegerMatrix Kmer_count_codon( SEXP seqs, int K, bool names ) {
  
  Rcpp::List strings(seqs);
  int num_strings = strings.length();
  IntegerMatrix X(num_strings, pow64int(K));
  int where = 0;
  std::vector<int> Where(K);
  
  for(int i=0; i<K; i++){
    Where[i] = pow64int(K-i-1);  // med K=3 blir Where lik 4096 og 64 og 1
  }
  
  for( int i=0; i<num_strings; i++ ) {        // looper over sekvenser
    SEXP s1 = strings[i];                     // s1 er sekvens i
    Rcpp::IntegerVector seq1(s1);             // konverterer til IntegerVector, seq1 er sekvens i
    int num_substr = seq1.length()-K+1;       // antall ord av lengde K i sekvens i
    for( int j=0; j<num_substr; j++ ) {       // looper over alle ord
      where = 0;
      for( int k=0; k<K; k++){                // looper over posisjon i ord
        where += seq1[j+k]*Where[k];          // where blir kolonnen til ordet i X, beregnet i 4-talls systemet
      }                                       // dette er alltid et tall fra 0 til (4^K)-1, med mindre
      if(where >= 0){                         // ett av elementene i sekvensen har verdien -4^K, da blir
        ++X(i, where);                        // where negativ, og ordet ignoreres 
      }
    }
  }
  
  if(names){
    int N = pow64int(K);
    Rcpp::CharacterVector ACGT = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","0","1","2","3","4","5","6","7","8","9","10","-","+"};
    Rcpp::CharacterVector ACGTs(N);
    std::vector< std::vector< std::string > > matr;
    matr.resize( K , std::vector<std::string>( N ) );
    Rcpp::CharacterVector cnms(N);
    for(int i=0; i<K; i++){
      ACGTs = rep(rep_each(ACGT, pow64int(K-i-1)), pow64int(i));
      matr[i] = Rcpp::as< std::vector< std::string > >(ACGTs);
    }
    for(int i=0; i<N; i++){
      std::stringstream ss;
      for(int j=0; j<K; j++){
        ss << matr[j][i];
      }
      cnms[i] = ss.str();      
    }
    Rcpp::List dimnms = Rcpp::List::create(R_NilValue, cnms);
    X.attr("dimnames") = dimnms;
  }
  
  return X;
}
