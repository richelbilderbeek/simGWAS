#include<Rcpp.h>
using namespace Rcpp;

/* convert 0|1 to haps */

// [[Rcpp::export]]
IntegerMatrix vcf2haps( const CharacterMatrix& x) {
  int n = x.nrow();
  int m = x.ncol();
  IntegerMatrix y(n,2*m);
  for(int i=0; i<n; i++) {
    for(int j=0; j<m; j++) {
      // Rcout << x(i,j) << ' ' << x(i,j)[0] << "/" << x(i,j)[2] <<  "\n";
      y(i,j) = x(i,j)[0] - '0';
      y(i,j+m) = x(i,j)[2] - '0';
    }
  }
  return(y);
}
