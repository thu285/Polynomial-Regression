# include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
String dance338(String x) {
  String sent = x += " keeps the dance floor silly.";
  return sent;
}
