#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int cat1(double value, NumericVector prob) {
  int res=prob.length();
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector z(prob.nrow());
  
  for(int i=0; i<prob.nrow();i++){
    z[i]=cat1(randu[i],prob(i,_));
  }
  return z;
}

//' This function calculates distance between two sets of coordinates
// [[Rcpp::export]]
NumericMatrix GetDistance(NumericMatrix AcCoord,NumericMatrix GridCoord, int Ngrid, int Nac) {
  
  NumericMatrix res(Ngrid,Nac);
  double x2;
  double y2;
  for(int i=0; i<Ngrid;i++){
    for(int j=0; j<Nac; j++){
      x2=pow(GridCoord(i,0)-AcCoord(j,0),2.0);
      y2=pow(GridCoord(i,1)-AcCoord(j,1),2.0);
      res(i,j)=sqrt(x2+y2);
    }
  }
  return res;
}

//' This function gets the maximum across rows
// [[Rcpp::export]]
NumericVector GetMaxRows(NumericMatrix mat) {
  
  NumericVector res(mat.nrow());
  double tmp;
  for(int i=0; i<mat.nrow();i++){
    tmp=R_NegInf;
    for(int j=0; j<mat.ncol(); j++){
      if (mat(i,j)>tmp) tmp=mat(i,j);
    }
    res[i]=tmp;
  }
  return res;
}

//' This function helps calculate mloglikel
// [[Rcpp::export]]
NumericMatrix mloglikel(int ntsegm, int nac, int ngrid, NumericMatrix lprob,IntegerMatrix dat,
                        NumericVector LogTheta){
  NumericMatrix res(ntsegm,nac);
  for(int k=0; k<nac;k++){
    for (int i=0; i<ntsegm; i++){
      for (int j=0; j<ngrid; j++){
        res(i,k)=res(i,k)+dat(i,j)*lprob(k,j);
      }
    }
    res(_,k)=res(_,k)+LogTheta[k];
  }
  return res;
}

//' This function helps calculate theta starting from stick-breaking v's
// [[Rcpp::export]]
NumericVector GetTheta(NumericVector v, int nac){
  NumericVector theta(nac);
  theta[0]=v[0];
  double tmp=1-v[0];
  for(int i=1; i<nac;i++){
    theta[i]=v[i]*tmp;
    tmp=tmp*(1-v[i]);
  }
  return theta;
}