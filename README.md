# opls

Based on algorithm of Orthogonal Polynomial Curve Fitting from "Orthogonal Polynomial Curve Fitting" by Jeff Reid Orthogonal Polynomial Curve Fitting November 28, 1990

input:
  vector<double> x,y;
  int polynom_degree;
 
 output:
  vector<double> b(polynom_degree+1,0), a(polynom_degree+1,0); ///must be create with correct sizes

run:
  opls Regreess(polynom_degree, x.size());
  Regress.manage_opls(x, y, b, a);
 
 vector a[] gives you regression coefficient of polynom with corresponding index degree
  
