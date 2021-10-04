#ifndef OSKOSERIE_H_INCLUDED
#define OSKOSERIE_H_INCLUDED
#include<vector>
#include<complex>

#endif // OSKOSERIE_H_INCLUDED

using namespace std;

void oskoserie(const vector<double> &zentro, const vector<double> &xpuntu, const vector<double> &normal, const double zo, const double kp, const double ks, const double young, const double poisson, const int nhmoz, vector<vector<complex<double> > > &serie, vector<vector<complex<double> > > &dserie, vector<vector<complex<double> > > &tserie);
