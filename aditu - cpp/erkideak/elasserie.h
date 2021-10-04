#ifndef ELASSERIE_H_INCLUDED
#define ELASSERIE_H_INCLUDED
#include<vector>
#include<complex>

#endif // ELASSERIE_H_INCLUDED

using namespace std;

void elasserie(const vector<double> &zentro, const vector<double> &xpuntu, const double &kp, const double &ks, const double &young, const double &poisson, const int nhmoz, vector<vector<complex<double> > > &serie, vector<vector<complex<double> > > &dserie, vector<vector<complex<double> > > &tserie);
void esfcarserie(const vector<complex<double> > &sbek1, const vector<complex<double> > &sbek2, const vector<vector<double> > &matc, vector<vector<complex<double> > > &serie);
void esfcardserie(const vector<complex<double> > &sbek1, const vector<complex<double> > &sbek2, const vector<complex<double> > &ssbek1, const vector<complex<double> > &ssbek2, const vector<vector<double> > &matc, const vector<vector<vector<double> > >&bimat, const double young, const double poisson, vector<vector<complex<double> > > &dserie, vector<vector<complex<double> > > &tserie);
