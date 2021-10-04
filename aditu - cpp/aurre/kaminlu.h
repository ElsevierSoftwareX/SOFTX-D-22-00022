#ifndef KAMINLU_H_INCLUDED
#define KAMINLU_H_INCLUDED
#include<vector>
#include<complex>

#endif // KAMINLU_H_INCLUDED

using namespace std;

void kaminlu(const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, vector<complex<double> > &ebazpena, float &errorea);
void mirabiderm (const vector<vector<complex<double> > > &amat, vector<vector<complex<double> > > &aia);
void mirabiderb (const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<complex<double> > &aib);
void hondarkalk (const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, const vector<complex<double> > &ebazpena, float &errorea);
double bnorm (const vector<complex<double> > &b);
