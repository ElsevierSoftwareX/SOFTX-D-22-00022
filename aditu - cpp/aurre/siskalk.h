#ifndef SISKALK_H_INCLUDED
#define SISKALK_H_INCLUDED
#include<vector>
#include<complex>
#include "../erkideak/FourTrans.h"
#include "../erkideak/elasserie.h"
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
#include "../erkideak/ObjektuPuntu.h"

#endif // SISKALK_H_INCLUDED

using namespace std;

void siskalk(const vector<ObjektuPuntu> &puntu, const unsigned int &objkop, const vector<vector<unsigned int> > &ekkop, const int harmonikoa, const unsigned int maizkont, const string taula, vector<vector<complex<double> > > &amat, vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, unsigned int &prokoefkop);
