#ifndef FourTrans_H_INCLUDED
#define FourTrans_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include<string>
#include<complex>

using namespace std;

class FourTrans
{
public:
    string proiektua;
    unsigned int maizkop;
    unsigned int lerrokop;
    vector<double> espektrua; /** Proiektuko maiztasun-espektrua **/
    vector<vector<complex<double> > > taula; /** Proiektuko maiztasun-Fourier Transformatu taula **/

    /** Eraikitzaileak **/

    FourTrans(string);  /** char: GiDek sortutako .tau testu-fitxategia **/

    /** Desegilea **/

    ~FourTrans();
};

#endif /** FourTrans_H_INCLUDED **/
