#include <iostream>
#include <vector>
#include <complex>
#include "eragiketak.h"
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_sf_bessel.h>
using namespace std;

const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

template <typename T, typename U>
inline complex<T> operator*(complex<T> lhs, const U& rhs)
{
    return lhs *= rhs;
}

void neumann(const vector<double> &xnodo, const double k, const int nhmoz, complex<double> &hbal)
{
double r, kos;
complex<double> am;
vector<double> besselfun1; /** r=1 denerako **/
vector<double> dbesselfun1;

double lp;
int m, kont;

hbal=0.0+0.0*I;

r=0;

for (kont=0; kont<3; kont++)
{
    r=r+xnodo[kont]*xnodo[kont];
}
r=sqrt(r);

kos=xnodo[2]/r;

/** 0 ordenako Bessel funtzioa **/
besselfun1.push_back(sin(k*r)/(k*r));

for (m=0; m<=nhmoz; m++)
{
    /** m+1 ordenako Bessel funtzioa **/
    besselfun1.push_back(gsl_sf_bessel_jl(m+1+0.5, k*r));

    /** m ordenako Bessel funtzioaren deribatua **/
    dbesselfun1.push_back((m/r)*besselfun1[m]-k*besselfun1[m+1]);

    /** Pm Legendre Polinomioak **/
    lp=gsl_sf_legendre_Pl(m, kos);
    am=-pow(I,m)*(2*m+1)*lp*dbesselfun1[m];
    hbal=hbal+am;
}

}
