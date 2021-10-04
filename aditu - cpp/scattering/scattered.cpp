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

void scattered(vector<double> &xnodo, double k, int nhmoz, complex<double> &press)
{
double r0=1.0;
double r, kos;
complex<double> am;
vector<complex<double> > hankelfun1; /** r=r0 denerako **/
vector<complex<double> > dhankelfun1;
vector<complex<double> > hankelfunr;
vector<complex<double> > dhankelfunr;

double lp;
int m, kont;

press=0.0+0.0*I;

r=0;

for (kont=0; kont<3; kont++)
{
    r=r+(xnodo[kont]*xnodo[kont]);
}
r=sqrt(r);

kos=xnodo[2]/r;

/** 0 ordenako Hankel funtzioak **/
hankelfunr.push_back(-I*exp(I*k*r)/(k*r));
hankelfun1.push_back(-I*exp(I*k*r0)/(k*r0));

for (m=0; m<=nhmoz; m++)
{
    /** m+1 ordenako Hankel funtzioak **/

    hankelfunr.push_back(gsl_sf_bessel_jl(m+1+0.5, k*r) + I*gsl_sf_bessel_yl(m+1+0.5, k*r));
    hankelfun1.push_back(gsl_sf_bessel_jl(m+1+0.5, k*r0)+ I*gsl_sf_bessel_yl(m+1+0.5, k*r0));

    /** m ordenako Hankel funtzioen deribatuak **/
    dhankelfunr.push_back((m/r)*hankelfunr[m]-k*hankelfunr[m+1]);
    dhankelfun1.push_back((m/r)*hankelfun1[m]-k*hankelfun1[m+1]);

    /** Pm Legendre Polinomioak **/
    lp=gsl_sf_legendre_Pl(m, kos);
    am=-pow(I,m)*(2*m+1)*real(dhankelfun1[m]);
    am=am/dhankelfun1[m];
    press=press-hankelfunr[m]*lp*am;
}

}
