#include <complex>
#include <complex.h>
#include <vector>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_sf_bessel.h>
#ifndef PI
    #define PI 3.14159265358979323846
#endif

using namespace std;

void gerdes(vector<double> &zentro, vector<double> &xnodo, double k, int nhmoz, complex<double> &hbal)
{

double koordc[3];
double r, theta;
double cbal;
vector<double> besselfun1; /** r=1 denerako **/
vector<double> dbesselfun1;

double lp;
int m, kont;

hbal=0.0+0.0*I;

r=0;

for (kont=0; kont<3; kont++)
{
    r=r+(xnodo[kont]-zentro[kont])*(xnodo[kont]-zentro[kont]);
    koordc[kont]=(xnodo[kont]-zentro[kont]);
}
r=sqrt(r);

theta=acos(koordc[2]/r);

/** 0 ordenako Bessel funtzioak **/
cbal=gsl_sf_bessel_jl(0.5, k);
besselfun1.push_back(cbal);

for (m=0; m<=nhmoz; m++)
{
    /** m+1 ordenako Bessel funtzioak **/
    cbal=gsl_sf_bessel_jl(m+1+0.5, k);
    besselfun1.push_back(cbal);

    /** m ordenako Bessel funtzioen deribatuak **/
    cbal=besselfun1[m]*m-k*besselfun1[m+1];
    dbesselfun1.push_back(cbal);

    /** Pm Legendre Polinomioak **/

    lp=gsl_sf_legendre_Pl(m, cos(theta));
    hbal=hbal+((2*m+1)*pow(I,m)*lp*dbesselfun1[m]);
}

}
