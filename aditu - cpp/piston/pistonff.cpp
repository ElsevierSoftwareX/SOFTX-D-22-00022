#include <iostream>
#include <vector>
#include <complex>
#include <boost/math/special_functions/bessel.hpp>
using namespace std;

const complex<double> I{0.0, 1.0};

void pistonff( vector<double> &xnodo, double k, double a, complex<double> &baliozehatz )
{
   double r = sqrt(xnodo[0]*xnodo[0]+xnodo[1]*xnodo[1]+xnodo[2]* xnodo[2]);
   double theta = acos(xnodo[0]/r);

    if (theta>0.001)
    {
    double besselfun1 = boost::math::cyl_bessel_j(1,k*a*sin(theta));     // 1st order Bessel function
    baliozehatz = -1.21*pow(a,2)*exp(I*k*r)*2.0*besselfun1/(k*a*sin(theta))/(2*r);
    }
    else
    {
    baliozehatz = -1.21*pow(a,2)*exp(I*k*r)/(2*r);
    }
}
