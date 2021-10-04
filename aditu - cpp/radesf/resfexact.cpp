#include <iostream>
#include <vector>
#include <complex>           // <=====
using namespace std;

const double PI = 3.14159265358979323846;
const complex<double> I{ 0.0, 1.0 };           // <=====

// azelerazioa=1/k^2;

void resfexact( vector<double> xpuntu, double k, complex<double> &pbal )
{
   double r0 = 1.0;
   double r = sqrt(xpuntu[0]*xpuntu[0] + xpuntu[1]*xpuntu[1] + xpuntu[2]*xpuntu[2]);

   pbal = pow(r0,2)*(exp(-I*k*r0)/(r*(1.0-I*k*r0)))*exp(I*k*r); /** rho*azelerazioa = 1.0 **/

   //pbal = -1.21*pow(r0,2)*exp(I*k*(r-r0))/(r*(1.0-I*k*r0));
   //pbal = -1.21*pow(r0,2)*exp(I*k*(r-r0))/(r*(1.0-I*k*r0)*pow(k,2));
   }
