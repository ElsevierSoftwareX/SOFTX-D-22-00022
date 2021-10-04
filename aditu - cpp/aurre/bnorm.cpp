#include "bmorm.h"
#include <complex>
const complex<double> I{0.0, 1.0};

using namespace std;

void bnorm (vector<complex<double> > &b, &double norm)
/**
** Input: b
** Output: dot =||b||
**/
{
    double dot=0;
    complex<double> dk=0.0+0.0*I;
    unsigned int i;
    unsigned int k = b.size();

    for (i=0; i < k; i++)
    {
        dk=dk+b[i]*conj(b[i]);
    }

    norm=sqrt(abs(dk));
}
