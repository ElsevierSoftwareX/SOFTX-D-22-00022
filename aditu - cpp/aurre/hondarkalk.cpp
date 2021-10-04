#include "hondarkalk.h"
#include <complex>
const complex<double> I{0.0, 1.0};

using namespace std;

void hondarkalk (vector<vector<complex<double> > > &amat, vector<complex<double> > &bbek, vector<complex<double> > &errbek, , vector<complex<double> > &ebazpena, vector<complex<double> > &hondarra)
{
/**
** Inputs: amat, bbek, ebazpena
** Output: errbek, hondarra
** bh= [amat(1:m,1:n)] * [ebazpena(1:n)]
**/
    unsigned int m=amat.size();
    unsigned int n=amat[0].size();
    unsigned int i, j;
    complex<double> dk;

    hondarra.clear();
    hondarra.resize(m, 0.0+0.0*I);

    errbek.clear();
    errbek.resize(m, 0.0+0.0*I);

    for (i=0; i < m; i++)
    {
        dk=0.0+0.0*I;
        for (j=0; j < n; j++)
        {
            dk=dk+amat[i][j]*ebazpena[j];
            errbek[i]=bbek[i];
        }
        if (abs(bbek[i]) < 1.e-6)
        {
            j=0;
            while (abs(amat[i][j]*ebazpena[j]) < 1.e-6 && j < m)
            {
                j++;
            }
            dk=dk-amat[i][j]*ebazpena[j];
            errbek[i]=bbek[i]-amat[i][j]*ebazpena[j];
        }
        else
        {
            errbek[i]=bbek[i];
        }
        hondarra.push_back(errbek[i]-dk);
    }
}
