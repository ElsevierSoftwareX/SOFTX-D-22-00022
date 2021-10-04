#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include "kaminlu.h"

const complex<double> I{0.0, 1.0};

using namespace std;

void kaminlu(const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, vector<complex<double> > &ebazpena, float &errorea)
{
/**
** Inputs: amat, bbek
** Output: ebazpena
** amat*ebazpena=bbek sistemaren Karratu Minimoren araberako ebazpen
** hurbildua L deskonposizioaren bidez
** L*U*ebazpena=amat =>
** U*y
    **/

vector<complex<double> > ck;

/**
** Karratu minimoetako matrize eta bektorea
**/

vector<vector<complex<double> > > aia;
vector<complex<double> > aib;
vector<complex<double> > y;
complex<double> lag;
double handiena;

/** LU faktorizaziokoak **/

vector<vector<complex<double> > > L;
vector<vector<complex<double> > > U;

unsigned int n, i, j, k, s;

n=amat[0].size();

vector<complex<double> > errbek;
aia.clear();
aib.clear();
lag=0.0+0.0*I;
ebazpena.clear();
ebazpena.resize(n, 0.0+0.0*I);

/**
** [amat(1:m,1:n)]T * [amat(1:m,1:n] = aia[n,n]
**/
mirabiderm (amat, aia);

/**
** [amat(1:m,1:n)]T * [b(1:m,1)] = aib[n]
**/
mirabiderb (amat, bbek, aib);

/**
** LU with partial decomposition
** Lehenik pibotaje partziala aia matrizean
**/

for(j=0; j<n; j++)
{
    handiena=0;
    k=j;
    for (i=j; i<n; i++)
    {
        if (abs(aia[i][j]) > handiena)
        {
            k=i;
            handiena=abs(aia[i][j]);
        }
    }
    if (j!=k)
    { /** aia matrizeko j eta k errenkadak trukatu **/
        for (s=0; s<n; s++)
        {
            lag=aia[j][s];
            aia[j][s]=aia[k][s];
            aia[k][s]=lag;
        }
        /** aib bektoreko j eta k osagaiak trukatu **/
        lag=aib[j];
        aib[j]=aib[k];
        aib[k]=lag;
    }
}

    /**
    ** aia matrizearen LU faktorizazioa
    **/

    L.clear();
    L.resize(n, vector<complex<double> >(n, 0.0+0.0*I));
    U.clear();
    U.resize(n, vector<complex<double> >(n, 0.0+0.0*I));

    y.clear();
    y.resize(n, 0.0+0.0*I);

    for(j=0; j<n; j++)
    {
        for(i=0; i<n; i++)
        {
            if(i<=j)
            {
                U[i][j]=aia[i][j];
                for(k=0; k+1<=i; k++) /** k unsigned delako egiten dut k+1<=i hau **/
                {
                    U[i][j]=U[i][j]-L[i][k]*U[k][j];
                }
                if(i==j)
                {
                    L[i][j]=1;
                }
                else
                {
                    L[i][j]=0;
                }
            } /** i<=j begiztaren amaiera **/
            else
            { /** i>j begiztaren hasiera **/
                L[i][j]=aia[i][j];
                for(k=0; k+1<=j; k++) /** k unsigned delako egiten dut k+1<=j hau **/
                {
                    L[i][j]=L[i][j]-L[i][k]*U[k][j];
                }
                L[i][j]=L[i][j]/U[j][j];
                U[i][j]=0;
            } /** i>j baldintzaren amaiera **/

        } /** errenkaden i begiztaren amaiera **/

    } /** zutabeen j begiztaren amaiera **/

    /**
    ** Aurreranzko ebazpena y lortzeko
    ** L*y=aib
    **/

    for (i=0; i < n; i++)
    {
    y[i]=aib[i];
    for (j=0; j < i; j++)
        {
            y[i]=y[i]-L[i][j]*y[j];
        }
    }

    /**
    ** Atzeranzko ebazpena ebazpena lortzeko
    ** u*ebazpena=y
    **/

    for (i=n-1; i+1 >0; i--) /** i unsigned delako egiten dut i+1 >0 hau **/
    {
        ebazpena[i]=y[i];
        for (j=n-1; j > i; j--)
        {
            ebazpena[i]=ebazpena[i]-U[i][j]*ebazpena[j];
        }
        ebazpena[i]=ebazpena[i]/U[i][i];
    }

    hondarkalk(amat, bbek, zutabeak, ebazpena, errorea);

} /** kaminlu.cpp-ren amaiera **/

/*****************************************************************************/

void mirabiderm (const vector<vector<complex<double> > > &amat, vector<vector<complex<double> > > &aia)
/**
** Inputs: amat
** Output: aia
** [amat(1:m,1:n)]T * [amat(1:m,1:n)] = aia[n,n]
**/
{
    unsigned int m=amat.size();
    unsigned int n=amat[0].size();
    unsigned int i, j, s;
    vector<complex<double> > blag;
    complex<double> lag;

    for (i=0; i<n; i++)
    {
        blag.clear();
        for (j=0; j<n; j++)
        {
            lag=0.0+0.0*I;
            for (s=0; s < m; s++)
                {
                   lag=lag+amat[s][i]*amat[s][j];
                }
                blag.push_back(lag);
        }
        aia.push_back(blag);
    }
}

void mirabiderb (const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<complex<double> > &aib)
/**
** Inputs: amat, bbek
** Output: aib
** aib[neurria] = [amat(1:m,1:n)]T * [b(1:m,1)]
**/
{
    unsigned int m=amat.size();
    unsigned int n=amat[0].size();
    unsigned int i, s;
    complex<double> dk;

    for (i=0; i<n; i++)
    {
        dk=0.0+0.0*I;
        for (s=0; s < m; s++)
        {
            dk=dk+amat[s][i]*bbek[s];
        }
        aib.push_back(dk);
    }
}

void hondarkalk (const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, const vector<complex<double> > &ebazpena, float &errorea)
{
    /**
    ** Inputs: amat, bbek, ebazpena
    ** Output: errbek, hondarra
    ** bh= [amat(1:m,1:n)] * [ebazpena(1:n)]
    **/
    unsigned int m=amat.size();
    unsigned int n=amat[0].size();
    unsigned int i, j;
    complex<double> dk, b1, b2;

    vector<complex<double> > hondarra; /** Hondar-bektorea **/
    vector<complex<double> > bbekb; /** Gai askeen bektore aldatua **/

    hondarra.clear();

    for (i=0; i < m; i++)
    {
        dk=0.0+0.0*I;
        b1=0.0+0.0*I;
        b2=0.0+0.0*I;
        if (abs(bbek[i]) < 1.e-25)
        {
            for (j = zutabeak[0][0]; j <= zutabeak[0][1]; j++)
            {
                b1=b1+amat[i][j]*ebazpena[j];
            }

            for (j = zutabeak[1][0]; j <= zutabeak[1][1]; j++)
            {
                b2=b2+amat[i][j]*ebazpena[j];
            }

            hondarra.push_back(b1+b2);
            bbekb.push_back(b1);
        }
        else
        {
            for (j=0; j < n; j++)
            {
                dk=dk+amat[i][j]*ebazpena[j];
            }
            hondarra.push_back(bbek[i]-dk);
            bbekb.push_back(bbek[i]);
        }
    }
    errorea=100*bnorm(hondarra)/bnorm(bbekb);
}

double bnorm (const vector<complex<double> > &b)
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

    dot=sqrt(abs(dk));
    return dot;
}

