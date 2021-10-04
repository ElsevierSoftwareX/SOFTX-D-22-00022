#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include "pcgnorm.h"

const complex<double> I{0.0, 1.0};

using namespace std;

void pcgnorm(const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, vector<complex<double> > &ebazpena, float &errorea)
{
/**
** Inputs: amat, bbek
** Output: ebazpena
** amat*ebazpena=bbek sistemaren Karratu Minimoren araberako ebazpen
** hurbildua "Preconditioned Conjugate Gradient (PCG) metodoaren bidez
    **/

vector<complex<double> > ck;

/**
** Matrizek eta bektoreak
**/
vector<vector<complex<double> > > aia;
vector<complex<double> > aib;
vector<complex<double> > d, r, q, aq, br;
complex<double> rgnorm, alpha, rnorm, beta;
double aibnorm;

/** Iterazioaren kontrola **/

vector<float> erroreak; /** Urrats guztietako erroreak **/
float abiadura, abiaduramax=0.05; /** azken urratseko errorearen aldaketaren abiadura **/
float erroremax=0.05;
erroreak.clear();
unsigned int urraskop=1;
float aurrekoerroreak[urraskop]; /** aurreko urratsetako erroreak **/

/** Indizeak **/

unsigned int n, i, iter, itermax;

n=amat[0].size();

aia.clear();
aib.clear();
ck.clear();
ebazpena.clear();

/**
** [amat(1:m,1:n)]T * [amat(1:m,1:n] = aia[n,n]
**/
mirabiderm (amat, aia);

/**
** [amat(1:m,1:n)]T * [b(1:m,1)] = aib[n]
**/
mirabiderb (amat, bbek, aib);
aibnorm=bnorm(aib);

for (i=0; i<n; i++)
{
    d.push_back((1.0+0.0*I)/aia[i][i]);
    //ebazpena.push_back(0.0+0.0*I);
    ebazpena.push_back(d[i]*aib[i]);
}

/**
** [amat(1:m,1:n)]T * [amat(1:m,1:n)] [ebazpena(1:n,1)] = ck[n]
**/
mbiderb (aia, ebazpena, ck);

for (i=0; i<n; i++)
{
    r.push_back(aib[i]-ck[i]);
    q.push_back(r[i]*d[i]);
}

rgnorm=dot(r,q);
errorea=100*bnorm(r)/aibnorm;
erroreak.push_back(errorea);

/**
** Iterazioen hasiera
**/

for (i=0; i < urraskop; i++)
{
    aurrekoerroreak[i]=50000.+i*1000000.;
}

abiadura=(aurrekoerroreak[urraskop-1]-errorea)/urraskop;
iter=0;

while (errorea > erroremax && abs(abiadura) > abiaduramax)
{
    iter++;
    aq.clear();
    br.clear();

    mbiderb(aia, q, aq);
    alpha=rgnorm/dot(q,aq);

    for (i=0; i<n; i++)
    {
        ebazpena[i]=ebazpena[i]+alpha*q[i];
        r[i]=r[i]-alpha*aq[i];
        br.push_back(r[i]*d[i]);
    }
    rnorm=dot(r,br);
    beta=rnorm/rgnorm;

    for (i=0; i<n; i++)
    {
        q[i]=br[i]+beta*q[i];
    }
    rgnorm=rnorm;
    errorea=100*bnorm(r)/aibnorm;

    erroreak.push_back(errorea);
    abiadura=(aurrekoerroreak[urraskop-1]-errorea)/urraskop;

    for (i=urraskop-1; i > 0; i--)
        {
        aurrekoerroreak[i]=aurrekoerroreak[i-1];
        }
    aurrekoerroreak[0]=errorea;

} /** Iterazioaren amaiera**/

cout << "1. PCG: errorea: % " << fixed << setprecision(3) << errorea << " ; abiadura= " << abiadura << "; iterazioak= " << iter << endl;
hondarkalk(amat, bbek, zutabeak, ebazpena, errorea);

/** Errorerik txikiena egin duen itermax iterazio kopuru metatua aukeratu eta berriteratu **/

if (iter>0)
{
    itermax=0;
    erroremax=1000.0;

    for (i=0;i<erroreak.size(); i++)
    {
        if (isnan(erroreak[i]))
        {
            erroreak[i]=1000.0;
        }
        else
        {
        if (erroreak[i] < erroremax)
        {
            itermax=i;
            erroremax=erroreak[i];
        }
        }
    }

    /** Berriteratu itermax bider **/

    erroreak.clear();
    d.clear();
    r.clear();
    q.clear();
    ck.clear();
    ebazpena.clear();

    for (i=0; i<n; i++)
    {
        d.push_back((1.0+0.0*I)/aia[i][i]);
        //ebazpena.push_back(0.0+0.0*I);
        ebazpena.push_back(d[i]*aib[i]);
    }

    /**
    ** [amat(1:m,1:n)]T * [amat(1:m,1:n)] [ebazpena(1:n,1)] = ck[n]
    **/
    mbiderb (aia, ebazpena, ck);

    for (i=0; i<n; i++)
    {
        r.push_back(aib[i]-ck[i]);
        q.push_back(r[i]*d[i]);
    }

    rgnorm=dot(r,q);
    errorea=100*bnorm(r)/aibnorm;
    erroreak.push_back(errorea);

    /**
    ** Iterazioen hasiera
    **/

    for (i=0; i < urraskop; i++)
    {
        aurrekoerroreak[i]=50000.+i*1000000.;
    }

    abiadura=(aurrekoerroreak[urraskop-1]-errorea)/urraskop;
    iter=0;

    while (iter < itermax)
    {
        iter++;
        aq.clear();
        br.clear();

        mbiderb(aia, q, aq);
        alpha=rgnorm/dot(q,aq);

        for (i=0; i<n; i++)
        {
            ebazpena[i]=ebazpena[i]+alpha*q[i];
            r[i]=r[i]-alpha*aq[i];
            br.push_back(r[i]*d[i]);
        }
        rnorm=dot(r,br);
        beta=rnorm/rgnorm;

        for (i=0; i<n; i++)
        {
            q[i]=br[i]+beta*q[i];
        }
        rgnorm=rnorm;
        errorea=100*bnorm(r)/aibnorm;

        erroreak.push_back(errorea);
        abiadura=(aurrekoerroreak[urraskop-1]-errorea)/urraskop;

        for (i=urraskop-1; i > 0; i--)
            {
            aurrekoerroreak[i]=aurrekoerroreak[i-1];
            }
        aurrekoerroreak[0]=errorea;

    } /** Berriterazioaren amaiera**/

    cout << "2. PCG: errorea: % " << fixed << setprecision(3) << errorea << " ; abiadura= " << abiadura << "; iterazioak= " << iter << endl;
    hondarkalk(amat, bbek, zutabeak, ebazpena, errorea);

} /** of if iter>0**/

} /** pcgnorm.cpp-ren amaiera **/

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

void mbiderb (const vector<vector<complex<double> > > &amat, const vector<complex<double> > &bbek, vector<complex<double> > &aib)
/**
** Inputs: amat, bbek
** Output: aib
** aib[neurria] = [amat(1:m,1:n)] * [b(1:n,1)]
**/
{
    unsigned int m=amat.size();
    unsigned int n=amat[0].size();
    unsigned int i, s;
    complex<double> dk;

    for (i=0; i<m; i++)
    {
        dk=0.0+0.0*I;
        for (s=0; s < n; s++)
        {
            dk=dk+amat[i][s]*bbek[s];
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

complex<double> dot (const vector<complex<double> > &a, const vector<complex<double> > &b)
/**
** Input: a, b
** Output: dot =||b||
**/
{
    complex<double> dot=0.0+0.0*I;
    unsigned int i;
    unsigned int k = a.size();

    for (i=0; i < k; i++)
    {
        dot=dot+a[i]*b[i];
    }
    return dot;
}
