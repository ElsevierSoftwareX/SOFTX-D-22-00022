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

void ffinserie(const vector<double> &zentro, const vector<double> &xpuntu, const double k, const int nhmoz,  vector<vector<complex<double> > > &serie, vector<vector<complex<double> > > &dserie)
{
complex<double> cbal; /** serieko gaien deribatuak txertatzeko aldagai iragankorra **/
vector<complex<double> > sbek /** deribatuak esferikotan **/, cbek; /** deribatuak cartesiarretan **/

double koordc[3];
double r, theta, fi, kos, ksk;
vector<complex<double> > hankelfun1;
vector<complex<double> > dhankelfun1;
vector<complex<double> > hankelfun2;
vector<complex<double> > dhankelfun2;
vector<double> hannorm;
vector<vector<double>> matc(3, vector<double> (3, 0.0));

double lp, lpp, dlp;
double knorm, rnorm=0.8;
int kont, m, n;

r=0;

for (kont=0; kont<3; kont++)
{
    r=r+(xpuntu[kont]-zentro[kont])*(xpuntu[kont]-zentro[kont]);
    koordc[kont]=(xpuntu[kont]-zentro[kont]);
}
r=sqrt(r);

if (r<0.00001)
{
    r=0.00001;
}

/** 0 ordenako Hankel funtzioak **/

hannorm.push_back(1.0);
//hannorm.push_back(abs(-I*exp(I*k*rnorm)/(k*rnorm)));
cbal= (1/hannorm[0])*-I*exp(I*k*r)/(k*r);
hankelfun1.push_back(cbal);
cbal= (1/hannorm[0])*I*exp(-I*k*r)/(k*r);
hankelfun2.push_back(cbal);

/** Seriearen kalkulua **/

theta=acos(koordc[2]/r);

if (fabs(cos(theta)-1)<0.001) /** Ipar-lurmuturra x ardatzeko esferikoak **/
{
    kos=1.0;
    ksk=0.0;
    fi=0.0;

    matc[0][0]=0.0;
    matc[0][1]=1.0/koordc[2];
    matc[0][2]=0.0;
    matc[1][0]=0.0;
    matc[1][1]=0.0;
    matc[1][2]=1.0/koordc[2];
    matc[2][0]=1.0;
    matc[2][1]=0.0;
    matc[2][2]=0.0;
}
else if (fabs(cos(theta)+1)<0.001) /** Hego-lurmuturra x ardatzeko esferikoak **/
{
    kos=-1.0;
    ksk=0.0;
    fi=0.0;

    matc[0][0]=0.0;
    matc[0][1]=-1.0/koordc[2];
    matc[0][2]=0.0;
    matc[1][0]=0.0;
    matc[1][1]=0.0;
    matc[1][2]=-1.0/koordc[2];
    matc[2][0]=-1.0;
    matc[2][1]=0.0;
    matc[2][2]=0.0;
}
else /** z ardatzeko esferikoak **/
{
    kos=koordc[2]/r;
    ksk=1/sin(theta);
    fi=atan2(koordc[1],koordc[0]);

    matc[0][0]=koordc[0]/r;
    matc[0][1]=(koordc[0]*koordc[2])/(pow(r,2)*sqrt(pow(koordc[0],2)+pow(koordc[1],2)));
    matc[0][2]=-koordc[1]/(pow(koordc[0],2)+pow(koordc[1],2));
    matc[1][0]=koordc[1]/r;
    matc[1][1]=(koordc[1]*koordc[2])/(pow(r,2)*sqrt(pow(koordc[0],2)+pow(koordc[1],2)));
    matc[1][2]=koordc[0]/(pow(koordc[0],2)+pow(koordc[1],2));
    matc[2][0]=koordc[2]/r;
    matc[2][1]=-sqrt(pow(koordc[0],2)+pow(koordc[1],2))/pow(r,2);
    matc[2][2]=0;

} /** end of if kos **/

for (m=0; m<=nhmoz; m++)
    {

    /** m+1 ordenako Hankel funtzioak **/

    hannorm.push_back(1.0);
    //hannorm.push_back(abs(gsl_sf_bessel_jl(m+1+0.5, k*rnorm)+ I*gsl_sf_bessel_yl(m+1+0.5, k*rnorm)));
    cbal=(1/hannorm[m+1])*(gsl_sf_bessel_jl(m+1+0.5, k*r)+ I*gsl_sf_bessel_yl(m+1+0.5, k*r));
    hankelfun1.push_back(cbal);
    cbal=(1/hannorm[m+1])*(gsl_sf_bessel_jl(m+1+0.5, k*r)- I*gsl_sf_bessel_yl(m+1+0.5, k*r));
    hankelfun2.push_back(cbal);

    /** m ordenako Hankel funtzioen deribatuak **/
    cbal=(m/r)*hankelfun1[m]-k*hankelfun1[m+1]*(hannorm[m+1]/hannorm[m]);
    dhankelfun1.push_back(cbal);
    cbal=(m/r)*hankelfun2[m]-k*hankelfun2[m+1]*(hannorm[m+1]/hannorm[m]);
    dhankelfun2.push_back(cbal);

    /** n=0 kasua **/

    /** Pm Legendre Polinomioak **/

    lp=gsl_sf_legendre_Pl(m, kos);
    lpp=gsl_sf_legendre_Pl(m+1, kos);
    dlp=-(1+m)*ksk*(kos*lp-lpp);

    knorm=sqrt((2*m+1)/(4*PI));

    /** am0**/

    sbek.clear();
    sbek.push_back(knorm*hankelfun1[m]*lp); /** p **/
    serie.push_back(sbek);

    sbek.clear();
    sbek.push_back(knorm*hankelfun2[m]*lp); /** p **/
    serie.push_back(sbek);

    sbek.clear();
    sbek.push_back(knorm*dhankelfun1[m]*lp); /** p_r**/
    sbek.push_back(knorm*hankelfun1[m]*dlp); /** p_theta**/
    sbek.push_back(0.0+0.0*I);               /** p_fi**/

    cbek.clear();
    cbek.push_back(sbek[0]*matc[0][0]+sbek[1]*matc[0][1]+sbek[2]*matc[0][2]); /** p_x **/
    cbek.push_back(sbek[0]*matc[1][0]+sbek[1]*matc[1][1]+sbek[2]*matc[1][2]); /** p_y **/
    cbek.push_back(sbek[0]*matc[2][0]+sbek[1]*matc[2][1]+sbek[2]*matc[2][2]); /** p_z **/
    dserie.push_back(cbek);

    sbek.clear();
    sbek.push_back(knorm*dhankelfun2[m]*lp); /** p_r**/
    sbek.push_back(knorm*hankelfun2[m]*dlp); /** p_theta**/
    sbek.push_back(0.0+0.0*I);               /** p_fi**/

    cbek.clear();
    cbek.push_back(sbek[0]*matc[0][0]+sbek[1]*matc[0][1]+sbek[2]*matc[0][2]); /** p_x **/
    cbek.push_back(sbek[0]*matc[1][0]+sbek[1]*matc[1][1]+sbek[2]*matc[1][2]); /** p_y **/
    cbek.push_back(sbek[0]*matc[2][0]+sbek[1]*matc[2][1]+sbek[2]*matc[2][2]); /** p_z **/
    dserie.push_back(cbek);

    for (n=1; n<=m; n++)
        {

        /** Pmn Legendre Polinomioak **/

        lp=gsl_sf_legendre_Plm(m, n, kos);
        lpp=gsl_sf_legendre_Plm(m+1, n, kos);
        dlp=-(1+m)*kos*ksk*lp+(1+m-n)*ksk*lpp;

        /** Normalizazio konstantea **/

        knorm=sqrt(((2*m+1)*faktorial(m-n))/(2*PI*faktorial(m+n)));

        /** amn**/

        sbek.clear();
        sbek.push_back(knorm*hankelfun1[m]*lp*cos(n*fi)); /** p **/
        serie.push_back(sbek);

        sbek.clear();
        sbek.push_back(knorm*hankelfun2[m]*lp*cos(n*fi)); /** p **/
        serie.push_back(sbek);

        sbek.clear();
        sbek.push_back(knorm*dhankelfun1[m]*lp*cos(n*fi));      /** p_r**/
        sbek.push_back(knorm*hankelfun1[m]*dlp*cos(n*fi));      /** p_theta**/
        sbek.push_back(knorm*hankelfun1[m]*lp*(-n)*sin(n*fi));  /** p_fi**/

        cbek.clear();
        cbek.push_back(sbek[0]*matc[0][0]+sbek[1]*matc[0][1]+sbek[2]*matc[0][2]); /** p_x **/
        cbek.push_back(sbek[0]*matc[1][0]+sbek[1]*matc[1][1]+sbek[2]*matc[1][2]); /** p_y **/
        cbek.push_back(sbek[0]*matc[2][0]+sbek[1]*matc[2][1]+sbek[2]*matc[2][2]); /** p_z **/
        dserie.push_back(cbek);

        sbek.clear();
        sbek.push_back(knorm*dhankelfun2[m]*lp*cos(n*fi));      /** p_r**/
        sbek.push_back(knorm*hankelfun2[m]*dlp*cos(n*fi));      /** p_theta**/
        sbek.push_back(knorm*hankelfun2[m]*lp*(-n)*sin(n*fi));  /** p_fi**/

        cbek.clear();
        cbek.push_back(sbek[0]*matc[0][0]+sbek[1]*matc[0][1]+sbek[2]*matc[0][2]); /** p_x **/
        cbek.push_back(sbek[0]*matc[1][0]+sbek[1]*matc[1][1]+sbek[2]*matc[1][2]); /** p_y **/
        cbek.push_back(sbek[0]*matc[2][0]+sbek[1]*matc[2][1]+sbek[2]*matc[2][2]); /** p_z **/
        dserie.push_back(cbek);

        /** bmn**/

        sbek.clear();
        sbek.push_back(knorm*hankelfun1[m]*lp*sin(n*fi)); /** p **/
        serie.push_back(sbek);

        sbek.clear();
        sbek.push_back(knorm*hankelfun2[m]*lp*sin(n*fi)); /** p **/
        serie.push_back(sbek);

        sbek.clear();
        sbek.push_back(knorm*dhankelfun1[m]*lp*sin(n*fi));  /** p_r**/
        sbek.push_back(knorm*hankelfun1[m]*dlp*sin(n*fi));  /** p_theta**/
        sbek.push_back(knorm*hankelfun1[m]*lp*n*cos(n*fi)); /** p_fi**/

        cbek.clear();
        cbek.push_back(sbek[0]*matc[0][0]+sbek[1]*matc[0][1]+sbek[2]*matc[0][2]); /** p_x **/
        cbek.push_back(sbek[0]*matc[1][0]+sbek[1]*matc[1][1]+sbek[2]*matc[1][2]); /** p_y **/
        cbek.push_back(sbek[0]*matc[2][0]+sbek[1]*matc[2][1]+sbek[2]*matc[2][2]); /** p_z **/
        dserie.push_back(cbek);

        sbek.clear();
        sbek.push_back(knorm*dhankelfun2[m]*lp*sin(n*fi));  /** p_r**/
        sbek.push_back(knorm*hankelfun2[m]*dlp*sin(n*fi));  /** p_theta**/
        sbek.push_back(knorm*hankelfun2[m]*lp*n*cos(n*fi)); /** p_fi**/

        cbek.clear();
        cbek.push_back(sbek[0]*matc[0][0]+sbek[1]*matc[0][1]+sbek[2]*matc[0][2]); /** p_x **/
        cbek.push_back(sbek[0]*matc[1][0]+sbek[1]*matc[1][1]+sbek[2]*matc[1][2]); /** p_y **/
        cbek.push_back(sbek[0]*matc[2][0]+sbek[1]*matc[2][1]+sbek[2]*matc[2][2]); /** p_z **/
        dserie.push_back(cbek);

        } /**  end of for n=1:m **/

    } /** end of for m=0;nhmoz **/

} /** end of function ffinserie **/
