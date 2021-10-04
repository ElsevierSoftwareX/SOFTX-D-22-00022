#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include "elasserie.h"
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

void elasserie(const vector<double> &zentro, const vector<double> &xpuntu, const double &kp, const double &ks, const double &lambda, const double &mu, const int nhmoz, vector<vector<complex<double> > > &serie, vector<vector<complex<double> > > &dserie, vector<vector<complex<double> > > &tserie)
{

complex<double> cbal; /** serieko gaien deribatuak txertatzeko aldagai iragankorra **/
vector<complex<double> > sbek1, sbek2; /** lehenengo deribatuak esferikotan **/
vector<complex<double> > ssbek1, ssbek2; /** bigarren deribatuak esferikotan **/

double koordc[3];
double r, theta, fi, kos, ksk, kot;

/** Hankel funtzioak P uhinerako (gradientetik datorren irrotational wave) **/

vector<complex<double> > hankelPw1;
vector<complex<double> > dhankelPw1;
vector<complex<double> > d2hankelPw1;
vector<complex<double> > hankelPw2;
vector<complex<double> > dhankelPw2;
vector<complex<double> > d2hankelPw2;

vector<double> hannormP;

/** Hankel funtzioak S uhinerako (errotazionaletik datorren equivoluminal wave) **/

vector<complex<double> > hankelSw1;
vector<complex<double> > dhankelSw1;
vector<complex<double> > d2hankelSw1;
vector<complex<double> > hankelSw2;
vector<complex<double> > dhankelSw2;
vector<complex<double> > d2hankelSw2;

vector<double> hannormS;

vector<vector<double>> matc(3, vector<double> (3, 0.0));
vector<vector<vector<double>>> bimat(3, vector<vector<double>>(3, vector<double>(3, 0.0)));

double lp, lpp, lppp, dlp, d2lp; /** Legendre polinomioa eta deribatuak **/
double knorm, rnorm=0.1;
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

/** 0 ordenako P uhineko (kp uhin-zenbakia) Hankel funtzioak **/

hannormP.push_back(abs(-I*exp(I*kp*rnorm)/(kp*rnorm)));
cbal= (1/hannormP[0])*-I*exp(I*kp*r)/(kp*r);
hankelPw1.push_back(cbal);
cbal= (1/hannormP[0])*I*exp(-I*kp*r)/(kp*r);
hankelPw2.push_back(cbal);

/** 0 ordenako S uhineko (ks uhin-zenbakia) Hankel funtzioak **/

hannormS.push_back(abs(-I*exp(I*ks*rnorm)/(ks*rnorm)));
cbal= (1/hannormS[0])*-I*exp(I*ks*r)/(ks*r);
hankelSw1.push_back(cbal);
cbal= (1/hannormS[0])*I*exp(-I*ks*r)/(ks*r);
hankelSw2.push_back(cbal);

/** Seriearen kalkulua **/

theta=acos(koordc[2]/r);

if (fabs(cos(theta)-1)<0.0001) /** Ipar-lurmuturra **/
{
    kos=1.0;
    ksk=0.0;
    kot=kos*ksk;
    fi=0.0;

    matc[0][0]=0;
    matc[0][1]=0;
    matc[0][2]=1;
    matc[1][0]=0;
    matc[1][1]=1;
    matc[1][2]=0;
    matc[2][0]=1;
    matc[2][1]=0;
    matc[2][2]=0;

    bimat[0][0][0]=0;
    bimat[0][0][1]=0;
    bimat[0][0][2]=0;
    bimat[0][1][0]=0;
    bimat[0][1][1]=0;
    bimat[0][1][2]=0;
    bimat[0][2][0]=0;
    bimat[0][2][1]=0;
    bimat[0][2][2]=0;

    bimat[1][0][0]=0;
    bimat[1][0][1]=0;
    bimat[1][0][2]=0;
    bimat[1][1][0]=0;
    bimat[1][1][1]=0;
    bimat[1][1][2]=0;
    bimat[1][2][0]=0;
    bimat[1][2][1]=0;
    bimat[1][2][2]=0;

    bimat[2][0][0]=0;
    bimat[2][0][1]=0;
    bimat[2][0][2]=0;
    bimat[2][1][0]=0;
    bimat[2][1][1]=0;
    bimat[2][1][2]=0;
    bimat[2][2][0]=0;
    bimat[2][2][1]=0;
    bimat[2][2][2]=0;
}
else if (fabs(cos(theta)+1)<0.0001) //Hego-lurburua
{

    kos=-1.0;
    ksk=0.0;
    kot=kos*ksk;
    fi=0.0;

    matc[0][0]=0;
    matc[0][1]=0;
    matc[0][2]=0;
    matc[1][0]=0;
    matc[1][1]=0;
    matc[1][2]=0;
    matc[2][0]=-1;
    matc[1][1]=0;
    matc[1][2]=0;

    bimat[0][0][0]=0;
    bimat[0][0][1]=0;
    bimat[0][0][2]=0;
    bimat[0][1][0]=0;
    bimat[0][1][1]=0;
    bimat[0][1][2]=0;
    bimat[0][2][0]=0;
    bimat[0][2][1]=0;
    bimat[0][2][2]=0;

    bimat[1][0][0]=0;
    bimat[1][0][1]=0;
    bimat[1][0][2]=0;
    bimat[1][1][0]=0;
    bimat[1][1][1]=0;
    bimat[1][1][2]=0;
    bimat[1][2][0]=0;
    bimat[1][2][1]=0;
    bimat[1][2][2]=0;

    bimat[2][0][0]=0;
    bimat[2][0][1]=0;
    bimat[2][0][2]=0;
    bimat[2][1][0]=0;
    bimat[2][1][1]=0;
    bimat[2][1][2]=0;
    bimat[2][2][0]=0;
    bimat[2][2][1]=0;
    bimat[2][2][2]=0;
}
else
{
    kos=koordc[2]/r;
    ksk=1/sin(theta);
    kot=kos*ksk;
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

    bimat[0][0][0]=(pow(koordc[1],2)+pow(koordc[2],2))/pow(r,3);
    bimat[0][0][1]=-koordc[0]*koordc[1]/pow(r,3);
    bimat[0][0][2]=-koordc[0]*koordc[2]/pow(r,3);
    bimat[0][1][0]=bimat[0][0][1];
    bimat[0][1][1]=(pow(koordc[0],2)+pow(koordc[2],2))/pow(r,3);
    bimat[0][1][2]=-koordc[1]*koordc[2]/pow(r,3);
    bimat[0][2][0]=bimat[0][0][2];
    bimat[0][2][1]=bimat[0][1][2];
    bimat[0][2][2]=(pow(koordc[0],2)+pow(koordc[1],2))/pow(r,3);

    bimat[1][0][0]=koordc[2]*(-2*pow(koordc[0],4)-pow(koordc[0],2)*pow(koordc[1],2)+pow(koordc[1],4)+pow(koordc[1],2)*pow(koordc[2],2))/(pow(pow(koordc[0],2)+pow(koordc[1],2),1.5)*pow(r,4));
    bimat[1][0][1]=-koordc[0]*koordc[1]*koordc[2]*(3*(pow(koordc[0],2)+pow(koordc[1],2))+pow(koordc[2],2))/(pow(pow(koordc[0],2)+pow(koordc[1],2),1.5)*pow(r,4));
    bimat[1][0][2]=koordc[0]*(pow(koordc[0],2)+pow(koordc[1],2)-pow(koordc[2],2))/(pow(pow(koordc[0],2)+pow(koordc[1],2),0.5)*pow(r,4));
    bimat[1][1][0]=bimat[1][0][1];
    bimat[1][1][1]=koordc[2]*(pow(koordc[0],4)-2*pow(koordc[1],4)+pow(koordc[0],2)*(pow(koordc[2],2)-pow(koordc[1],2)))/(pow(pow(koordc[0],2)+pow(koordc[1],2),1.5)*pow(r,4));
    bimat[1][1][2]=koordc[1]*(pow(koordc[0],2)+pow(koordc[1],2)-pow(koordc[2],2))/(pow(pow(koordc[0],2)+pow(koordc[1],2),0.5)*pow(r,4));
    bimat[1][2][0]=bimat[1][0][2];
    bimat[1][2][1]=bimat[1][1][2];
    bimat[1][2][2]=2*koordc[2]*pow(pow(koordc[0],2)+pow(koordc[1],2),0.5)/pow(r,4);

    bimat[2][0][0]=2*koordc[0]*koordc[1]/pow(pow(koordc[0],2)+pow(koordc[1],2),2);
    bimat[2][0][1]=(-pow(koordc[0],2)+pow(koordc[1],2))/pow(pow(koordc[0],2)+pow(koordc[1],2),2);
    bimat[2][0][2]=0.0;
    bimat[2][1][0]= bimat[2][0][1];
    bimat[2][1][1]= -bimat[2][0][0];
    bimat[2][1][2]=0.0;
    bimat[2][2][0]=0.0;
    bimat[2][2][1]=0.0;
    bimat[2][2][2]=0.0;

} /** end of if kos **/

/***************************************************************************************************************************************/
/**************************************************** Seriearen kalkulua ***************************************************************/
/***************************************************************************************************************************************/

for (m=0; m<=nhmoz; m++)
    {
    /** Pm0 Legendre Polinomioak eta deribatuak **/

    /** m+1 ordenako P uhineko (kp uhin-zenbakia) Hankel funtzioak **/

    hannormP.push_back(abs(gsl_sf_bessel_jl(m+1+0.5, kp*rnorm)+ I*gsl_sf_bessel_yl(m+1+0.5, kp*rnorm)));
    cbal=(1/hannormP[m+1])*(gsl_sf_bessel_jl(m+1+0.5, kp*r)+ I*gsl_sf_bessel_yl(m+1+0.5, kp*r));
    hankelPw1.push_back(cbal);
    cbal=(1/hannormP[m+1])*(gsl_sf_bessel_jl(m+1+0.5, kp*r)- I*gsl_sf_bessel_yl(m+1+0.5, kp*r));
    hankelPw2.push_back(cbal);

    /** m ordenako P uhineko (kp uhin-zenbakia) Hankel funtzioen deribatuak **/

    cbal=(m/r)*hankelPw1[m]-kp*hankelPw1[m+1]*(hannormP[m+1]/hannormP[m]);
    dhankelPw1.push_back(cbal);
    cbal=(m/r)*hankelPw2[m]-kp*hankelPw2[m+1]*(hannormP[m+1]/hannormP[m]);
    dhankelPw2.push_back(cbal);

    /** m ordenako P uhineko (kp uhin-zenbakia) Hankel funtzioen bigarren deribatuak **/

    cbal=(((-1+m)*m-pow(kp*r,2))*hankelPw1[m]+2*kp*r*hankelPw1[m+1]*(hannormP[m+1]/hannormP[m]))/pow(r,2);
    d2hankelPw1.push_back(cbal);
    cbal=(((-1+m)*m-pow(kp*r,2))*hankelPw2[m]+2*kp*r*hankelPw2[m+1]*(hannormP[m+1]/hannormP[m]))/pow(r,2);
    d2hankelPw2.push_back(cbal);

    /** m+1 ordenako S uhineko (ks uhin-zenbakia) Hankel funtzioak **/

    hannormS.push_back(abs(gsl_sf_bessel_jl(m+1+0.5, ks*rnorm)+ I*gsl_sf_bessel_yl(m+1+0.5, ks*rnorm)));
    cbal=(1/hannormS[m+1])*(gsl_sf_bessel_jl(m+1+0.5, ks*r)+ I*gsl_sf_bessel_yl(m+1+0.5, ks*r));
    hankelSw1.push_back(cbal);
    cbal=(1/hannormS[m+1])*(gsl_sf_bessel_jl(m+1+0.5, ks*r)- I*gsl_sf_bessel_yl(m+1+0.5, ks*r));
    hankelSw2.push_back(cbal);

    /** m ordenako S uhineko (ks uhin-zenbakia) Hankel funtzioen deribatuak **/

    cbal=(m/r)*hankelSw1[m]-ks*hankelSw1[m+1]*(hannormS[m+1]/hannormS[m]);
    dhankelSw1.push_back(cbal);
    cbal=(m/r)*hankelSw2[m]-ks*hankelSw2[m+1]*(hannormS[m+1]/hannormS[m]);
    dhankelSw2.push_back(cbal);

    /** m ordenako S uhineko (kp uhin-zenbakia) Hankel funtzioen bigarren deribatuak **/

    cbal=(((-1+m)*m-pow(ks*r,2))*hankelSw1[m]+2*ks*r*hankelSw1[m+1]*(hannormS[m+1]/hannormS[m]))/pow(r,2);
    d2hankelSw1.push_back(cbal);
    cbal=(((-1+m)*m-pow(ks*r,2))*hankelSw2[m]+2*ks*r*hankelSw2[m+1]*(hannormS[m+1]/hannormS[m]))/pow(r,2);
    d2hankelSw2.push_back(cbal);

    /**************************************************************************************************************************/
    /************************************************************** n=0 kasua *************************************************/
    /**************************************************************************************************************************/

    /** Pm Legendre Polinomioak eta deribatuak **/
    lp=gsl_sf_legendre_Pl(m, kos);
    lpp=gsl_sf_legendre_Pl(m+1, kos);
    lppp=gsl_sf_legendre_Pl(m+2, kos);
    dlp=-(1+m)*ksk*(kos*lp-lpp);
    d2lp=0.5*(1+m)*pow(ksk,2)*((3+m+(1+m)*cos(2*theta))*lp-2*(2+m)*(2*kos*lpp-lppp));

    /** Normalizazio konstantea **/

    knorm=sqrt((2*m+1)/(4*PI));

    /************************************************************** a_m0^0 *************************************************/

    /******************************************************* Desplazamenduak *************************************************/

    sbek1.clear();
    /** P uhineko 1 motako Hankel funtzioarekin **/
    sbek1.push_back(knorm*dhankelPw1[m]*lp);        /** u_r **/
    sbek1.push_back(knorm*(1/r)*hankelPw1[m]*dlp);  /** u_theta **/
    sbek1.push_back(0.0+0.0*I);                     /** u_fi **/

    sbek2.clear();
    /** P uhineko 1 motako Hankel funtzioarekin **/
    sbek2.push_back(knorm*dhankelPw2[m]*lp);        /** u_r **/
    sbek2.push_back(knorm*(1/r)*hankelPw2[m]*dlp);  /** u_theta **/
    sbek2.push_back(0.0+0.0*I);                     /** u_fi **/

    esfcarserie(sbek1, sbek2, matc, serie);

    /******************************************************* Deformazioak eta tentsioak *************************************************/

    ssbek1.clear();
    /** P uhineko 1 motako Hankel funtzioarekin **/
    ssbek1.push_back(knorm*d2hankelPw1[m]*lp);                                      /** 0: u_r,r **/
    ssbek1.push_back(knorm*dhankelPw1[m]*dlp);                                      /** 1: u_r,theta**/
    ssbek1.push_back(0.0+0.0*I);                                                    /** 2: u_r,fi**/
    ssbek1.push_back(knorm*((1/pow(r,2))*(-hankelPw1[m]+r*dhankelPw1[m])*dlp));     /** 3: u_theta,r **/
    ssbek1.push_back(knorm*(1/r)*hankelPw1[m]*d2lp);                                /** 4: u_theta,theta**/
    ssbek1.push_back(0.0+0.0*I);                                                    /** 5: u_theta,fi**/
    ssbek1.push_back(0.0+0.0*I);                                                    /** 6: u_fi,r **/
    ssbek1.push_back(0.0+0.0*I);                                                    /** 7: u_fi,theta**/
    ssbek1.push_back(knorm*(-(ksk/r))*hankelPw1[m]*lp*pow(n,2));                    /** 8: u_fi,fi**/

    ssbek2.clear();
    /** P uhineko 2 motako Hankel funtzioarekin **/
    ssbek2.push_back(knorm*d2hankelPw2[m]*lp);                                      /** 0: u_r,r **/
    ssbek2.push_back(knorm*dhankelPw2[m]*dlp);                                      /** 1: u_r,theta**/
    ssbek2.push_back(0.0+0.0*I);                                                    /** 2: u_r,fi**/
    ssbek2.push_back(knorm*((1/pow(r,2))*(-hankelPw2[m]+r*dhankelPw2[m])*dlp));     /** 3: u_theta,r **/
    ssbek2.push_back(knorm*(1/r)*hankelPw2[m]*d2lp);                                /** 4: u_theta,theta**/
    ssbek2.push_back(0.0+0.0*I);                                                    /** 5: u_theta,fi**/
    ssbek2.push_back(0.0+0.0*I);                                                    /** 6: u_fi,r **/
    ssbek2.push_back(0.0+0.0*I);                                                    /** 7: u_fi,theta**/
    ssbek2.push_back(knorm*(-(ksk/r))*hankelPw2[m]*lp*pow(n,2));                    /** 8: u_fi,fi**/

    esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

    /************************************************************** b_m0^0 denak 0   *************************************************/

    /************************************************************** a_m0^r *************************************************/

    /******************************************************* Desplazamenduak *************************************************/

    sbek1.clear();
    /** S uhineko 1 motako Hankel funtzioarekin **/
    sbek1.push_back(0.0+0.0*I);                     /** u_r **/
    sbek1.push_back(0.0+0.0*I);                     /** u_theta **/
    sbek1.push_back(knorm*(-1/r)*hankelSw1[m]*dlp); /** u_fi **/

    sbek2.clear();
    /** S uhineko 2 motako Hankel funtzioarekin **/
    sbek2.push_back(0.0+0.0*I);                     /** u_r **/
    sbek2.push_back(0.0+0.0*I);                     /** u_theta **/
    sbek2.push_back(knorm*(-1/r)*hankelSw2[m]*dlp); /** u_fi **/

    esfcarserie(sbek1, sbek2, matc, serie);

    /******************************************************* Deformazioak eta tentsioak *************************************************/

    ssbek1.clear();
    /** S uhineko 1 motako Hankel funtzioarekin **/
    ssbek1.push_back(0.0+0.0*I);                                                /** 0: u_r,r **/
    ssbek1.push_back(0.0+0.0*I);                                                /** 1: u_r,theta**/
    ssbek1.push_back(0.0+0.0*I);                                                /** 2: u_r,fi**/
    ssbek1.push_back(0.0+0.0*I);                                                /** 3: u_theta,r **/
    ssbek1.push_back(0.0+0.0*I);                                                /** 4: u_theta,theta**/
    ssbek1.push_back(knorm*(-ksk/r)*hankelSw1[m]*lp*pow(n,2));                  /** 5: u_theta,fi**/
    ssbek1.push_back(knorm*((1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m])*dlp));  /** 6: u_fi,r **/
    ssbek1.push_back(knorm*(-1/r)*hankelSw1[m]*d2lp);                           /** 7: u_fi,theta**/
    ssbek1.push_back(0.0+0.0*I);                                                /** 8: u_fi,fi**/

    ssbek2.clear();
    /** S uhineko 2 motako Hankel funtzioarekin **/
    ssbek2.push_back(0.0+0.0*I);                                                /** 0: u_r,r **/
    ssbek2.push_back(0.0+0.0*I);                                                /** 1: u_r,theta**/
    ssbek2.push_back(0.0+0.0*I);                                                /** 2: u_r,fi**/
    ssbek2.push_back(0.0+0.0*I);                                                /** 3: u_theta,r **/
    ssbek2.push_back(0.0+0.0*I);                                                /** 4: u_theta,theta**/
    ssbek2.push_back(knorm*(-ksk/r)*hankelSw2[m]*lp*pow(n,2));                  /** 5: u_theta,fi**/
    ssbek2.push_back(knorm*((1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m])*dlp));  /** 6: u_fi,r **/
    ssbek2.push_back(knorm*(-1/r)*hankelSw2[m]*d2lp);                           /** 7: u_fi,theta**/
    ssbek2.push_back(0.0+0.0*I);                                                /** 8: u_fi,fi**/

    esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

    /************************************************************** b_m0^r denak 0 *************************************************/

    /************************************************************** a_m0^theta *************************************************/

    /******************************************************* Desplazamenduak *************************************************/

    sbek1.clear();
    /** S uhineko 1 motako Hankel funtzioarekin **/
    sbek1.push_back(0.0+0.0*I);                                     /** u_r **/
    sbek1.push_back(0.0+0.0*I);                                     /** u_theta **/
    sbek1.push_back(knorm*(1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp); /** u_fi **/

    sbek2.clear();
    /** S uhineko 2 motako Hankel funtzioarekin **/
    sbek2.push_back(0.0+0.0*I);                                     /** u_r **/
    sbek2.push_back(0.0+0.0*I);                                     /** u_theta **/
    sbek2.push_back(knorm*(1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp); /** u_fi **/

    esfcarserie(sbek1, sbek2, matc, serie);

    /******************************************************* Deformazioak eta tentsioak *************************************************/

    ssbek1.clear();
    /** S uhineko 1 motako Hankel funtzioarekin **/
    ssbek1.push_back(0.0+0.0*I);                                                                        /** 0: u_r,r **/
    ssbek1.push_back(0.0+0.0*I);                                                                        /** 1: u_r,theta **/
    ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*lp*pow(n,2));                                           /** 2: u_r,fi **/
    ssbek1.push_back(0.0+0.0*I);                                                                        /** 3: u_theta,r **/
    ssbek1.push_back(0.0+0.0*I);                                                                        /** 4: u_theta,theta **/
    ssbek1.push_back(0.0+0.0*I);                                                                        /** 5: u_theta,fi **/
    ssbek1.push_back(knorm*(1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m]+pow(r,2)*d2hankelSw1[m])*lp);    /** 6: u_fi,r **/
    ssbek1.push_back(knorm*((1/r)*(hankelSw1[m]+r*dhankelSw1[m])*dlp));                                 /** 7: u_fi,theta **/
    ssbek1.push_back(0.0+0.0*I);                                                                        /** 8: u_fi,fi **/

    ssbek2.clear();
    /** S uhineko 2 motako Hankel funtzioarekin **/
    ssbek2.push_back(0.0+0.0*I);                                                                        /** 0: u_r,r **/
    ssbek2.push_back(0.0+0.0*I);                                                                        /** 1: u_r,theta **/
    ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*lp*pow(n,2));                                           /** 2: u_r,fi **/
    ssbek2.push_back(0.0+0.0*I);                                                                        /** 3: u_theta,r **/
    ssbek2.push_back(0.0+0.0*I);                                                                        /** 4: u_theta,theta **/
    ssbek2.push_back(0.0+0.0*I);                                                                        /** 5: u_theta,fi **/
    ssbek2.push_back(knorm*(1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m]+pow(r,2)*d2hankelSw2[m])*lp);    /** 6: u_fi,r **/
    ssbek2.push_back(knorm*((1/r)*(hankelSw2[m]+r*dhankelSw2[m])*dlp));                                 /** 7: u_fi,theta **/
    ssbek2.push_back(0.0+0.0*I);                                                                        /** 8: u_fi,fi **/

    esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

    /************************************************************** b_m0^theta denak 0 *************************************************/

    /************************************************************** a_m0^fi *************************************************/

    /******************************************************* Desplazamenduak *************************************************/

    sbek1.clear();
    /** S uhineko 1 motako Hankel funtzioarekin **/
    sbek1.push_back(knorm*(1/r)*hankelSw1[m]*(kot*lp+dlp));           /** u_r **/
    sbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp);  /** u_theta **/
    sbek1.push_back(0.0+0.0*I);                                       /** u_fi **/

    sbek2.clear();
    /** S uhineko 2 motako Hankel funtzioarekin **/
    sbek2.push_back(knorm*(1/r)*hankelSw2[m]*(kot*lp+dlp));           /** u_r **/
    sbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp);  /** u_theta **/
    sbek2.push_back(0.0+0.0*I);                                       /** u_fi **/

    esfcarserie(sbek1, sbek2, matc, serie);

    /******************************************************* Deformazioak eta tentsioak *************************************************/

    ssbek1.clear();
    /** S uhineko 1 motako Hankel funtzioarekin **/
    ssbek1.push_back(knorm*(1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m])*(kot*lp+dlp));              /** 0: u_r,r **/
    ssbek1.push_back(knorm*(1/r)*hankelSw1[m]*(-ksk*ksk*lp+kot*dlp+d2lp));                          /** 1: u_r,theta **/
    ssbek1.push_back(0.0+0.0*I);                                                                    /** 2: u_r,fi **/
    ssbek1.push_back(knorm*(1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m]-pow(r,2)*d2hankelSw1[m])*lp); /** 3: u_theta,r **/
    ssbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*dlp);                              /** 4: u_theta,theta **/
    ssbek1.push_back(0.0+0.0*I);                                                                    /** 5: u_theta,fi **/
    ssbek1.push_back(0.0+0.0*I);                                                                    /** 6: u_fi,r **/
    ssbek1.push_back(0.0+0.0*I);                                                                    /** 7: u_fi,theta **/
    ssbek1.push_back(0.0+0.0*I);                                                                    /** 8: u_fi,fi **/

    ssbek2.clear();
    /** S uhineko 2 motako Hankel funtzioarekin **/
    ssbek2.push_back(knorm*(1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m])*(kot*lp+dlp));              /** 0: u_r,r **/
    ssbek1.push_back(knorm*(1/r)*hankelSw2[m]*(-ksk*ksk*lp+kot*dlp+d2lp));                          /** 1: u_r,theta **/
    ssbek2.push_back(0.0+0.0*I);                                                                    /** 2. u_r,fi **/
    ssbek2.push_back(knorm*(1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m]-pow(r,2)*d2hankelSw2[m])*lp); /** 3: u_theta,r **/
    ssbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*dlp);                              /** 4: u_theta,theta **/
    ssbek2.push_back(0.0+0.0*I);                                                                    /** 5: u_theta,fi **/
    ssbek2.push_back(0.0+0.0*I);                                                                    /** 6: u_fi,r **/
    ssbek2.push_back(0.0+0.0*I);                                                                    /** 7: u_fi,theta **/
    ssbek2.push_back(0.0+0.0*I);                                                                    /** 8: u_fi,fi**/

    esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

    /************************************************************** b_m0^fi denak 0 *************************************************/

    for (n=1; n<=m; n++)
        {
        /** Pm Legendre Polinomioak eta deribatuak **/
        lp=gsl_sf_legendre_Plm(m, n, kos);
        lpp=gsl_sf_legendre_Plm(m+1, n, kos);
        lppp=gsl_sf_legendre_Plm(m+2, n, kos);
        dlp=-(1+m)*(ksk*cos(theta))*lp+(1+m-n)*ksk*lpp;
        d2lp=0.5*pow(ksk,2)*((1+m)*((3+m+(1+m)*cos(2*theta)))*lp+2*(1+m-n)*(-2*(2+m)*kos*lpp+(2+m-n)*lppp));

        /** Normalizazio konstantea **/

        knorm=sqrt(((2*m+1)*faktorial(m-n))/(2*PI*faktorial(m+n)));

        /************************************************************** a_mn^0 *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** P uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(knorm*dhankelPw1[m]*lp*cos(n*fi));              /** u_r **/
        sbek1.push_back(knorm*(1/r)*hankelPw1[m]*dlp*cos(n*fi));        /** u_theta **/
        sbek1.push_back(knorm*(-ksk/r)*hankelPw1[m]*lp*n*sin(n*fi));    /** u_fi **/

        sbek2.clear();
        /** P uhineko 1 motako Hankel funtzioarekin **/
        sbek2.push_back(knorm*dhankelPw2[m]*lp*cos(n*fi));              /** u_r **/
        sbek2.push_back(knorm*(1/r)*hankelPw2[m]*dlp*cos(n*fi));        /** u_theta **/
        sbek2.push_back(knorm*(-ksk/r)*hankelPw2[m]*lp*n*sin(n*fi));    /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** P uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(knorm*d2hankelPw1[m]*lp*cos(n*fi));                                    /** 0: u_r,r **/
        ssbek1.push_back(knorm*dhankelPw1[m]*dlp*cos(n*fi));                                    /** 1: u_r,theta**/
        ssbek1.push_back(knorm*dhankelPw1[m]*lp*(-n)*sin(n*fi));                                /** 2: u_r,fi**/
        ssbek1.push_back(knorm*((1/pow(r,2))*(-hankelPw1[m]+r*dhankelPw1[m])*dlp*cos(n*fi)));   /** 3: u_theta,r **/
        ssbek1.push_back(knorm*(1/r)*hankelPw1[m]*d2lp*cos(n*fi));                              /** 4: u_theta,theta**/
        ssbek1.push_back(knorm*(-1/r)*hankelPw1[m]*dlp*n*sin(n*fi));                            /** 5: u_theta,fi**/
        ssbek1.push_back(knorm*((ksk/pow(r,2))*(hankelPw1[m]-r*dhankelPw1[m])*lp*n*sin(n*fi))); /** 6: u_fi,r **/
        ssbek1.push_back(knorm*((ksk/r)*hankelPw1[m]*(kot*lp-dlp)*n*sin(n*fi)));                /** 7: u_fi,theta**/
        ssbek1.push_back(knorm*(-(ksk/r))*hankelPw1[m]*lp*pow(n,2)*cos(n*fi));                  /** 8: u_fi,fi**/

        ssbek2.clear();
        /** P uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(knorm*d2hankelPw2[m]*lp*cos(n*fi));                                    /** 0: u_r,r **/
        ssbek2.push_back(knorm*dhankelPw2[m]*dlp*cos(n*fi));                                    /** 1: u_r,theta**/
        ssbek2.push_back(knorm*dhankelPw2[m]*lp*(-n)*sin(n*fi));                                /** 2: u_r,fi**/
        ssbek2.push_back(knorm*((1/pow(r,2))*(-hankelPw2[m]+r*dhankelPw2[m])*dlp*cos(n*fi)));   /** 3: u_theta,r **/
        ssbek2.push_back(knorm*(1/r)*hankelPw2[m]*d2lp*cos(n*fi));                              /** 4: u_theta,theta**/
        ssbek2.push_back(knorm*(-1/r)*hankelPw2[m]*dlp*n*sin(n*fi));                            /** 5: u_theta,fi**/
        ssbek2.push_back(knorm*((ksk/pow(r,2))*(hankelPw2[m]-r*dhankelPw2[m])*lp*n*sin(n*fi))); /** 6: u_fi,r **/
        ssbek2.push_back(knorm*((ksk/r)*hankelPw2[m]*(kot*lp-dlp)*n*sin(n*fi)));                /** 7: u_fi,theta**/
        ssbek2.push_back(knorm*(-(ksk/r))*hankelPw2[m]*lp*pow(n,2)*cos(n*fi));                  /** 8: u_fi,fi**/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** b_mn^0 *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** P uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(knorm*dhankelPw1[m]*lp*sin(n*fi));          /** u_r **/
        sbek1.push_back(knorm*(1/r)*hankelPw1[m]*dlp*sin(n*fi));    /** u_theta **/
        sbek1.push_back(knorm*(ksk/r)*hankelPw1[m]*lp*n*cos(n*fi)); /** u_fi **/

        sbek2.clear();
        /** P uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(knorm*dhankelPw2[m]*lp*sin(n*fi));          /** u_r **/
        sbek2.push_back(knorm*(1/r)*hankelPw2[m]*dlp*sin(n*fi));    /** u_theta **/
        sbek2.push_back(knorm*(ksk/r)*hankelPw2[m]*lp*n*cos(n*fi)); /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** P uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(knorm*d2hankelPw1[m]*lp*sin(n*fi));                                                /** 0: u_r,r **/
        ssbek1.push_back(knorm*dhankelPw1[m]*dlp*sin(n*fi));                                                /** 1: u_r,theta**/
        ssbek1.push_back(knorm*dhankelPw1[m]*lp*n*cos(n*fi));                                               /** 2: u_r,fi**/
        ssbek1.push_back(knorm*((1/pow(r,2))*(-hankelPw1[m]+r*dhankelPw1[m])*dlp*sin(n*fi)));               /** 3: u_theta,r **/
        ssbek1.push_back(knorm*(1/r)*hankelPw1[m]*d2lp*sin(n*fi));                                          /** 4: u_theta,theta**/
        ssbek1.push_back(knorm*(1/r)*hankelPw1[m]*dlp*n*cos(n*fi));                                         /** 5: u_theta,fi**/
        ssbek1.push_back(knorm*((ksk/pow(r,2))*(-hankelPw1[m]+r*dhankelPw1[m])*lp*n*cos(n*fi)));            /** 6: u_fi,r **/
        ssbek1.push_back(knorm*((ksk/r)*hankelPw1[m]*(-kot*lp+dlp)*n*cos(n*fi)));                           /** 7: u_fi,theta**/
        ssbek1.push_back(knorm*((-ksk/r)*hankelPw1[m]*lp*pow(n,2)*sin(n*fi)));                              /** 8: u_fi,fi**/

        ssbek2.clear();
        /** P uhineko 1 motako Hankel funtzioarekin **/
        ssbek2.push_back(knorm*d2hankelPw2[m]*lp*sin(n*fi));                                                /** 0: u_r,r **/
        ssbek2.push_back(knorm*dhankelPw2[m]*dlp*sin(n*fi));                                                /** 1: u_r,theta**/
        ssbek2.push_back(knorm*dhankelPw2[m]*lp*n*cos(n*fi));                                               /** 2: u_r,fi**/
        ssbek2.push_back(knorm*((1/pow(r,2))*(-hankelPw2[m]+r*dhankelPw2[m])*dlp*sin(n*fi)));               /** 3: u_theta,r **/
        ssbek2.push_back(knorm*(1/r)*hankelPw2[m]*d2lp*sin(n*fi));                                          /** 4: u_theta,theta**/
        ssbek2.push_back(knorm*(1/r)*hankelPw2[m]*dlp*n*cos(n*fi));                                         /** 5: u_theta,fi**/
        ssbek2.push_back(knorm*((ksk/pow(r,2))*(-hankelPw2[m]+r*dhankelPw2[m])*lp*n*cos(n*fi)));            /** 6: u_fi,r **/
        ssbek2.push_back(knorm*((ksk/r)*hankelPw2[m]*(-kot*lp+dlp)*n*cos(n*fi)));                           /** 7: u_fi,theta**/
        ssbek2.push_back(knorm*((-ksk/r)*hankelPw2[m]*lp*pow(n,2)*sin(n*fi)));                              /** 8: u_fi,fi**/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** a_mn^r *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(0.0+0.0*I);                                     /** u_r **/
        sbek1.push_back(knorm*(-ksk/r)*hankelSw1[m]*lp*n*sin(n*fi));    /** u_theta **/
        sbek1.push_back(knorm*(-1/r)*hankelSw1[m]*dlp*cos(n*fi));       /** u_fi **/

        sbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(0.0+0.0*I);                                     /** u_r **/
        sbek2.push_back(knorm*(-ksk/r)*hankelSw2[m]*lp*n*sin(n*fi));    /** u_theta **/
        sbek2.push_back(knorm*(-1/r)*hankelSw2[m]*dlp*cos(n*fi));       /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(0.0+0.0*I);                                                                        /** 0: u_r,r **/
        ssbek1.push_back(0.0+0.0*I);                                                                        /** 1: u_r,theta**/
        ssbek1.push_back(0.0+0.0*I);                                                                        /** 2: u_r,fi**/
        ssbek1.push_back(knorm*ksk*(1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m])*lp*n*sin(n*fi));             /** 3: u_theta,r **/
        ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*((kot*lp-dlp)*n*sin(n*fi)));                            /** 4: u_theta,theta**/
        ssbek1.push_back(knorm*(-ksk/r)*hankelSw1[m]*lp*pow(n,2)*cos(n*fi));                                /** 5: u_theta,fi**/
        ssbek1.push_back(knorm*((1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m])*dlp*cos(n*fi)));                /** 6: u_fi,r **/
        ssbek1.push_back(knorm*(-1/r)*hankelSw1[m]*d2lp*cos(n*fi));                                         /** 7: u_fi,theta**/
        ssbek1.push_back(knorm*(1/r)*hankelSw1[m]*dlp*n*sin(n*fi));                                         /** 8: u_fi,fi**/

        ssbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(0.0+0.0*I);                                                                        /** 0: u_r,r **/
        ssbek2.push_back(0.0+0.0*I);                                                                        /** 1: u_r,theta**/
        ssbek2.push_back(0.0+0.0*I);                                                                        /** 2: u_r,fi**/
        ssbek2.push_back(knorm*ksk*(1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m])*lp*n*sin(n*fi));             /** 3: u_theta,r **/
        ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*((kot*lp-dlp)*n*sin(n*fi)));                            /** 4: u_theta,theta**/
        ssbek2.push_back(knorm*(-ksk/r)*hankelSw2[m]*lp*pow(n,2)*cos(n*fi));                                /** 5: u_theta,fi**/
        ssbek2.push_back(knorm*((1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m])*dlp*cos(n*fi)));                /** 6: u_fi,r **/
        ssbek2.push_back(knorm*(-1/r)*hankelSw2[m]*d2lp*cos(n*fi));                                         /** 7: u_fi,theta**/
        ssbek2.push_back(knorm*(1/r)*hankelSw2[m]*dlp*n*sin(n*fi));                                         /** 8: u_fi,fi**/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** b_mn^r *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(0.0+0.0*I);                                     /** u_r **/
        sbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*lp*n*cos(n*fi));     /** u_theta **/
        sbek1.push_back(knorm*(-1/r)*hankelSw1[m]*dlp*sin(n*fi));       /** u_fi **/

        sbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(0.0+0.0*I);                                     /** u_r **/
        sbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*lp*n*cos(n*fi));     /** u_theta **/
        sbek2.push_back(knorm*(-1/r)*hankelSw2[m]*dlp*sin(n*fi));       /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(0.0+0.0*I);                                                                        /** 0: u_r,r **/
        ssbek1.push_back(0.0+0.0*I);                                                                        /** 1: u_r,theta**/
        ssbek1.push_back(0.0+0.0*I);                                                                        /** 2: u_r,fi**/
        ssbek1.push_back(knorm*ksk*((1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m])*lp*n*cos(n*fi)));          /** 3: u_theta,r **/
        ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*(-kot*lp+dlp)*n*cos(n*fi));                             /** 4: u_theta,theta**/
        ssbek1.push_back(knorm*(-ksk/r)*hankelSw1[m]*lp*pow(n,2)*sin(n*fi));                                /** 5: u_theta,fi**/
        ssbek1.push_back(knorm*((1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m])*dlp*sin(n*fi)));                /** 6: u_fi,r **/
        ssbek1.push_back(knorm*(-1/r)*hankelSw1[m]*d2lp*sin(n*fi));                                         /** 7: u_fi,theta**/
        ssbek1.push_back(knorm*(-1/r)*hankelSw1[m]*dlp*n*cos(n*fi));                                        /** 8: u_fi,fi**/

        ssbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(0.0+0.0*I);                                                                        /** 0: u_r,r **/
        ssbek2.push_back(0.0+0.0*I);                                                                        /** 1: u_r,theta**/
        ssbek2.push_back(0.0+0.0*I);                                                                        /** 2: u_r,fi**/
        ssbek2.push_back(knorm*ksk*((1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m])*lp*n*cos(n*fi)));          /** 3: u_theta,r **/
        ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*(-kot*lp+dlp)*n*cos(n*fi));                             /** 4: u_theta,theta**/
        ssbek2.push_back(knorm*(-ksk/r)*hankelSw2[m]*lp*pow(n,2)*sin(n*fi));                                /** 5: u_theta,fi**/
        ssbek2.push_back(knorm*((1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m])*dlp*sin(n*fi)));                /** 6: u_fi,r **/
        ssbek2.push_back(knorm*(-1/r)*hankelSw2[m]*d2lp*sin(n*fi));                                         /** 7: u_fi,theta**/
        ssbek2.push_back(knorm*(-1/r)*hankelSw2[m]*dlp*n*cos(n*fi));                                        /** 8: u_fi,fi**/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** a_mn^theta *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*lp*n*sin(n*fi));                 /** u_r **/
        sbek1.push_back(0.0+0.0*I);                                                 /** u_theta **/
        sbek1.push_back(knorm*(1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*cos(n*fi));   /** u_fi **/

        sbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*lp*n*sin(n*fi));                 /** u_r **/
        sbek2.push_back(0.0+0.0*I);                                                 /** u_theta **/
        sbek2.push_back(knorm*(1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*cos(n*fi));   /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(knorm*(ksk/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m])*lp*n*sin(n*fi));                      /** 0: u_r,r **/
        ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*(-kot*lp+dlp)*n*sin(n*fi));                                     /** 1: u_r,theta **/
        ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*lp*pow(n,2)*cos(n*fi));                                         /** 2: u_r,fi **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 3: u_theta,r **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 4: u_theta,theta **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 5: u_theta,fi **/
        ssbek1.push_back(knorm*(1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m]+pow(r,2)*d2hankelSw1[m])*lp*cos(n*fi));  /** 6: u_fi,r **/
        ssbek1.push_back(knorm*((1/r)*(hankelSw1[m]+r*dhankelSw1[m])*dlp*cos(n*fi)));                               /** 7: u_fi,theta **/
        ssbek1.push_back(knorm*((-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*n*sin(n*fi)));                             /** 8: u_fi,fi **/

        ssbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(knorm*(ksk/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m])*lp*n*sin(n*fi));                      /** 0: u_r,r **/
        ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*(-kot*lp+dlp)*n*sin(n*fi));                                     /** 1: u_r,theta **/
        ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*lp*pow(n,2)*cos(n*fi));                                         /** 2: u_r,fi **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 3: u_theta,r **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 4: u_theta,theta **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 5: u_theta,fi **/
        ssbek2.push_back(knorm*(1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m]+pow(r,2)*d2hankelSw2[m])*lp*cos(n*fi));  /** 6: u_fi,r **/
        ssbek2.push_back(knorm*((1/r)*(hankelSw2[m]+r*dhankelSw2[m])*dlp*cos(n*fi)));                               /** 7: u_fi,theta **/
        ssbek2.push_back(knorm*((-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*n*sin(n*fi)));                             /** 8: u_fi,fi **/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** b_mn^theta *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(-knorm*(ksk/r)*hankelSw1[m]*lp*n*cos(n*fi));                /** u_r **/
        sbek1.push_back(0.0+0.0*I);                                                 /** u_theta **/
        sbek1.push_back(knorm*(1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*sin(n*fi));   /** u_fi **/

        sbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(-knorm*(ksk/r)*hankelSw2[m]*lp*n*cos(n*fi));                /** u_r **/
        sbek2.push_back(0.0+0.0*I);                                                 /** u_theta **/
        sbek2.push_back(knorm*(1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*sin(n*fi));   /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(knorm*(ksk/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m])*lp*n*cos(n*fi));                       /** 0: u_r,r **/
        ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*(kot*lp-dlp)*n*cos(n*fi));                                      /** 1: u_r,theta **/
        ssbek1.push_back(knorm*(ksk/r)*hankelSw1[m]*lp*pow(n,2)*sin(n*fi));                                         /** 2: u_r,fi **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 3: u_theta,r **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 4: u_theta,theta **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 5: u_theta,fi **/
        ssbek1.push_back(knorm*(1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m]+pow(r,2)*d2hankelSw1[m])*lp*sin(n*fi));  /** 6: u_fi,r **/
        ssbek1.push_back(knorm*((1/r)*(hankelSw1[m]+r*dhankelSw1[m]))*dlp*sin(n*fi));                               /** 7: u_fi,theta **/
        ssbek1.push_back(knorm*((1/r)*hankelSw1[m]+r*dhankelSw1[m])*lp*n*cos(n*fi));                                /** 8: u_fi,fi **/

        ssbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(knorm*(ksk/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m])*lp*n*cos(n*fi));                       /** 0: u_r,r **/
        ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*(kot*lp-dlp)*n*cos(n*fi));                                      /** 1: u_r,theta **/
        ssbek2.push_back(knorm*(ksk/r)*hankelSw2[m]*lp*pow(n,2)*sin(n*fi));                                         /** 2: u_r,fi **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 3: u_theta,r **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 4: u_theta,theta **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 5: u_theta,fi **/
        ssbek2.push_back(knorm*(1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m]+pow(r,2)*d2hankelSw2[m])*lp*cos(n*fi));  /** 6: u_fi,r **/
        ssbek2.push_back(knorm*((1/r)*(hankelSw2[m]+r*dhankelSw2[m]))*dlp*sin(n*fi));                               /** 7: u_fi,theta **/
        ssbek2.push_back(knorm*((1/r)*hankelSw2[m]+r*dhankelSw2[m])*lp*n*cos(n*fi));                                /** 8: u_fi,fi **/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** a_mn^fi *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(knorm*(1/r)*hankelSw1[m]*(kot*lp+dlp)*cos(n*fi));           /** u_r **/
        sbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*cos(n*fi));  /** u_theta **/
        sbek1.push_back(0.0+0.0*I);                                                 /** u_fi **/

        sbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(knorm*(1/r)*hankelSw2[m]*(kot*lp+dlp)*cos(n*fi));           /** u_r **/
        sbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*cos(n*fi));  /** u_theta **/
        sbek2.push_back(0.0+0.0*I);                                                 /** u_fi **/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(knorm*(1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m])*(kot*lp+dlp)*cos(n*fi));                /** 0: u_r,r **/
        ssbek1.push_back(knorm*(1/r)*hankelSw1[m]*(-ksk*ksk*lp+kot*dlp+d2lp)*cos(n*fi));                            /** 1: u_r,theta **/
        ssbek1.push_back(knorm*(-1/r)*hankelSw1[m]*(kot*lp+dlp)*n*sin(n*fi));                                       /** 2: u_r,fi **/
        ssbek1.push_back(knorm*(1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m]-pow(r,2)*d2hankelSw1[m])*lp*cos(n*fi));   /** 3: u_theta,r **/
        ssbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*dlp*cos(n*fi));                                /** 4: u_theta,theta **/
        ssbek1.push_back(knorm*(1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*n*sin(n*fi));                                /** 5: u_theta,fi **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 6: u_fi,r **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 7: u_fi,theta **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 8: u_fi,fi **/

        ssbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(knorm*(1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m])*(kot*lp+dlp)*cos(n*fi));                /** 0: u_r,r **/
        ssbek1.push_back(knorm*(1/r)*hankelSw2[m]*(-ksk*ksk*lp+kot*dlp+d2lp)*cos(n*fi));                            /** 1: u_r,theta **/
        ssbek2.push_back(knorm*(-1/r)*hankelSw2[m]*(kot*lp+dlp)*n*sin(n*fi));                                       /** 2. u_r,fi **/
        ssbek2.push_back(knorm*(1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m]-pow(r,2)*d2hankelSw2[m])*lp*cos(n*fi));   /** 3: u_theta,r **/
        ssbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*dlp*cos(n*fi));                                /** 4: u_theta,theta **/
        ssbek2.push_back(knorm*(1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*n*sin(n*fi));                                /** 5: u_theta,fi **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 6: u_fi,r **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 7: u_fi,theta **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 8: u_fi,fi**/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        /************************************************************** b_mn^fi *************************************************/

        /******************************************************* Desplazamenduak *************************************************/

        sbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        sbek1.push_back(knorm*(1/r)*hankelSw1[m]*(kot*lp+dlp)*sin(n*fi));           /** u_r **/
        sbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*sin(n*fi));  /** u_theta **/
        sbek1.push_back(0.0+0.0*I);                                                 /** u_fi **/

        sbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        sbek2.push_back(knorm*(1/r)*hankelSw2[m]*(kot*lp+dlp)*sin(n*fi));           /** u_r **/
        sbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*sin(n*fi));  /** u_theta**/
        sbek2.push_back(0.0+0.0*I);                                                 /** u_fi**/

        esfcarserie(sbek1, sbek2, matc, serie);

        /******************************************************* Deformazioak eta tentsioak *************************************************/

        ssbek1.clear();
        /** S uhineko 1 motako Hankel funtzioarekin **/
        ssbek1.push_back(knorm*(1/pow(r,2))*(-hankelSw1[m]+r*dhankelSw1[m])*(kot*lp+dlp)*sin(n*fi));                /** 0: u_r,r **/
        ssbek1.push_back(knorm*(1/r)*hankelSw1[m]*(-ksk*ksk*lp+kot*dlp+d2lp)*sin(n*fi));                            /** 1: u_r,theta **/
        ssbek1.push_back(knorm*(1/r)*hankelSw1[m]*(kot*lp+dlp)*n*cos(n*fi));                                        /** 2: u_r,fi **/
        ssbek1.push_back(knorm*(1/pow(r,2))*(hankelSw1[m]-r*dhankelSw1[m]-pow(r,2)*d2hankelSw1[m])*lp*sin(n*fi));   /** 3: u_theta,r **/
        ssbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*dlp*sin(n*fi));                                /** 4: u_theta,theta **/
        ssbek1.push_back(knorm*(-1/r)*(hankelSw1[m]+r*dhankelSw1[m])*lp*n*cos(n*fi));                               /** 5: u_theta,fi **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 6: u_fi,r **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 7: u_fi,theta **/
        ssbek1.push_back(0.0+0.0*I);                                                                                /** 8: u_fi,fi **/

        ssbek2.clear();
        /** S uhineko 2 motako Hankel funtzioarekin **/
        ssbek2.push_back(knorm*(1/pow(r,2))*(-hankelSw2[m]+r*dhankelSw2[m])*(kot*lp+dlp)*sin(n*fi));                /** 0: u_r,r **/
        ssbek1.push_back(knorm*(1/r)*hankelSw2[m]*(-ksk*ksk*lp+kot*dlp+d2lp)*sin(n*fi));                            /** 1: u_r,theta **/
        ssbek2.push_back(knorm*(1/r)*hankelSw2[m]*(kot*lp+dlp)*n*cos(n*fi));                                        /** 2. u_r,fi **/
        ssbek2.push_back(knorm*(1/pow(r,2))*(hankelSw2[m]-r*dhankelSw2[m]-pow(r,2)*d2hankelSw2[m])*lp*sin(n*fi));   /** 3: u_theta,r **/
        ssbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*dlp*sin(n*fi));                                /** 4: u_theta,theta **/
        ssbek2.push_back(knorm*(-1/r)*(hankelSw2[m]+r*dhankelSw2[m])*lp*n*cos(n*fi));                               /** 5: u_theta,fi **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 6: u_fi,r **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 7: u_fi,theta **/
        ssbek2.push_back(0.0+0.0*I);                                                                                /** 8: u_fi,fi**/

        esfcardserie(sbek1, sbek2, ssbek1, ssbek2, matc, bimat, lambda, mu, dserie, tserie);

        } /**  end of for n=1:m **/

    } /** end of for m=1;nhmoz **/

}  /** end of function elasserie **/

/**************************************************************************************************************************/
/**************************************************************************************************************************/
/**************************************************************************************************************************/

void esfcarserie(const vector<complex<double> > &sbek1, const vector<complex<double> > &sbek2, const vector<vector<double> > &matc, vector<vector<complex<double> > > &serie)
/**
** 1.(u_r, u_theta, u_fi) -> (u_x, u_y, u_z)
** 2. add to kserie
**/
{

    complex<double> cbal1, cbal2;
    vector<complex<double> > cbek1, cbek2;
    unsigned int i, k;

    /** Desplazamenduak **/
    /** 1 aurreranzko uhina **/
    /** 2 atzeranzko uhina **/

    cbek1.clear();
    cbek2.clear();
    for (i=0; i<3; i++)
    {
        cbal1=0.0+0.0*I;
        cbal2=0.0+0.0*I;
        for (k=0; k<3; k++)
        {
            cbal1=cbal1+sbek1[k]*matc[i][k];
            cbal2=cbal2+sbek2[k]*matc[i][k];
        }
        cbek1.push_back(cbal1);
        cbek2.push_back(cbal2);
    }

    //cout << setprecision(10) << cbek1[0] << "  " <<  cbek1[1] << "  " <<  cbek1[2] << endl;

    /** 1 aurreranzko uhina **/
    serie.push_back(cbek1);

    /** 2 atzeranzko uhina **/
    serie.push_back(cbek2);
}

void esfcardserie(const vector<complex<double> > &sbek1, const vector<complex<double> > &sbek2, const vector<complex<double> > &ssbek1, const vector<complex<double> > &ssbek2, const vector<vector<double> > &matc, const vector<vector<vector<double> > >&bimat, const double lambda, const double mu, vector<vector<complex<double> > > &dserie, vector<vector<complex<double> > > &tserie)
/**
 1.(u_r,r, u_r,theta, u_r,fi, u_theta,r, u_theta,theta, u_theta,fi, u_fi,r, u_fi,theta, u_fi,fi) ->
--> (u_xx, u_xy, u_xz, u_yx, u_yy, u_yz, u_zx, u_zy, u_zz)

** 2. dserie serieari gehitu

** 3.(u_xx, u_xy, u_xz, u_yx, u_yy, u_yz, u_zx, u_zy, u_zz)->
--> (s_xx, s_xy, s_xz, s_yx, s_yy, s_yz, s_zx, s_zy, s_zz)

** 4. tserie serieari gehitu

**/
{

    complex<double> cbal1, cbal2;
    vector<complex<double> > ebek1, ebek2, tbek;
    unsigned int i, j, k, l, kont;

    /** Deformazioak **/
    /** 1 aurreranzko uhina **/
    /** 2 atzeranzko uhina **/

    ebek1.clear();
    ebek2.clear();

    for (i=0; i<3; i++)
    {
        for (j=0; j<3; j++)
        {
            cbal1=0.0+0.0*I;
            cbal2=0.0+0.0*I;
            for (k=0; k<3; k++)
            {
                cbal1=cbal1+sbek1[k]*bimat[k][i][j];
                cbal2=cbal2+sbek2[k]*bimat[k][i][j];
            }
            kont=0;
            for (k=0; k<3; k++)
            {
                for (l=0; l<3; l++)
                {
                    cbal1=cbal1+ssbek1[kont]*matc[i][k]*matc[j][l];
                    cbal2=cbal2+ssbek2[kont]*matc[i][k]*matc[j][l];
                    kont++;
                }
            }
            ebek1.push_back(cbal1);
            ebek2.push_back(cbal2);
        }
    }

    dserie.push_back(ebek1);
    dserie.push_back(ebek2);

    /** Tentsioak **/
    /** 1 aurreranzko uhina **/

    tbek.clear();
    cbal1=ebek1[0]+ebek1[4]+ebek1[8];

    tbek.push_back(lambda*cbal1+2*mu*ebek1[0]);     /** 0: s^1_x,x **/
    tbek.push_back(mu*(ebek1[1]+ebek1[3]));         /** 1: s^1_x,y **/
    tbek.push_back(mu*(ebek1[2]+ebek1[6]));         /** 2: s^1_x,z **/
    tbek.push_back(mu*(ebek1[3]+ebek1[1]));         /** 3: s^1_y,x **/
    tbek.push_back(lambda*cbal1+2*mu*ebek1[4]);     /** 4: s^1_y,y **/
    tbek.push_back(mu*(ebek1[5]+ebek1[7]));         /** 5: s^1_y,z **/
    tbek.push_back(mu*(ebek1[6]+ebek1[2]));         /** 6: s^1_z,x **/
    tbek.push_back(mu*(ebek1[7]+ebek1[5]));         /** 7: s^1_y,x **/
    tbek.push_back(lambda*cbal1+2*mu*ebek1[8]);     /** 8: s^1_z,z **/

    tserie.push_back(tbek);

    /** Tentsioak **/
    /** 2 atzeranzko uhina **/

    tbek.clear();
    cbal2=ebek2[0]+ebek2[4]+ebek2[8];

    tbek.push_back(lambda*cbal1+2*mu*ebek2[0]);     /** 0: s^1_x,x **/
    tbek.push_back(mu*(ebek2[1]+ebek2[3]));         /** 1: s^1_x,y **/
    tbek.push_back(mu*(ebek2[2]+ebek2[6]));         /** 2: s^1_x,z **/
    tbek.push_back(mu*(ebek2[3]+ebek2[1]));         /** 3: s^1_y,x **/
    tbek.push_back(lambda*cbal1+2*mu*ebek2[4]);     /** 4: s^1_y,y **/
    tbek.push_back(mu*(ebek2[5]+ebek2[7]));         /** 5: s^1_y,z **/
    tbek.push_back(mu*(ebek2[6]+ebek2[2]));         /** 6: s^1_z,x **/
    tbek.push_back(mu*(ebek2[7]+ebek2[5]));         /** 7: s^1_y,x **/
    tbek.push_back(lambda*cbal1+2*mu*ebek2[8]);     /** 8: s^1_z,z **/

    tserie.push_back(tbek);
}
