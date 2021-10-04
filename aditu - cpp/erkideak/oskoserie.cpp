#include <vector>
#include <complex>
using namespace std;

const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

template <typename T, typename U>
inline complex<T> operator*(complex<T> lhs, const U& rhs)
{
    return lhs *= rhs;
}

void oskoserie(const vector<double> &zentro, const vector<double> &xpuntu, const vector<double> &normal, const double zo, const double kp, const double ks, const double young, const double poisson, const int nhmoz, vector<vector<complex<double> > > &serie, vector<vector<complex<double> > > &dserie, vector<vector<complex<double> > > &tserie)
{
double norm;
double theta, tartea;
double kp0, kp1;
double ks0, ks1;
int m;

vector<complex<double>> ebek, tbek; /** oskol-sistemako bektoreak **/
vector<complex<double>> bek; /** sistema orokorreko bektorea **/

vector<double> xos, xl, b1, b2;

/** D matrizea **/

vector<vector<double>> dmat(5, vector<double> (5, 0.0));

dmat[0][0]=young/(1-pow(poisson,2));
dmat[0][1]=poisson*young/(1-pow(poisson,2));
dmat[0][2]=0.0;
dmat[0][3]=0.0;
dmat[0][4]=0.0;

dmat[1][0]=poisson*young/(1-pow(poisson,2));
dmat[1][1]=young/(1-pow(poisson,2));
dmat[1][2]=0.0;
dmat[1][3]=0.0;
dmat[1][4]=0.0;

dmat[2][0]=0.0;
dmat[2][1]=0.0;
dmat[2][2]=0.5*(1-pow(poisson,2))*young/(1+poisson);
dmat[2][3]=0.0;
dmat[2][4]=0.0;

dmat[3][0]=0.0;
dmat[3][1]=0.0;
dmat[3][2]=0.0;
dmat[3][3]=0.5*(5/6)*young/(1+poisson);
dmat[3][4]=0.0;

dmat[4][0]=0.0;
dmat[4][1]=0.0;
dmat[4][2]=0.0;
dmat[4][3]=0.0;
dmat[4][4]=0.5*(5/6)*young/(1+poisson);

/** Transformazioko bektoreen kalkulua **/

if (abs(normal[0])>0)
{
    b1.push_back(-normal[1]/normal[0]);
    b1.push_back(1.0);
    b1.push_back(0.0);

    b2.push_back(normal[1]*b1[2]-normal[2]*b1[1]);
    b2.push_back(normal[2]*b1[0]-normal[0]*b1[2]);
    b2.push_back(normal[0]*b1[1]-normal[1]*b1[0]);
}

else if (abs(normal[1])>0)
{
    b2.push_back(1.0);
    b2.push_back(-normal[0]/normal[1]);
    b2.push_back(0.0);

    b1.push_back(normal[2]*b2[1]-normal[1]*b2[2]);
    b1.push_back(normal[0]*b2[2]-normal[2]*b2[0]);
    b1.push_back(normal[1]*b2[0]-normal[0]*b2[1]);
}
else
{

    b1.push_back(1.0);
    b1.push_back(0.0);
    b1.push_back(0.0);

    b2.push_back(0.0);
    b2.push_back(1.0);
    b2.push_back(0.0);
}

norm=sqrt(b1[0]*b1[0]+b1[1]*b1[1]+b1[2]*b1[2]);
b1[0]=b1[0]/norm;
b1[1]=b1[1]/norm;
b1[2]=b1[2]/norm;

norm=sqrt(b2[0]*b2[0]+b2[1]*b2[1]+b2[2]*b2[2]);
b2[0]=b2[0]/norm;
b2[1]=b2[1]/norm;
b2[2]=b2[2]/norm;

xos.push_back(xpuntu[0]-zo*normal[0]);
xos.push_back(xpuntu[1]-zo*normal[1]);
xos.push_back(xpuntu[2]-zo*normal[2]);

xos[0]=xos[0]-zentro[0];
xos[1]=xos[1]-zentro[1];
xos[2]=xos[2]-zentro[2];

xl.push_back(b1[0]*xos[0]+b1[1]*xos[1]+b1[2]*xos[2]);
xl.push_back(b2[0]*xos[0]+b2[1]*xos[1]+b2[2]*xos[2]);
xl.push_back(normal[0]*xos[0]+normal[1]*xos[1]+normal[2]*xos[2]);

/***************************************************************************************************************************************/
/**************************************************** Seriearen kalkulua ***************************************************************/
/***************************************************************************************************************************************/

serie.clear();
dserie.clear();
tserie.clear();

theta=2*PI/(nhmoz+2);
tartea=2*PI/(nhmoz+2);

for (m=0; m<=nhmoz; m++)
    {

    kp0=kp*cos(theta);
    kp1=kp*sin(theta);

    ks0=ks*cos(theta);
    ks1=ks*sin(theta);

    /******************************************************* Desplazamenduak *************************************************/

    /** P1 = sin(kp0*xl[0])*sin(kp1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(kp0*cos(kp0*xl[0])*sin(kp1*xl[1]));   /** ulx **/
    ebek.push_back(kp1*sin(kp0*xl[0])*cos(kp1*xl[1]));   /** uly **/
    ebek.push_back(sin(kp0*xl[0])*sin(kp1*xl[1]));       /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** P2 = sin(kp0*xl[0])*cos(kp1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(kp0*cos(kp0*xl[0])*cos(kp1*xl[1]));   /** ulx **/
    ebek.push_back(-kp1*sin(kp0*xl[0])*sin(kp1*xl[1]));  /** uly **/
    ebek.push_back(sin(kp0*xl[0])*cos(kp1*xl[1]));       /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** P3 = cos(kp0*xl[0])*sin(kp1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(-kp0*sin(kp0*xl[0])*sin(kp1*xl[1])); /** ulx **/
    ebek.push_back(kp1*cos(kp0*xl[0])*cos(kp1*xl[1]));  /** uly **/
    ebek.push_back(cos(kp0*xl[0])*sin(kp1*xl[1]));      /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** P4 = cos(kp0*xl[0])*cos(kp1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(-kp0*sin(kp0*xl[0])*cos(kp1*xl[1])); /** ulx **/
    ebek.push_back(-kp1*cos(kp0*xl[0])*sin(kp1*xl[1])); /** uly **/
    ebek.push_back(cos(kp0*xl[0])*cos(kp1*xl[1]));      /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** S1 = sin(ks0*xl[0])*sin(ks1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(zo*ks0*cos(ks0*xl[0]*sin(ks1*xl[1])));   /** ulx **/
    ebek.push_back(zo*ks1*sin(ks0*xl[0]*cos(ks1*xl[1])));   /** uly **/
    ebek.push_back(0.0);                                    /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** S2 = sin(ks0*xl[0])*cos(ks1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(zo*ks0*cos(ks0*xl[0]*cos(ks1*xl[1])));   /** ulx **/
    ebek.push_back(-zo*ks1*sin(ks0*xl[0]*sin(ks1*xl[1])));  /** uly **/
    ebek.push_back(0.0);                                    /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** S3 = cos(ks0*xl[0])*sin(ks1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(-zo*ks0*sin(ks0*xl[0]*sin(ks1*xl[1]))); /** ulx **/
    ebek.push_back(zo*ks1*cos(ks0*xl[0]*cos(ks1*xl[1])));  /** uly **/
    ebek.push_back(0.0);                                   /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /** S4 = cos(ks0*xl[0])*cos(ks1*xl[1]) **/

    /** oskol-sistemako desplazamenduak **/
    ebek.clear();
    ebek.push_back(-zo*ks0*sin(ks0*xl[0]*cos(ks1*xl[1]))); /** ulx **/
    ebek.push_back(-zo*ks1*cos(ks0*xl[0]*sin(ks1*xl[1]))); /** uly **/
    ebek.push_back(0.0);                                   /** ulz **/

    /** Sistema orokorreko desplazamenduak **/
    bek.clear();
    bek.push_back(b1[0]*ebek[0]+b2[0]*ebek[1]+normal[0]*ebek[2]);
    bek.push_back(b1[1]*ebek[0]+b2[1]*ebek[1]+normal[1]*ebek[2]);
    bek.push_back(b1[2]*ebek[0]+b2[2]*ebek[1]+normal[2]*ebek[2]);
    serie.push_back(bek);

    /******************************************************* Deformazioak eta tentsioak *************************************************/

    /** P1 = sin(kp0*xl[0])*sin(kp1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-pow(kp0,2)*sin(kp0*xl[0])*sin(kp1*xl[1]));  /** exx **/
    ebek.push_back(-pow(kp1,2)*sin(kp0*xl[0])*sin(kp1*xl[1]));  /** eyy **/
    ebek.push_back(2*kp0*kp1*cos(kp0*xl[0])*cos(kp1*xl[1]));    /** exy **/
    ebek.push_back(kp0*cos(kp0*xl[0])*sin(kp1*xl[1]));          /** exz **/
    ebek.push_back(kp1*sin(kp0*xl[0])*cos(kp1*xl[1]));          /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** P2 = sin(kp0*xl[0])*cos(kp1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-pow(kp0,2)*sin(kp0*xl[0])*cos(kp1*xl[1]));  /** exx **/
    ebek.push_back(-pow(kp1,2)*sin(kp0*xl[0])*cos(kp1*xl[1]));  /** eyy **/
    ebek.push_back(-2*kp0*kp1*cos(kp0*xl[0])*sin(kp1*xl[1]));   /** exy **/
    ebek.push_back(kp0*cos(kp0*xl[0])*cos(kp1*xl[1]));          /** exz **/
    ebek.push_back(-kp1*sin(kp0*xl[0])*sin(kp1*xl[1]));         /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** P3 = cos(kp0*xl[0])*sin(kp1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-pow(kp0,2)*cos(kp0*xl[0])*sin(kp1*xl[1]));  /** exx **/
    ebek.push_back(-pow(kp1,2)*cos(kp0*xl[0])*sin(kp1*xl[1]));  /** eyy **/
    ebek.push_back(-2*kp0*kp1*sin(kp0*xl[0])*cos(kp1*xl[1]));   /** exy **/
    ebek.push_back(-kp0*sin(kp0*xl[0])*sin(kp1*xl[1]));         /** exz **/
    ebek.push_back(kp1*cos(kp0*xl[0])*cos(kp1*xl[1]));          /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** P4 = cos(kp0*xl[0])*cos(kp1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-pow(kp0,2)*cos(kp0*xl[0])*cos(kp1*xl[1]));  /** exx **/
    ebek.push_back(-pow(kp1,2)*cos(kp0*xl[0])*cos(kp1*xl[1]));  /** eyy **/
    ebek.push_back(2*kp0*kp1*sin(kp0*xl[0])*sin(kp1*xl[1]));    /** exy **/
    ebek.push_back(-kp0*sin(kp0*xl[0])*cos(kp1*xl[1]));         /** exz **/
    ebek.push_back(-kp1*cos(kp0*xl[0])*cos(kp1*xl[1]));         /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** S1 = sin(ks0*xl[0])*sin(ks1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-zo*pow(ks0,2)*sin(ks0*xl[0])*sin(ks1*xl[1]));  /** exx **/
    ebek.push_back(-zo*pow(ks1,2)*sin(ks0*xl[0])*sin(ks1*xl[1]));  /** eyy **/
    ebek.push_back(2*zo*ks0*ks1*cos(ks0*xl[1])*cos(ks1*xl[1]));    /** exy **/
    ebek.push_back(ks0*cos(ks0*xl[0])*sin(ks1*xl[1]));             /** exz **/
    ebek.push_back(ks1*sin(ks0*xl[0])*cos(ks1*xl[1]));             /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** S2 = sin(ks0*xl[0])*cos(ks1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-zo*pow(ks0,2)*sin(ks0*xl[0])*cos(ks1*xl[1]));  /** exx **/
    ebek.push_back(-zo*pow(ks1,2)*sin(ks0*xl[0])*cos(ks1*xl[1]));  /** eyy **/
    ebek.push_back(-2*zo*ks0*ks1*cos(ks0*xl[1])*sin(ks1*xl[1]));   /** exy **/
    ebek.push_back(ks0*cos(ks0*xl[0])*cos(ks1*xl[1]));             /** exz **/
    ebek.push_back(-ks1*sin(ks0*xl[0])*sin(ks1*xl[1]));            /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** S3 = cos(ks0*xl[0])*sin(ks1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-zo*pow(ks0,2)*cos(ks0*xl[0])*sin(ks1*xl[1]));  /** exx **/
    ebek.push_back(-zo*pow(ks1,2)*cos(ks0*xl[0])*sin(ks1*xl[1]));  /** eyy **/
    ebek.push_back(-2*zo*ks0*ks1*sin(ks0*xl[1])*cos(ks1*xl[1]));   /** exy **/
    ebek.push_back(-ks0*sin(ks0*xl[0])*sin(ks1*xl[1]));            /** exz **/
    ebek.push_back(ks1*cos(ks0*xl[0])*cos(ks1*xl[1]));             /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    /** S4 = cos(ks0*xl[0])*cos(ks1*xl[1]) **/

    /** oskol-sistemako deformazioak **/
    ebek.clear();
    ebek.push_back(-zo*pow(ks0,2)*cos(ks0*xl[0])*cos(ks1*xl[1]));  /** exx **/
    ebek.push_back(-zo*pow(ks1,2)*cos(ks0*xl[0])*cos(ks1*xl[1]));  /** eyy **/
    ebek.push_back(2*zo*ks0*ks1*sin(ks0*xl[1])*sin(ks1*xl[1]));    /** exy **/
    ebek.push_back(-ks0*sin(ks0*xl[0])*cos(ks1*xl[1]));            /** exz **/
    ebek.push_back(-ks1*cos(ks0*xl[0])*cos(ks1*xl[1]));            /** exz **/

    /** oskol-sistemako tentsioak **/
    tbek.clear();
    tbek.push_back(dmat[0][0]*ebek[0]+dmat[0][1]*ebek[1]+dmat[0][2]*ebek[2]+dmat[0][3]*ebek[3]+dmat[0][4]*ebek[4]);  /** Sxx **/
    tbek.push_back(dmat[1][0]*ebek[0]+dmat[1][1]*ebek[1]+dmat[1][2]*ebek[2]+dmat[1][3]*ebek[3]+dmat[1][4]*ebek[4]);  /** Syy **/
    tbek.push_back(dmat[2][0]*ebek[0]+dmat[2][1]*ebek[1]+dmat[2][2]*ebek[2]+dmat[2][3]*ebek[3]+dmat[2][4]*ebek[4]);  /** Sxy **/
    tbek.push_back(dmat[3][0]*ebek[0]+dmat[3][1]*ebek[1]+dmat[3][2]*ebek[2]+dmat[3][3]*ebek[3]+dmat[3][4]*ebek[4]);  /** Sxz **/
    tbek.push_back(dmat[4][0]*ebek[0]+dmat[4][1]*ebek[1]+dmat[4][2]*ebek[2]+dmat[4][3]*ebek[3]+dmat[4][4]*ebek[4]);  /** Syz **/

    /** Sistema orokorreko deformazioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*ebek[1]+b1[0]*(b1[0]*ebek[0]+2*b2[0]*ebek[2]+2*normal[0]*ebek[3])+2*b2[0]*normal[0]*ebek[4]);                                                 /** exx **/
    bek.push_back(b2[1]*b2[1]*ebek[1]+b1[1]*(b1[1]*ebek[0]+2*b2[1]*ebek[2]+2*normal[1]*ebek[3])+2*b2[1]*normal[1]*ebek[4]);                                                 /** eyy **/
    bek.push_back(b2[2]*b2[2]*ebek[1]+b1[2]*(b1[2]*ebek[0]+2*b2[2]*ebek[2]+2*normal[2]*ebek[3])+2*b2[2]*normal[2]*ebek[4]);                                                 /** ezz **/
    bek.push_back(b1[0]*(b1[1]*ebek[0]+b2[1]*ebek[2]+normal[1]*ebek[3])+normal[0]*(b1[1]*ebek[3]+b2[1]*ebek[4])+b2[0]*(b2[1]*ebek[1]+b1[1]*ebek[2]+normal[1]*ebek[4]));     /** exy **/
    bek.push_back(b1[0]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[0]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[0]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** exz **/
    bek.push_back(b1[1]*(b1[2]*ebek[0]+b2[2]*ebek[2]+normal[2]*ebek[3])+normal[1]*(b1[2]*ebek[3]+b2[2]*ebek[4])+b2[1]*(b2[2]*ebek[1]+b1[2]*ebek[2]+normal[2]*ebek[4]));     /** eyz **/
    dserie.push_back(bek);

    /** Sistema orokorreko tentsioak **/
    bek.clear();
    bek.push_back(b2[0]*b2[0]*tbek[1]+b1[0]*(b1[0]*tbek[0]+2*b2[0]*tbek[2]+2*normal[0]*tbek[3])+2*b2[0]*normal[0]*tbek[4]);                                                 /** Sxx **/
    bek.push_back(b2[1]*b2[1]*tbek[1]+b1[1]*(b1[1]*tbek[0]+2*b2[1]*tbek[2]+2*normal[1]*tbek[3])+2*b2[1]*normal[1]*tbek[4]);                                                 /** Syy **/
    bek.push_back(b2[2]*b2[2]*tbek[1]+b1[2]*(b1[2]*tbek[0]+2*b2[2]*tbek[2]+2*normal[2]*tbek[3])+2*b2[2]*normal[2]*tbek[4]);                                                 /** Szz **/
    bek.push_back(b1[0]*(b1[1]*tbek[0]+b2[1]*tbek[2]+normal[1]*tbek[3])+normal[0]*(b1[1]*tbek[3]+b2[1]*tbek[4])+b2[0]*(b2[1]*tbek[1]+b1[1]*tbek[2]+normal[1]*tbek[4]));     /** Sxy **/
    bek.push_back(b1[0]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[0]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[0]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Sxz **/
    bek.push_back(b1[1]*(b1[2]*tbek[0]+b2[2]*tbek[2]+normal[2]*tbek[3])+normal[1]*(b1[2]*tbek[3]+b2[2]*tbek[4])+b2[1]*(b2[2]*tbek[1]+b1[2]*tbek[2]+normal[2]*tbek[4]));     /** Syz **/
    tserie.push_back(bek);

    } /** end of for m=1;nhmoz **/

    theta+=tartea;

}  /** end of function oskolserie **/

