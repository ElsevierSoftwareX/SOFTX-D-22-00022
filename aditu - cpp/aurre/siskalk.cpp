#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include "siskalk.h"
#include "../erkideak/FourTrans.h"
#include "../erkideak/elasserie.h"
#include "../erkideak/oskoserie.h"
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
#include "../scattering/neumann.h"

const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

using namespace std;

void siskalk(const vector<ObjektuPuntu> &puntu, const unsigned int &objkop, const vector<vector<unsigned int> > &ekkop, const int harmonikoa, const unsigned int maizkont, const string taula, vector<vector<complex<double> > > &amat, vector<complex<double> > &bbek, vector<vector <unsigned int> > &zutabeak, unsigned int &prokoefkop)
{
FourTrans ft(taula);

vector<complex<double> > lerroa; /** amat matrizeari eransteko errenkada **/

vector<unsigned int> itemp; /** Laguntzaileak **/

unsigned int objin, objout, objkont;
unsigned int puntua, puntukop;
unsigned int kont=0, kont2=0;

vector<vector<complex<double> > >  seriein;  /** Serie konplexua **/
vector<vector<complex<double> > > dseriein;  /** Serie konplexuaren deribatuak cartesiarretan **/
vector<vector<complex<double> > > tseriein;  /** Serie konplexuaren tentsioak cartesiarretan **/

vector<vector<complex<double> > >  serieout; /** Serie konplexua **/
vector<vector<complex<double> > > dserieout; /** Serie konplexuaren deribatuak cartesiarretan **/
vector<vector<complex<double> > > tserieout; /** Serie konplexuaren deribatuak cartesiarretan **/

complex<double> cbal;
vector<complex<double> > cbek; /** bektore laguntzailea **/

/** Egitura elastikoen propietateak **/

double lodiera, youngin, youngout, poissonin, poissonout, rho, lambda, mu; /** Young moduluak, Poisson ratioak, Dentsitatea... **/
double cp, cs; /** Uhin-abiadurak **/
double kpin, kpout, ksin, ksout; /** P eta S uhinen uhin-zenbakiak k=2*pi*f/c **/
double zoin, zoout; /** Oskol sistemako z koordenatua (t/2 edo -t/2)**/

/** Fluidoaren propietateak**/

double kin, kout;

amat.clear();
bbek.clear();
zutabeak.clear();
prokoefkop=0;

for (objkont=0; objkont<objkop; objkont++)
{
    itemp.clear();
    itemp.push_back(kont2);

    if (puntu[objkont].mota=="ELAS")
    {
        if (harmonikoa < puntu[objkont].harmonikomax)
        {
            prokoefkop=prokoefkop+8*(harmonikoa+1);
            kont2=kont2+8*(harmonikoa+1);
        }
        else
        {
            prokoefkop=prokoefkop+8*(puntu[objkont].harmonikomax+1);
            kont2=kont2+8*(puntu[objkont].harmonikomax+1);
        }
    } /** ELAS **/

    if (puntu[objkont].mota=="OSKO")
    {
        if (harmonikoa < puntu[objkont].harmonikomax)
        {
            prokoefkop=prokoefkop+8*(harmonikoa+1);
            kont2=kont2+8*(harmonikoa+1);
        }
        else
        {
            prokoefkop=prokoefkop+8*(puntu[objkont].harmonikomax+1);
            kont2=kont2+8*(puntu[objkont].harmonikomax+1);
        }
    } /** OSKO **/

    if (puntu[objkont].mota=="FFIN")
    {
        if (harmonikoa < puntu[objkont].harmonikomax)
        {
            prokoefkop=prokoefkop+2*pow(harmonikoa+1,2);
            kont2=kont2+2*pow(harmonikoa+1,2);
        }
        else
        {
            prokoefkop=prokoefkop+2*pow(puntu[objkont].harmonikomax+1,2);
            kont2=kont2+2*pow(puntu[objkont].harmonikomax+1,2);
        }
    } /** FFIN **/

    if (puntu[objkont].mota=="FINF")
    {
        if (harmonikoa < puntu[objkont].harmonikomax)
        {
            prokoefkop=prokoefkop+pow(harmonikoa+1,2);
            kont2=kont2+pow(harmonikoa+1,2);
        }
        else
        {
            prokoefkop=prokoefkop+pow(puntu[objkont].harmonikomax+1,2);
            kont2=kont2+pow(puntu[objkont].harmonikomax+1,2);
        }
    } /** FINF **/

    itemp.push_back(kont2-1);
    zutabeak.push_back(itemp);

} /** end of for obj **/


vector<vector<unsigned int> > ekkont;
ekkont.clear();
ekkont.resize(objkop, vector<unsigned int>(objkop,0));

/************************* OBJEKTUEN PROZESAMENDUA ********************/

for (objin=0; objin<objkop; objin++)
{
    /************************* ALDAGAI KOMUNAK OBJEKTURAKO ********************/

    if (puntu[objin].mota.compare("ELAS")==0)

    /************************* EGITURA ELASTIKOA ********************/

    { /** Egitura elastiko objektuaren prozesamendua **/

    puntukop=puntu[objin].puntukop;

    youngin=puntu[objin].mat[0];      /** Young modulua **/
    poissonin=puntu[objin].mat[1];    /** Poisson ratioa **/
    rho=puntu[objin].mat[2];        /** Dentsitatea **/

    lambda=(youngin*poissonin)/((1-2*poissonin)*(1+poissonin)); /** Lameren 1. parametroa **/
    mu=youngin/(2*(1+poissonin));                           /** Lameren 2. parametroa **/
    cp=sqrt((lambda+2*poissonin)/rho);                    /** 1. uhin-abiadura **/
    kpin=2*PI*ft.espektrua[maizkont]/cp;                  /** P uhinaren uhin-zenbakia k1=2*pi*f/cp **/
    cs=sqrt(mu/rho);                                    /** 2. uhin-abiadura **/
    ksin=2*PI*ft.espektrua[maizkont]/cs;                  /** S uhinaren uhin-zenbakia k1=2*pi*f/cs **/

    for (puntua=0; puntua<puntukop; puntua++)
    {
        seriein.clear();
        dseriein.clear();
        tseriein.clear();
        elasserie(puntu[objin].zentroa, puntu[objin].puntuak[puntua], kpin, ksin, lambda, mu, harmonikoa, seriein, dseriein, tseriein);

        switch (puntu[objin].mbzenb[puntua])
        {
        case 1: /** Bermatua-lerroa **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of while ekkont **/
            break;

        case 2: /** Bermatua-gainazala **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of while ekkont **/
            break;

        case 3: /** Landatua **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of while ekkont **/
            break;

        case 4: /** Labainkorra **/

            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0]*puntu[objin].normalak[puntua][0]+seriein[kont2][1]*puntu[objin].normalak[puntua][1]+seriein[kont2][2]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 5: /** Tentsio-bektorea **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** T_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    /** Tx = sxx*nx + sxy*ny* sxz*nz **/
                    lerroa[kont]=(1/lambda)*(tseriein[kont2][0]*puntu[objin].normalak[puntua][0]+tseriein[kont2][1]*puntu[objin].normalak[puntua][1]+tseriein[kont2][2]*puntu[objin].normalak[puntua][2]);
                    //lerroa[kont]=tseriein[kont2][0]*puntu[objin].normalak[puntua][0]+tseriein[kont2][1]*puntu[objin].normalak[puntua][1]+tseriein[kont2][2]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back((1/lambda)*(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]));  /** Taulako balioa **/
                    //bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

            /** T_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    /** Ty = syx*nx + syy*ny* syz*nz **/
                    lerroa[kont]=(1/lambda)*(tseriein[kont2][3]*puntu[objin].normalak[puntua][0]+tseriein[kont2][4]*puntu[objin].normalak[puntua][1]+tseriein[kont2][5]*puntu[objin].normalak[puntua][2]);
                    //lerroa[kont]=tseriein[kont2][3]*puntu[objin].normalak[puntua][0]+tseriein[kont2][4]*puntu[objin].normalak[puntua][1]+tseriein[kont2][5]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back((1/lambda)*(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]));  /** Taulako balioa **/
                    //bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** T_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    /** Tz = szx*nx + szy*ny* szz*nz **/
                    //lerroa[kont]=(1/lambda)*(tseriein[kont2][6]*puntu[objin].normalak[puntua][0]+tseriein[kont2][7]*puntu[objin].normalak[puntua][1]+tseriein[kont2][8]*puntu[objin].normalak[puntua][2]);
                    lerroa[kont]=tseriein[kont2][6]*puntu[objin].normalak[puntua][0]+tseriein[kont2][7]*puntu[objin].normalak[puntua][1]+tseriein[kont2][8]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    //bbek.push_back((1/lambda)*(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]));  /** Taulako balioa **/
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 6: /** Egitura-akoplamendua-Dirichlet **/

            cout << "ELAS:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 7: /** Egitura-akoplamendua-labainkorra **/

            cout << "ELAS:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 8: /** Fluido-finitu-akoplamendua **/

            cout << "ELAS:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 9: /** Fluido-infinitu-akoplamendua **/

            cout << "ELAS:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        } /** switch mbmota switcharen amaiera **/

    } /** Egitura Elastiko (ELAS) objektuaren Mugalde Baldintzen prozesamenduaren amaiera **/

    } /** Egitura Elastiko (ELAS) objektuaren amaiera **/

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    else if (puntu[objin].mota.compare("OSKO")==0)

    /************************* OSKOL ELASTIKOA ********************/

    { /** Oskol elastiko objektuaren prozesamendua **/

    puntukop=puntu[objin].puntukop;

    lodiera=puntu[objin].mat[0];    /** t lodiera **/
    youngin=puntu[objin].mat[1];    /** Young modulua **/
    poissonin=puntu[objin].mat[2];  /** Poisson ratioa **/
    rho=puntu[objin].mat[3];      /** Dentsitatea **/

    lambda=(youngin*poissonin)/((1+poissonin)*(1-2*poissonin));
    mu=youngin/(2*(1+poissonin));

    cp=sqrt((lambda+2*mu)/rho);                           /** 1. uhin-abiadura **/
    kpin=2*PI*ft.espektrua[maizkont]/cp;                  /** P uhinaren uhin-zenbakia k1=2*pi*f/cp **/
    cs=sqrt(mu/rho);                                      /** 2. uhin-abiadura **/
    ksin=2*PI*ft.espektrua[maizkont]/cs;                  /** S uhinaren uhin-zenbakia k1=2*pi*f/cs **/

    for (puntua=0; puntua<puntukop; puntua++)
    {
        seriein.clear();
        dseriein.clear();
        tseriein.clear();
        /** Oskol sistemako z koordenatua (t/2 edo -t/2)**/
        zoin=0.5*lodiera;
        oskoserie(puntu[objin].zentroa, puntu[objin].puntuak[puntua], puntu[objin].normalak[puntua], zoin, kpin, ksin, youngin, poissonin, harmonikoa, seriein, dseriein, tseriein);

        cout << setprecision(20) << "serie= " << seriein[0][0] << " , " << seriein[1][0] << " , " << seriein[2][0] << " , " << seriein[3][0] << endl;


        switch (puntu[objin].mbzenb[puntua])
        {
        case 1: /** Bermatua **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of while ekkont **/
            break;


        case 2: /** Landatua **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** u_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of while ekkont **/
            break;

        case 3: /** Labainkorra **/

            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0]*puntu[objin].normalak[puntua][0]+seriein[kont2][1]*puntu[objin].normalak[puntua][1]+seriein[kont2][2]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 4: /** Tentsio-bektorea **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** T_x **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    /** Tx = sxx*nx + sxy*ny* sxz*nz **/
                    lerroa[kont]=(1/lambda)*(tseriein[kont2][0]*puntu[objin].normalak[puntua][0]+tseriein[kont2][1]*puntu[objin].normalak[puntua][1]+tseriein[kont2][2]*puntu[objin].normalak[puntua][2]);
                    //lerroa[kont]=tseriein[kont2][0]*puntu[objin].normalak[puntua][0]+tseriein[kont2][1]*puntu[objin].normalak[puntua][1]+tseriein[kont2][2]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

            /** T_y **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    /** Ty = syx*nx + syy*ny* syz*nz **/
                    //lerroa[kont]=tseriein[kont2][3]*puntu[objin].normalak[puntua][0]+tseriein[kont2][4]*puntu[objin].normalak[puntua][1]+tseriein[kont2][5]*puntu[objin].normalak[puntua][2];
                    lerroa[kont]=(1/lambda)*(tseriein[kont2][3]*puntu[objin].normalak[puntua][0]+tseriein[kont2][4]*puntu[objin].normalak[puntua][1]+tseriein[kont2][5]*puntu[objin].normalak[puntua][2]);
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][1]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][1]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                /** T_z **/

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    /** Tz = szx*nx + szy*ny* szz*nz **/
                    //lerroa[kont]=tseriein[kont2][6]*puntu[objin].normalak[puntua][0]+tseriein[kont2][7]*puntu[objin].normalak[puntua][1]+tseriein[kont2][8]*puntu[objin].normalak[puntua][2];
                    lerroa[kont]=(1/lambda)*(tseriein[kont2][6]*puntu[objin].normalak[puntua][0]+tseriein[kont2][7]*puntu[objin].normalak[puntua][1]+tseriein[kont2][8]*puntu[objin].normalak[puntua][2]);
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][2]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    //bbek.push_back((1/lambda)*(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]));  /** Taulako balioa **/
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][2]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 5: /** Egitura-akoplamendua-Dirichlet **/

            cout << "OSKO:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 6: /** Egitura-akoplamendua-labainkorra **/

            cout << "OSKO:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 7: /** Fluido-finitu-akoplamendua **/

            cout << "OSKO:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 8: /** Fluido-infinitu-akoplamendua **/

            cout << "OSKO:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        } /** switch mbmota switcharen amaiera **/

    } /** Egitura Elastiko (ELAS) objektuaren Mugalde Baldintzen prozesamenduaren amaiera **/

    } /** Egitura Elastiko (ELAS) objektuaren amaiera **/

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    else if (puntu[objin].mota.compare("FFIN")==0)

    /************************* FLUIDO FINITUA ********************/

    { /** Fluido Finitu (FFIN) objektuaren   prozesamendua**/

    puntukop=puntu[objin].puntukop;
    kin=2*PI*ft.espektrua[maizkont]/puntu[objin].mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

    for (puntua=0; puntua<puntukop; puntua++)
    {
        seriein.clear();
        dseriein.clear();
        ffinserie(puntu[objin].zentroa, puntu[objin].puntuak[puntua], kin, harmonikoa, seriein, dseriein);

        switch (puntu[objin].mbzenb[puntua])
        {
        case 1: /** Dirichlet-ezarria **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of while ekkont **/
            break;

        case 2: /** Neumann-ezarria **/

            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][0]*puntu[objin].normalak[puntua][0]+dseriein[kont2][1]*puntu[objin].normalak[puntua][1]+dseriein[kont2][2]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    //neumann(puntu[objin].puntuak[puntua], kin, 15, cbal);
                    //bbek.push_back(cbal);  /** Scattering **/
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 3: /** Egitura-akoplamendua **/
            cout << "FINF:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 4: /** Fluido-finitu-akoplamendua **/

            objout=puntu[objin].mbobj[puntua][0];
            kout=2*ft.espektrua[maizkont]/puntu[objout].mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

            if (ekkont[objin][objout] < ekkop[objin][objout])
            {
                serieout.clear();
                dserieout.clear();
                ffinserie(puntu[objout].zentroa, puntu[objin].puntuak[puntua], kout, harmonikoa, serieout, dserieout);

                /** Continuity of acoustic pressure **/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-serieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_x**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_y**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][1];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_z**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][2];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

            } /** end of if ekkont **/
            break;

        case 5: /** Fluido-infinitu-akoplamendua **/

            objout=puntu[objin].mbobj[puntua][0];
            kout=2*PI*ft.espektrua[maizkont]/puntu[objout].mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

            if (ekkont[objin][objout] < ekkop[objin][objout])
            {
                serieout.clear();
                dserieout.clear();
                finfserie(puntu[objout].zentroa, puntu[objin].puntuak[puntua], kout, harmonikoa, serieout, dserieout);

                /** Continuity of acoustic pressure **/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-serieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_x**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_y**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][1];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_z**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][2];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

            } /** end of if ekkont **/

            break;

            } /** switch mbmota switcharen amaiera **/

        } /** Fluido Infinitu (FINF) objektuaren Mugalde Baldintzen prozesamenduaren amaiera **/

    } /** Fluido Finitu (FFIN) objektuaren amaiera **/

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    else if (puntu[objin].mota.compare("FINF")==0)

    /************************* FLUIDO INFINITUA ********************/

    { /** Fluido Infinitu (FINF) objektuaren prozesamendua**/

    puntukop=puntu[objin].puntukop;
    kin=2*PI*ft.espektrua[maizkont]/puntu[objin].mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

    for (puntua=0; puntua<puntukop; puntua++)
        {
        seriein.clear();
        dseriein.clear();
        finfserie(puntu[objin].zentroa, puntu[objin].puntuak[puntua], kin, harmonikoa, seriein, dseriein);

        switch (puntu[objin].mbzenb[puntua])
        {
        case 1: /** Dirichlet-ezarria **/
            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 2: /** Neumann-ezarria **/

            if (ekkont[objin][objin] < ekkop[objin][objin])
            {
                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][0]*puntu[objin].normalak[puntua][0]+dseriein[kont2][1]*puntu[objin].normalak[puntua][1]+dseriein[kont2][2]*puntu[objin].normalak[puntua][2];
                    kont2++;
                }
                amat.push_back(lerroa);

                if (puntu[objin].mbobj[puntua][0]==-1)
                {
                    bbek.push_back(0.0+0.0*I);  /** MB homogeneoa **/
                }
                else
                {
                    //neumann(puntu[objin].puntuak[puntua], kin, 15, cbal);
                    //bbek.push_back(cbal);  /** Scattering **/
                    bbek.push_back(ft.taula[maizkont][puntu[objin].mbobj[puntua][0]]);  /** Taulako balioa **/
                }

                ekkont[objin][objin]++;

            } /** end of if ekkont **/
            break;

        case 3: /** Egitura-akoplamendua **/
            cout << "FINF:  "<< puntua << " baldintza Egitura akoplamendua motakoa da" << endl;
            break;

        case 4: /** Fluido-finitu-akoplamendua **/

            objout=puntu[objin].mbobj[puntua][0];
            kout=2*PI*ft.espektrua[maizkont]/puntu[objout].mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

            if (ekkont[objin][objout] < ekkop[objin][objout])
            {

                serieout.clear();
                dserieout.clear();
                ffinserie(puntu[objout].zentroa, puntu[objin].puntuak[puntua], kout, harmonikoa, serieout, dserieout);

                /** Continuity of acoustic pressure **/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-serieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_x**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_y**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][1];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_z**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][2];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

            } /** end of if ekkont **/
            break;

        case 5: /** Fluido-infinitu-akoplamendua **/

            objout=puntu[objin].mbobj[puntua][0];
            kout=2*PI*ft.espektrua[maizkont]/puntu[objout].mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

            if (ekkont[objin][objout] < ekkop[objin][objout])
            {

                serieout.clear();
                dserieout.clear();
                finfserie(puntu[objout].zentroa, puntu[objin].puntuak[puntua], kout, harmonikoa, serieout, dserieout);

                /** Continuity of acoustic pressure **/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=seriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-serieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_x**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][0];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][0];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_y**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][1];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][1];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

                /** Continuity of the gradient of acoustic pressure **/
                /** Grad_z**/

                lerroa.clear();
                lerroa.resize(prokoefkop, 0.0+0.0*I);

                kont2=0;
                for (kont = zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
                {
                    lerroa[kont]=dseriein[kont2][2];
                    kont2++;
                }

                kont2=0;
                for (kont = zutabeak[objout][0]; kont <= zutabeak[objout][1]; kont++)
                {
                    lerroa[kont]=-dserieout[kont2][2];
                    kont2++;
                }

                amat.push_back(lerroa);
                bbek.push_back(0.0+0.0*I);

                ekkont[objin][objout]++;

            } /** end of if ekkont **/
            break;

            } /** switch mbmota switcharen amaiera **/

        } /** Fluido Infinitu (FINF) objektuaren Mugalde Baldintzen prozesamenduaren amaiera **/

    } /** Fluido Infinitu (FINF) objektuaren prozesamenduaren amaiera **/

} /** objin < objkop **/

}
