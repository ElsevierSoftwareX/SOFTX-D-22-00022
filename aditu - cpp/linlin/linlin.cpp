#include <iomanip>
#include <sstream>
#include <math.h>
#include <complex>
#include "../erkideak/elasserie.h"
#include "wsquare.h"
#include "../erkideak/ObjektuPuntu.h"
#include "../erkideak/ObjektuAdi.h"
#include "../erkideak/postmsh.h"
const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

using namespace std;

int main()
{

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
string adi, pun, postres;

unsigned int maizkont, kont, puntua;

double cp, kp, cs, ks, f, omega, young, poisson, lambda, mu, rho, d, h=0.002, a=1.0, b=1.0, wsq;
complex<double> gezia; /** Nodoko gezia sartzeko **/
complex<double> pbal;

unsigned int nhmoz;

/**
ama ObjektuPost objektua sortu adi fitxategitik
**/

cout << "Objektuaren .adi fitxategiaren izena? " << endl;
cin >> adi;

ObjektuAdi ama(adi); /** Objektu sortzailea **/

unsigned int maizkop=ama.maizkop;

/**
kume ObjektuAurre objektua sortu mub fitxategitik
**/

cout << "Bolumeneko puntu hodeia daukan .pun fitxategiaren izena? " << endl;
cin >> pun;
//ObjektuPuntu kume(pun, 0, 1, -1, 0); /** Neuk sortutako .pun fitxategia **/
ObjektuPuntu kume(pun, 0, 1, -1); /** GiDek sortutako .pun fitxategia **/

postmsh(pun);

/**
post.res fitxategia idatzi
**/

postres = bidea + pun + ".post.res";

vector<vector<complex<double> > > serie; /** Serie konplexua **/
vector<vector<complex<double> > > dserie; /** Serie konplexuaren deribatuak cartesiarretan **/
vector<vector<complex<double> > > tserie; /** Serie konplexuaren deribatuak cartesiarretan (tentsioak) **/

if (ama.mota.compare("ELAS")==0)

/************************* EGITURA ELASTIKOA ********************/
{ /** Egitura elastiko objektuaren prozesamendua **/

ofstream resfitx(postres.c_str());
resfitx << "GiD Post Results File 1.0" << endl;

    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        nhmoz=ama.shkop[maizkont];
        f=ama.espektrua[maizkont];
        omega=2*PI*f;

        young=ama.mat[0]; /** Young modulua **/
        poisson=ama.mat[1]; /** Poisson ratioa **/
        rho=ama.mat[2]; /** Dentsitatea **/
        d=young*pow(h,3)/(12*(1-pow(poisson,2)));

        lambda=(young*poisson)/((1-2*poisson)*(1+poisson)); /** Lameren 1. parametroa **/
        mu=young/(1+poisson);                               /** Lameren 2. parametroa **/
        cp=sqrt((lambda+2*poisson)/(rho)); /** 1. uhin-abiadura **/
        kp=2*PI*f/cp; /** P uhinaren uhin-zenbakia k1=2*pi*f/cp **/
        cs=sqrt(mu/rho); /** 2. uhin-abiadura **/
        ks=2*PI*f/cs; /** S uhinaren uhin-zenbakia k1=2*pi*f/cs **/

        resfitx << "ResultGroup" << " " << (char)34 << pun << (char)34 << " " << maizkont+1 << " OnNodes" << endl;
        resfitx << "ResultDescription" << " " << (char)34 << "Nodal Displacements" << (char)34 << " " << "Vector" << endl;
        resfitx << "ComponentNames ";
        resfitx << (char)34 << "uz zehatza   f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "uz hurbildua f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "Errorea      f= " << setprecision(2) << f << " Hz." << (char)34 << endl;

        resfitx << "Values" << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            wsquare(kume.puntuak[puntua], omega, d, h, rho, a, b, wsq);

            serie.clear();
            dserie.clear();
            tserie.clear();
            elasserie(ama.zentroa, kume.puntuak[puntua], kp, ks, young, poisson, nhmoz, serie, dserie, tserie);

            gezia=0.0;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                gezia=gezia+abs(ama.koef[maizkont][kont]*serie[kont][2]);
            }

            resfitx << setprecision(10) << scientific << puntua+1 << "\t" << abs(wsq) << "\t" << abs(gezia) << "\t" << abs(wsq-gezia) << endl;
        }
        resfitx << "end values" << endl;
    }
    resfitx.close();

} /** Egitura Elastiko (ELAS) objektuaren amaiera **/

return 0;
}
