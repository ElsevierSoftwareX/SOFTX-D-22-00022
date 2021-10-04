#include <iomanip>
#include <sstream>
#include <math.h>
#include <complex>
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
#include "../erkideak/ObjektuPuntu.h"
#include "../erkideak/ObjektuAdi.h"
#include "../scattering/neumann.h"
const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

using namespace std;

int main()
{

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
string adi, pun, flaviares;

unsigned int maizkont, kont, puntua;

double f, k;
vector<complex<double> > gradientea; /** Nodoko Gradientea sartzeko **/
complex<double> hbal;
/**
ama ObjektuPost objektua sortu adi fitxategitik
**/

cout << "Objektuaren .adi fitxategiaren izena? " << endl;
cin >> adi;

ObjektuAdi ama(adi); /** Objektu sortzailea **/

unsigned int maizkop=ama.maizkop;
unsigned int nhmoz;

/**
kume ObjektuAurre objektua sortu mub fitxategitik
**/

cout << "Bolumeneko puntu hodeia daukan .pun fitxategiaren izena? " << endl;
cin >> pun;
//ObjektuPuntu kume(pun);
ObjektuPuntu kume(pun, 0, 1, -1);

/**
flavia.res fitxategia idatzi
**/

flaviares = bidea + pun + ".flavia.res";

vector<vector<complex<double> > > kserie; /** Serie konplexua **/
vector<vector<complex<double> > > dkserie; /** Serie konplexuaren deribatuak cartesiarretan **/

if (ama.mota.compare("ELAS")==0)

/************************* EGITURA ELASTIKOA ********************/
{ /** Egitura elastiko objektuaren prozesamendua **/


} /** Egitura Elastiko (ELAS) objektuaren amaiera **/

else if (ama.mota.compare("FFIN")==0)
{ /** Fluido Finitu (FFIN) objektuaren prozesamendua**/

    ofstream resfitx(flaviares.c_str());

    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        nhmoz=ama.shkop[maizkont];
        f=ama.espektrua[maizkont];
        k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/
        resfitx << setfill('0');
        resfitx << "dp/dn" << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        resfitx << fixed << setprecision(2) << "Modu  f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Real  f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Imag  f = " << f << "Hz." << endl;

//        resfitx << "Gradient" << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
//        resfitx << fixed << setprecision(2) << "Gradx f = " << f << "Hz." << endl;
//        resfitx << fixed << setprecision(2) << "Grady f = " << f << "Hz." << endl;
//        resfitx << fixed << setprecision(2) << "Gradz f = " << f << "Hz." << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            kserie.clear();
            dkserie.clear();
            ffinserie(ama.zentroa, kume.puntuak[puntua], k, nhmoz, kserie, dkserie);

            gradientea.clear();
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);

            hbal=0.0+0.0*I;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                // cout << "KOEF[" << maizkont << "," << kont << "]= " << ama.koef[maizkont][kont] << endl;

                gradientea[0]=gradientea[0]+ama.koef[maizkont][kont]*dkserie[kont][0];
                gradientea[1]=gradientea[1]+ama.koef[maizkont][kont]*dkserie[kont][1];
                gradientea[2]=gradientea[2]+ama.koef[maizkont][kont]*dkserie[kont][2];
                hbal=hbal+ama.koef[maizkont][kont]*(dkserie[kont][0]*kume.normalak[puntua][0]+dkserie[kont][1]*kume.normalak[puntua][1]+dkserie[kont][2]*kume.normalak[puntua][2]);
            }

//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << real(gradientea[0]-hbal*normal[0]) << "\t" << real(gradientea[1]-hbal*normal[1]) << "\t" << real(gradientea[2]-hbal*normal[2]) << endl;
//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << -real(gradientea[0]) << "\t" << -real(gradientea[1]) << "\t" << -real(gradientea[2]) << endl;
            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(hbal) << "\t" << real(hbal) << "\t" << imag(hbal) << endl;
//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << -normal[0] << "\t" << -normal[1] << "\t" << -normal[2] << endl;
        }

    }
    resfitx.close();

} /** Fluido Finitu (FFIN) objektuaren amaiera**/

else if (ama.mota.compare("FINF")==0)
{ /** Fluido Infinitu (FINF) objektuaren prozesamendua**/

    ofstream resfitx(flaviares.c_str());

    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        nhmoz=ama.shkop[maizkont];
        f=ama.espektrua[maizkont];
        k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/
        resfitx << setfill('0');

        resfitx << "dp/dn" << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        resfitx << fixed << setprecision(2) << "Modu  f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Real  f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Imag  f = " << f << "Hz." << endl;

//        resfitx << "Gradient" << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
//        resfitx << fixed << setprecision(2) << "Gradx f = " << f << "Hz." << endl;
//        resfitx << fixed << setprecision(2) << "Grady f = " << f << "Hz." << endl;
//        resfitx << fixed << setprecision(2) << "Gradz f = " << f << "Hz." << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            //cout << "n(" << kume.normalak[puntua][0] << "," << kume.normalak[puntua][1] << "," << kume.normalak[puntua][2] << ") " << endl;
            kserie.clear();
            dkserie.clear();
            finfserie(ama.zentroa, kume.puntuak[puntua], k, nhmoz, kserie, dkserie);

            gradientea.clear();
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);

            hbal=0.0+0.0*I;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                // cout << "KOEF[" << maizkont << "," << kont << "]= " << ama.koef[maizkont][kont] << endl;

                gradientea[0]=gradientea[0]+ama.koef[maizkont][kont]*dkserie[kont][0];
                gradientea[1]=gradientea[1]+ama.koef[maizkont][kont]*dkserie[kont][1];
                gradientea[2]=gradientea[2]+ama.koef[maizkont][kont]*dkserie[kont][2];
                hbal=hbal+ama.koef[maizkont][kont]*(dkserie[kont][0]*kume.normalak[puntua][0]+dkserie[kont][1]*kume.normalak[puntua][1]+dkserie[kont][2]*kume.normalak[puntua][2]);
            }

//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << real(gradientea[0]-hbal*normal[0]) << "\t" << real(gradientea[1]-hbal*normal[1]) << "\t" << real(gradientea[2]-hbal*normal[2]) << endl;
//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << -real(gradientea[0]) << "\t" << -real(gradientea[1]) << "\t" << -real(gradientea[2]) << endl;
            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(hbal) << "\t" << real(hbal) << "\t" << imag(hbal) << endl;
//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << -normal[0] << "\t" << -normal[1] << "\t" << -normal[2] << endl;
        }

    }
    resfitx.close();
} /** Fluido Infinitu (FINF) objektuaren amaiera**/

return 0;
}
