#include <iomanip>
#include <sstream>
#include <math.h>
#include <complex>
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
#include "../erkideak/ObjektuPuntu.h"
#include "../erkideak/ObjektuAdi.h"
#include "scattered.h"
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
complex<double> balio;
complex<double> baliozehatz, baliohur;
double baliozehatzdb, baliohurdb, argzehatz, arghur;

/**
ama ObjektuNodo objektua sortu adi fitxategitik
**/

cout << "Objektuaren .adi fitxategiaren izena? " << endl;
cin >> adi;

ObjektuAdi ama(adi); /** Objektu sortzailea **/

unsigned int maizkop=ama.maizkop;
unsigned int nhmoz;

/**
kume ObjektuPuntu objektua sortu pun fitxategitik
**/

cout << "Bolumeneko puntu hodeia daukan .pun fitxategiaren izena? " << endl;
cin >> pun;
//ObjektuPuntu kume(pun, 0, 1, -1, 0); /** Neuk sortutako .pun fitxategia **/
ObjektuPuntu kume(pun, 0, 1, -1); /** GiDek sortutako .pun fitxategia **/

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
        resfitx << "Presioa " << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        //resfitx << "Presioa " << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        resfitx << fixed << setprecision(2) << "P(Ex) f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "P(Ap) f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Error f = " << f << "Hz." << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            kserie.clear();
            dkserie.clear();
            ffinserie(ama.zentroa, kume.puntuak[puntua], k, nhmoz, kserie, dkserie);
            balio=0.0;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                balio=balio+ama.koef[maizkont][kont]*kserie[kont][0];
            }

            baliohur=balio;
            scattered(kume.puntuak[puntua], k, 15, baliozehatz);

            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(baliozehatz) << "\t" << abs(baliohur) << "\t" << 100*abs(baliohur-baliozehatz)/abs(baliozehatz) << endl;
            //resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(baliohur) << "\t" << abs(baliozehatz) << "\t" << 100*abs(abs(baliohur)-abs(baliozehatz)/abs(baliozehatz)) << endl;
            //resfitx << fixed << setprecision(10) << nodo+1 << "\t" << (180/PI)*atan2(imag(baliohur),real(baliohur)) << "\t" << (180/PI)*atan2(imag(baliozehatz),real(baliozehatz)) << "\t" << 100*(atan2(imag(baliohur),real(baliohur))-atan2(imag(baliozehatz),real(baliozehatz)))/abs(atan2(imag(baliozehatz),real(baliozehatz))) << endl;
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
        resfitx << "Presioa " << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        //resfitx << "Presioa " << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        resfitx << fixed << setprecision(2) << "P(Ex) f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "P(Ap) f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Error f = " << f << "Hz." << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            kserie.clear();
            dkserie.clear();
            finfserie(ama.zentroa, kume.puntuak[puntua], k, nhmoz, kserie, dkserie);
            balio=0.0;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                balio=balio+ama.koef[maizkont][kont]*kserie[kont][0];
            }

            baliohur=balio;
            scattered(kume.puntuak[puntua], k, 15, baliozehatz);

            baliozehatzdb=20*log10(abs(baliozehatz)/2.e-5);
            baliohurdb   =20*log10(abs(baliohur)/2.e-5);

            argzehatz=(180/PI)*atan2(imag(baliozehatz),real(baliozehatz));
            arghur=(180/PI)*atan2(imag(baliohur),real(baliohur));

            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(baliozehatz) << "\t" << abs(baliohur) << "\t" << 100*abs(baliohur-baliozehatz)/abs(baliozehatz) << endl;
            //resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(baliohur) << "\t" << abs(baliozehatz) << "\t" << 100*abs(baliohur-baliozehatz)/abs(baliozehatz) << endl;
            //resfitx << fixed << setprecision(10) << puntua+1 << "\t" << baliozehatzdb << "\t" << baliohurdb << "\t" << 100*abs(abs(baliozehatzdb)-abs(baliohurdb))/abs(baliozehatzdb) << "\t" << endl;
            //resfitx << fixed << setprecision(10) << puntua+1 << "\t" << baliozehatzdb << "\t" << baliohurdb << "\t" << abs(baliozehatzdb-baliohurdb) << "\t" << endl;
//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << argzehatz << "\t" << arghur << "\t" << 100*(abs(argzehatz)-abs(arghur))/abs(argzehatz) << endl;
        }

    }
    resfitx.close();
} /** Fluido Infinitu (FINF) objektuaren amaiera**/

return 0;
}

