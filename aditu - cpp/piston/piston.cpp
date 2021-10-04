#include <iomanip>
#include <sstream>
#include <math.h>
#include <complex>
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
#include "../erkideak/ObjektuAurre.h"
#include "../erkideak/ObjektuPost.h"
#include "pistonff.h"
const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

using namespace std;

int main()
{

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
string adi, mub, flaviares;

unsigned int maizkont, kont, nodo;

complex<double> balio;
complex<double> baliozehatz, baliohur;
double lambda;

/**
ama ObjektuPost objektua sortu adi fitxategitik
**/

cout << "Objektuaren .adi fitxategiaren izena? " << endl;
cin >> adi;

ObjektuPost ama(adi); /** Objektu sortzailea **/

vector<double> zentro;
zentro.push_back(ama.zentroa[0]);
zentro.push_back(ama.zentroa[1]);
zentro.push_back(ama.zentroa[2]);

unsigned int maizkop=ama.maizkop;
unsigned int term=ama.luzera;
unsigned int nhmoz;

/**
kume ObjektuAurre objektua sortu mub fitxategitik
**/

cout << "Bolumeneko puntu hodeia daukan .mub fitxategiaren izena? " << endl;
cin >> mub;
ObjektuAurre kume(mub, 0, -1);

/**
flavia.res fitxategia idatzi
**/

flaviares = bidea + mub + ".flavia.res";

unsigned int nodokop=kume.nodoak.size();

vector<double> xnodo(3, 0.0); /** nodoaren koordenatuak sartzeko **/

vector<complex<double> > kserie; /** Serie konplexua **/
vector<vector<complex<double> > > dkserie; /** Serie konplexuaren deribatuak cartesiarretan **/

if (ama.mota.compare("ELAS")==0)

/************************* EGITURA ELASTIKOA ********************/
{ /** Egitura elastiko objektuaren prozesamendua **/


} /** Egitura Elastiko (ELAS) objektuaren amaiera **/

else if (ama.mota.compare("FFIN")==0)
{ /** Fluido Finitu (FFIN) objektuaren prozesamendua**/


} /** Fluido Finitu (FFIN) objektuaren amaiera**/

else if (ama.mota.compare("FINF")==0)
{ /** Fluido Infinitu (FINF) objektuaren prozesamendua**/

    nhmoz=floor(sqrt(term));
    ofstream resfitx(flaviares.c_str());
    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        double f=ama.espektrua[maizkont];
        double k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/
        lambda=ama.mat[2]/f;
        resfitx << setfill('0');
        resfitx << "Presioa " << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        //resfitx << "***** error = " << fixed << setprecision(2) << ama.honderl[maizkont] << " %" << endl;
        resfitx << fixed << setprecision(2) << "P(Ap) f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "P(Ex) f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Diff  f = " << f << "Hz." << endl;

        for (nodo=0; nodo<nodokop; nodo++)
        {
            xnodo[0]=kume.nodoak[nodo][0]*lambda*6;
            xnodo[1]=kume.nodoak[nodo][1]*lambda*6;
            xnodo[2]=kume.nodoak[nodo][2]*lambda*6;

            kserie.clear();
            dkserie.clear();
            finfserie(zentro, xnodo, k, nhmoz, kserie, dkserie);
            balio=0.0+0.0*I;

            for (kont=0; kont<term; kont++)
            {
                balio=balio+ama.koef[maizkont][kont]*kserie[kont];
            }

            baliohur=balio;
//            baliolog=20*log10(baliomod/2.e-5);
//            baliofase=(180/PI)*atan2(imag(balio),real(balio));

            pistonff(xnodo, k, 0.1, baliozehatz);

            resfitx << fixed << setprecision(10) << nodo+1 << "\t" << abs(baliohur) << "\t" << abs(baliozehatz) << "\t" << 100*abs(baliozehatz-baliohur)/abs(baliozehatz) << endl;
        }

    }
    resfitx.close();
} /** Fluido Infinitu (FINF) objektuaren amaiera**/

return 0;
}
