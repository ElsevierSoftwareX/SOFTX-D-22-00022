#include <iomanip>
#include <sstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <complex>
#ifndef PI
    #define PI 3.14159265358979323846
#endif
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
#include "../erkideak/ObjektuAurre.h"
#include "../erkideak/ObjektuPost.h"

using namespace std;

const complex<double> I{0.0, 1.0};

int main()
{

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
string adi, mub, flaviares;

unsigned int maizkont, kont, nodo;

vector<complex<double> > neumann; /** Nodoko Gradientea sartzeko **/
complex<double> hbal;
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
vector<double> normal(3, 0.0); /** nodoko normala sartzeko **/

vector<vector<complex<double> > >  kserie; /** Serie konplexua **/
vector<vector<complex<double> > > dkserie; /** Serie konplexuaren deribatuak cartesiarretan **/

if (ama.mota.compare("ELAS")==0)

/************************* EGITURA ELASTIKOA ********************/
{ /** Egitura elastiko objektuaren prozesamendua **/


} /** Egitura Elastiko (ELAS) objektuaren amaiera **/

else if (ama.mota.compare("FFIN")==0)
{ /** Fluido Finitu (FFIN) objektuaren prozesamendua**/

    nhmoz=floor(sqrt(term/2));
    ofstream resfitx(flaviares.c_str());
    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        double f=ama.espektrua[maizkont];
        double k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/
        resfitx << setfill('0');
        resfitx << "Normals" << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        //resfitx << "# error = " << fixed << setprecision(2) << ama.honderl[maizkont] << endl;
        resfitx << fixed << setprecision(2) << "Nx    f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Ny    f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Nz    f = " << f << "Hz." << endl;

        for (nodo=0; nodo<nodokop; nodo++)
        {
            xnodo[0]=kume.nodoak[nodo][0];
            xnodo[1]=kume.nodoak[nodo][1];
            xnodo[2]=kume.nodoak[nodo][2];

            normal[0]=kume.normalak[nodo][0];
            normal[1]=kume.normalak[nodo][1];
            normal[2]=kume.normalak[nodo][2];

            kserie.clear();
            dkserie.clear();
            ffinserie(zentro, xnodo, k, nhmoz, kserie, dkserie);

            neumann.clear();
            neumann.push_back(0.0+0.0*I);
            neumann.push_back(0.0+0.0*I);
            neumann.push_back(0.0+0.0*I);

            hbal=0.0+0.0*I;

            for (kont=0; kont<term; kont++)
            {
                //cout << "Grad = (" <<dkserie[kont][0] << " , " <<dkserie[kont][1] << " , " <<dkserie[kont][2] << " )" << endl;

                neumann[0]=neumann[0]+ama.koef[maizkont][kont]*dkserie[kont][0];
                neumann[1]=neumann[1]+ama.koef[maizkont][kont]*dkserie[kont][1];
                neumann[2]=neumann[2]+ama.koef[maizkont][kont]*dkserie[kont][2];
                hbal=hbal+ama.koef[maizkont][kont]*(dkserie[kont][0]*normal[0]+dkserie[kont][1]*normal[1]+dkserie[kont][2]*normal[2]);
            }

            resfitx << fixed << setprecision(10) << nodo+1 << "\t" << normal[0] << "\t" << normal[1] << "\t" << normal[2] << endl;
        }

    }
    resfitx.close();
} /** Fluido Finitu (FFIN) objektuaren amaiera**/

else if (ama.mota.compare("FINF")==0)
{ /** Fluido Infinitu (FINF) objektuaren prozesamendua**/

    nhmoz=floor(sqrt(term));
    ofstream resfitx(flaviares.c_str());
    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        double f=ama.espektrua[maizkont];
        double k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/
        resfitx << setfill('0');
        resfitx << "Normals" << fixed << setprecision(0) << setw(5) << f << "Hz" << fixed << setprecision(0) << "    3    1    2    1    1" << endl;
        //resfitx << "# error = " << fixed << setprecision(2) << ama.honderl[maizkont] << endl;
        resfitx << fixed << setprecision(2) << "Nx    f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Ny    f = " << f << "Hz." << endl;
        resfitx << fixed << setprecision(2) << "Nz    f = " << f << "Hz." << endl;

        for (nodo=0; nodo<nodokop; nodo++)
        {
            xnodo[0]=kume.nodoak[nodo][0];
            xnodo[1]=kume.nodoak[nodo][1];
            xnodo[2]=kume.nodoak[nodo][2];

            normal[0]=kume.normalak[nodo][0];
            normal[1]=kume.normalak[nodo][1];
            normal[2]=kume.normalak[nodo][2];

            kserie.clear();
            dkserie.clear();
            finfserie(zentro, xnodo, k, nhmoz, kserie, dkserie);

            neumann.clear();
            neumann.push_back(0.0+0.0*I);
            neumann.push_back(0.0+0.0*I);
            neumann.push_back(0.0+0.0*I);

            hbal=0.0+0.0*I;

            for (kont=0; kont<term; kont++)
            {

                neumann[0]=neumann[0]+ama.koef[maizkont][kont]*dkserie[kont][0];
                neumann[1]=neumann[1]+ama.koef[maizkont][kont]*dkserie[kont][1];
                neumann[2]=neumann[2]+ama.koef[maizkont][kont]*dkserie[kont][2];
                hbal=hbal+ama.koef[maizkont][kont]*(dkserie[kont][0]*normal[0]+dkserie[kont][1]*normal[1]+dkserie[kont][2]*normal[2]);

            }

            resfitx << fixed << setprecision(10) << nodo+1 << "\t" << normal[0] << "\t" << normal[1] << "\t" << normal[2] << endl;
        }

    }
    resfitx.close();
} /** Fluido Infinitu (FINF) objektuaren amaiera**/

return 0;
}
