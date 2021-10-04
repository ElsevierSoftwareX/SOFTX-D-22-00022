#include <iomanip>
#include <sstream>
#include <math.h>
#include <complex>
#include "../erkideak/elasserie.h"
#include "../erkideak/ffinserie.h"
#include "../erkideak/finfserie.h"
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

double f, k, young, poisson, mu, rho, cp, kp, cs, ks, lambda;
vector<complex<double> > gradientea; /** Nodoko Gradientea sartzeko **/
complex<double> pbal;
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

        young=ama.mat[0]; /** Young modulua **/
        poisson=ama.mat[1]; /** Poisson ratioa **/
        rho=ama.mat[2]; /** Dentsitatea **/

        lambda=(young*poisson)/((1-2*poisson)*(1+poisson)); /** Lameren 1. parametroa **/
        mu=young/(1+poisson);                               /** Lameren 2. parametroa **/
        cp=sqrt((lambda+2*poisson)/(rho)); /** 1. uhin-abiadura **/
        kp=2*PI*f/cp; /** P uhinaren uhin-zenbakia k1=2*pi*f/cp **/
        cs=sqrt(mu/rho); /** 2. uhin-abiadura **/
        ks=2*PI*f/cs; /** S uhinaren uhin-zenbakia k1=2*pi*f/cs **/

        resfitx << "ResultGroup" << " " << (char)34 << pun << (char)34 << " " << maizkont+1 << " OnNodes" << endl;
        resfitx << "ResultDescription" << " " << (char)34 << "Nodal Displacements" << (char)34 << " " << "Vector" << endl;
        resfitx << "ComponentNames ";
        resfitx << (char)34 << "ux f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "uy f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "uz f= " << setprecision(2) << f << " Hz." << (char)34 << ", " << endl;

        resfitx << "Values" << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            serie.clear();
            dserie.clear();
            tserie.clear();
            elasserie(ama.zentroa, kume.puntuak[puntua], kp, ks, young, poisson, nhmoz, serie, dserie, tserie);

            gradientea.clear();
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                gradientea[0]=gradientea[0]+ama.koef[maizkont][kont]*serie[kont][0];
                gradientea[1]=gradientea[1]+ama.koef[maizkont][kont]*serie[kont][1];
                gradientea[2]=gradientea[2]+ama.koef[maizkont][kont]*serie[kont][2];
            }

            resfitx << setprecision(10) << scientific << puntua+1 << "\t" << abs(gradientea[0]) << "\t" << abs(gradientea[1]) << "\t" << abs(gradientea[2]) << endl;
        }
        resfitx << "end values" << endl;
    }
    resfitx.close();

} /** Egitura Elastiko (ELAS) objektuaren amaiera **/

else if (ama.mota.compare("FFIN")==0)
{ /** Fluido Finitu (FFIN) objektuaren prozesamendua**/

    ofstream resfitx(postres.c_str());
    resfitx << "GiD Post Results File 1.0" << endl;

    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        nhmoz=ama.shkop[maizkont];
        f=ama.espektrua[maizkont];
        k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

        resfitx << "ResultGroup" << " " << (char)34 << "ADITU " << pun << (char)34 << " " << maizkont+1 << " OnNodes" << endl;
        resfitx << "ResultDescription" << " " << (char)34 << "Pressure" << (char)34 << " " << "Vector" << endl;
        resfitx << "ComponentNames ";
        resfitx << (char)34 << "P(Pa) f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "P(dB) f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "Phase f= " << setprecision(2) << f << " Hz." << (char)34 << ", " << endl;

        resfitx << "Values" << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            serie.clear();
            dserie.clear();
            ffinserie(ama.zentroa, kume.puntuak[puntua], k, nhmoz, serie, dserie);

            gradientea.clear();
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);

            pbal=0.0+0.0*I;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                gradientea[0]=gradientea[0]+ama.koef[maizkont][kont]*dserie[kont][0];
                gradientea[1]=gradientea[1]+ama.koef[maizkont][kont]*dserie[kont][1];
                gradientea[2]=gradientea[2]+ama.koef[maizkont][kont]*dserie[kont][2];
                pbal=pbal+ama.koef[maizkont][kont]*serie[kont][0];
            }

//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << real(gradientea[0]) << "\t" << real(gradientea[1]) << "\t" << real(gradientea[2]) << endl;
            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(pbal) << "\t" << 20*log10(abs(pbal/2.e-5)) << "\t" << 180/PI*atan2(imag(pbal),real(pbal)) << endl;
        }
        resfitx << "End Values" << endl;
    }
    resfitx.close();

} /** Fluido Finitu (FFIN) objektuaren amaiera**/

else if (ama.mota.compare("FINF")==0)
{ /** Fluido Infinitu (FINF) objektuaren prozesamendua**/

    ofstream resfitx(postres.c_str());
    resfitx << "GiD Post Results File 1.0" << endl;

    for (maizkont=0; maizkont<maizkop; maizkont++)
    {
        nhmoz=ama.shkop[maizkont];
        double f=ama.espektrua[maizkont];
        double k=2*PI*f/ama.mat[2]; /** uhin-zenbakia k=2*pi*f/c **/

        resfitx << "ResultGroup" << " " << (char)34 << "ADITU " << pun << (char)34 << " " << maizkont+1 << " OnNodes" << endl;
        resfitx << "ResultDescription" << " " << (char)34 << "Pressure" << (char)34 << " " << "Vector" << endl;
        resfitx << "ComponentNames ";
        resfitx << (char)34 << "P(Pa) f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "P(dB) f= " << setprecision(2) << f << " Hz." << (char)34 << ", ";
        resfitx << (char)34 << "Phase f= " << setprecision(2) << f << " Hz." << (char)34 << ", " << endl;

        resfitx << "Values" << endl;

        for (puntua=0; puntua<kume.puntukop; puntua++)
        {
            serie.clear();
            dserie.clear();
            finfserie(ama.zentroa, kume.puntuak[puntua], k, nhmoz, serie, dserie);

            gradientea.clear();
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);
            gradientea.push_back(0.0+0.0*I);

            pbal=0.0+0.0*I;

            for (kont=0; kont<ama.koefkop[maizkont]; kont++)
            {
                // cout << "KOEF[" << maizkont << "," << kont << "]= " << ama.koef[maizkont][kont] << endl;

                gradientea[0]=gradientea[0]+ama.koef[maizkont][kont]*dserie[kont][0];
                gradientea[1]=gradientea[1]+ama.koef[maizkont][kont]*dserie[kont][1];
                gradientea[2]=gradientea[2]+ama.koef[maizkont][kont]*dserie[kont][2];

                pbal=pbal+ama.koef[maizkont][kont]*serie[kont][0];
            }

//            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << real(gradientea[0]) << "\t" << real(gradientea[1]) << "\t" << real(gradientea[2]) << endl;
            resfitx << fixed << setprecision(10) << puntua+1 << "\t" << abs(pbal) << "\t" << 20*log10(abs(pbal/2.e-5)) << "\t" << 180/PI/atan2(imag(pbal),real(pbal)) << endl;
        }
        resfitx << "End Values" << endl;
    }

    resfitx.close();

} /** Fluido Infinitu (FINF) objektuaren amaiera**/

return 0;
}
