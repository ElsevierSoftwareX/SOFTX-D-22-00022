#include "ObjektuAdi.h"
#include <sstream>
#include <math.h>

using namespace std;

const complex<double> I{0.0, 1.0};

/** ObjektuAdi hutsa sortzeko eraikitzailea (binarioa iraultzeko) **/

ObjektuAdi::ObjektuAdi()
{

}

/** ADITUk sortutako .adi testu-fitxategi batetik eraikitzeko **/

ObjektuAdi::ObjektuAdi(string nom)
{
    string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
    string fitxategia;
    string lerroa; /** Irakurritako lerro osoa **/

    vector<double> moduluak; /** GiDek emandako bektore normalen moduluak (normalak unitario bihurtzeko) **/
    vector<unsigned int> itemp; /** Lerroak txertatzeko bitarteko bektorea **/
    vector<double> dtemp; /** Lerroak txertatzeko bitarteko bektorea **/

    izena = nom;

    string txartela="HASI";

    unsigned int ikok, jkok, tam; /** Irakurtzeko laguntzaileak **/
    double bik; /**Irakurtzeko double **/

    float err;
    complex<double> cik;

    fitxategia = bidea + izena + ".adi";
    ifstream adifitx(fitxategia.c_str(), ios::in);

    getline(adifitx,lerroa); /** Objektuaren proiektua **/
    istringstream ls(lerroa);
    ls >> proiektua;

    getline(adifitx,lerroa); /** Objektu mota **/
    istringstream ls1(lerroa);
    ls1 >> mota;

    getline(adifitx,lerroa); /** Maiztasun kopurua **/
    istringstream ls3(lerroa);
    ls3 >> maizkop;

    getline(adifitx,lerroa); /** Koefiziente kopuru maximoa **/
    istringstream ls4(lerroa);
    ls4 >> tam;

    espektrua.clear();
    espektrua.resize(maizkop, 0.0);
    honderl.clear();
    honderl.resize(maizkop, 0.0);
    koef.clear();
    koef.resize(maizkop, vector<complex<double>>(tam,0.0+0.0*I));

    while(txartela.compare("AMAI") != 0)
    {
    getline(adifitx,lerroa);
    istringstream ls(lerroa);
    ls >> txartela;

    if (txartela.compare("ZENT")==0)
        {
        ls >> bik;
        zentroa.push_back(bik);
        ls >> bik;
        zentroa.push_back(bik);
        ls >> bik;
        zentroa.push_back(bik);
        }

    else if (txartela.compare("MAT1")==0)
        {
        ls >> bik;
        mat.push_back(bik);
        ls >> bik;
        mat.push_back(bik);
        ls >> bik;
        mat.push_back(bik);
        ls >> bik;
        mat.push_back(bik);
        }

    else if (txartela.compare("FREQ")==0)
        {
        ls >> ikok;
        ls >> bik;
        espektrua[ikok-1]=bik;
        ls >> ikok;
        shkop.push_back(ikok);
        ls >> ikok;
        koefkop.push_back(ikok);
        ls >> err;
        honderl.push_back(err);
        }

    else if (txartela.compare("KOEF")==0)
        {
        ls >> ikok;
        ls >> jkok;
        ls >> cik;
        koef[ikok-1][jkok-1]=cik;
        }
    } /* of while txartela */

    adifitx.close();

} /** eraikitzailearen amaiera **/

/** Desegilea **/
ObjektuAdi::~ObjektuAdi()
{

}

