#include "ObjektuPuntu.h"
#include <sstream>
#include <math.h>

using namespace std;

const complex<double> I{0.0, 1.0};

/** ObjektuPuntu hutsa sotrtzeko eraikitzailea (binarioa iraultzeko) **/

ObjektuPuntu::ObjektuPuntu()
{

}
/** string: neuk sortutako .pun testu-fitxategia **/

ObjektuPuntu::ObjektuPuntu(string nom, unsigned int zenb, unsigned int objkop, int kok, int neuk)
{
    unsigned int kont;

    string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
    string fitxategia;
    string lerroa; /** Irakurritako lerro osoa **/

    vector<double> moduluak; /** GiDek emandako bektore normalen moduluak (normalak unitario bihurtzeko) **/
    vector<unsigned int> itemp; /** Lerroak txertatzeko bitarteko bektorea **/
    vector<double> dtemp; /** Lerroak txertatzeko bitarteko bektorea **/

    ordinala = zenb;

    akopla.resize(objkop, 0);

    kok1=kok+1;
    kok2=kok;

    izena = nom;

    string txartela="HASI";

    unsigned int etiketa; /** Irakurtzeko laguntzaileak **/
    unsigned int osoa, osoa2; /** Irakurtzeko int **/
    double bik; /**Irakurtzeko double **/

    fitxategia = bidea + izena + ".pun";
    ifstream adifitx(fitxategia.c_str(), ios::in);

    getline(adifitx,lerroa); /** ObjektuAurre-mota **/
    istringstream ls(lerroa);
    ls >> mota;

    getline(adifitx,lerroa); /** Nodo kopurua **/
    istringstream ls2(lerroa);
    ls2 >> puntukop;
    luzera=puntukop;

    mbzenb.resize(puntukop,0);
    mbmota.resize(puntukop,0);
    mbeku.resize(puntukop,0);
    mbobj.resize(puntukop, vector<int> (3, 0));

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

    else if (txartela.compare("GRID")==0)
        {
        ls >> etiketa;

        dtemp.clear();
        ls >> bik;
        dtemp.push_back(bik);
        ls >> bik;
        dtemp.push_back(bik);
        ls >> bik;
        dtemp.push_back(bik);
        puntuak.push_back(dtemp);
        }

    else if (txartela.compare("MBAL")==0)
        {
        ls >> etiketa;

        ls >> osoa;
        mbzenb[etiketa-1]=osoa;
        ls >> osoa;
        mbmota[etiketa-1]=osoa;
        ls >> osoa;
        mbeku[etiketa-1]=osoa;
        for (kont=0; kont< osoa; kont++)
        {
            ls >> osoa2;
            mbobj[etiketa-1][kont]=osoa2-1;
        }
        }

    else if (txartela.compare("NORM")==0)
        {
        ls >> etiketa;

        dtemp.clear();
        ls >> bik;
        dtemp.push_back(bik);
        ls >> bik;
        dtemp.push_back(bik);
        ls >> bik;
        dtemp.push_back(bik);
        normalak.push_back(dtemp);
        }

    } /** of while txartela **/

//    for (kont=0; kont < puntukop; kont++)
//        {
//        kok2=kok2+mbeku[kont];
//        }
//
//    /** Baldintza akoplamenduak matrizean gehitu **/
//
//    for (kont=0; kont < puntukop; kont++)
//    {
//        if (mbmota[kont]==0)
//        {
//            akopla[zenb]=akopla[zenb]+mbeku[zenb];
//        }
//        else
//        {
//            akopla[mbobj[kont][0]]=akopla[mbobj[kont][0]]+mbeku[zenb];
//        }
//    }

    adifitx.close();

} /** eraikitzailearen amaiera **/

/** string: GiDek sortutako .pun testu-fitxategia **/

ObjektuPuntu::ObjektuPuntu(string nom, unsigned int zenb, unsigned int objkop, int kok)
{
    unsigned int kont;
    double norma;

    string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
    string fitxategia;
    string lerroa; /** Irakurritako lerro osoa **/

    vector<double> moduluak; /** GiDek emandako bektore normalen moduluak (normalak unitario bihurtzeko) **/
    vector<unsigned int> itemp; /** Lerroak txertatzeko bitarteko bektorea **/
    vector<double> dtemp; /** Lerroak txertatzeko bitarteko bektorea **/

    ordinala = 0;

    akopla.resize(objkop, 0);

    kok1=kok+1;
    kok2=kok;

    izena = nom;

    string txartela="HASI";

    unsigned int etiketa, nodele; /** Irakurtzeko laguntzaileak **/
    unsigned int osoa, osoa2; /** Irakurtzeko int **/
    double bik; /**Irakurtzeko double **/

    fitxategia = bidea + izena + ".pun";
    ifstream adifitx(fitxategia.c_str(), ios::in);

    getline(adifitx,lerroa); /** ObjektuAurre-mota **/
    istringstream ls(lerroa);
    ls >> mota;

    getline(adifitx,lerroa); /** Nodo kopurua **/
    istringstream ls2(lerroa);
    ls2 >> puntukop;
    luzera=puntukop;

    mbzenb.resize(puntukop,0);
    mbmota.resize(puntukop,0);
    mbeku.resize(puntukop,0);
    mbobj.resize(puntukop, vector<int> (3, 0));

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

    else if (txartela.compare("GRID")==0)
        {
        ls >> etiketa;

        ls >> bik;
        dtemp.push_back(bik);
        ls >> bik;
        dtemp.push_back(bik);
        ls >> bik;
        dtemp.push_back(bik);
        puntuak.push_back(dtemp);
        dtemp.clear();

        moduluak.push_back(0.0);

        dtemp.push_back(0.0);
        dtemp.push_back(0.0);
        dtemp.push_back(0.0);
        normalak.push_back(dtemp);
        dtemp.clear();
        }

    else if (txartela.compare("MBAL")==0)
        {
        ls >> etiketa;

        ls >> osoa;
        mbzenb[etiketa-1]=osoa; /** Mugalde baldintzaren etiketa-zenbakia **/
        ls >> osoa;
        mbmota[etiketa-1]=osoa; /** Mugalde baldintza mota: 0: barnekoa; 1; akoplamendua **/
        ls >> osoa;
        mbeku[etiketa-1]=osoa;  /** Mugalde baldintzak ezartzen duen ekuazio kopurua **/

        for (kont=0; kont< osoa; kont++)
        {
            ls >> osoa2;
            mbobj[etiketa-1][kont]=osoa2-1;
        }

        }

    else if (txartela.compare("ELEM")==0)
        {
        ls >> etiketa;
        ls >> nodele;

        for (kont=0; kont < nodele; kont++)
            {
            ls >> etiketa;
            itemp.push_back(etiketa-1);
            }
        for (kont=0; kont < 3; kont++)
            {
            ls >> bik;
            dtemp.push_back(bik);
            }

        for (kont=0; kont < nodele; kont++)
        {
            /** Nodoak dagoeneko normalik badu**/
            norma=sqrt(normalak[itemp[kont]][0]*normalak[itemp[kont]][0]+normalak[itemp[kont]][1]*normalak[itemp[kont]][1]+normalak[itemp[kont]][2]*normalak[itemp[kont]][2]);
            if (norma > 0.5)
            {
                normalak[itemp[kont]][0]=normalak[itemp[kont]][0]/norma;
                normalak[itemp[kont]][1]=normalak[itemp[kont]][1]/norma;
                normalak[itemp[kont]][2]=normalak[itemp[kont]][2]/norma;

                /** Egiaztatu normal berria batezbestekoan sartuko den (biderkarura eskalarra > balioa bada) **/
                if (normalak[itemp[kont]][0]*dtemp[0]+normalak[itemp[kont]][1]*dtemp[1]+normalak[itemp[kont]][2]*dtemp[2] > 0.7)
                {
                    normalak[itemp[kont]][0]=normalak[itemp[kont]][0]+dtemp[0];
                    normalak[itemp[kont]][1]=normalak[itemp[kont]][1]+dtemp[1];
                    normalak[itemp[kont]][2]=normalak[itemp[kont]][2]+dtemp[2];
                }
            }
            else /** Lehenengo normala **/
            {
                normalak[itemp[kont]][0]=dtemp[0];
                normalak[itemp[kont]][1]=dtemp[1];
                normalak[itemp[kont]][2]=dtemp[2];
            }
        }
        itemp.clear();
        dtemp.clear();
        }

    else if (txartela.compare("AMAI")==0)
        {

        }
    } /* of while txartela */

    for (kont=0; kont < puntukop; kont++)
        {
        moduluak[kont]=sqrt(normalak[kont][0]*normalak[kont][0]+normalak[kont][1]*normalak[kont][1]+normalak[kont][2]*normalak[kont][2]);
        }

    for (kont=0; kont < puntukop; kont++)
        {
        normalak[kont][0]=normalak[kont][0]/moduluak[kont];
        normalak[kont][1]=normalak[kont][1]/moduluak[kont];
        normalak[kont][2]=normalak[kont][2]/moduluak[kont];
        }

//    for (kont=0; kont < puntukop; kont++)
//        {
//        kok2=kok2+mbeku[kont];
//        }
//
//    /** Baldintza akoplamenduak matrizean gehitu **/
//
//    for (kont=0; kont < puntukop; kont++)
//    {
//        if (mbmota[kont]==0)
//        {
//            akopla[zenb]=akopla[zenb]+mbeku[zenb];
//        }
//        else
//        {
//            akopla[mbobj[kont][0]]=akopla[mbobj[kont][0]]+mbeku[zenb];
//        }
//    }

    adifitx.close();

} /** eraikitzailearen amaiera **/

/** Desegilea **/
ObjektuPuntu::~ObjektuPuntu()
{

}

