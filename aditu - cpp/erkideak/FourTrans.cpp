#include "FourTrans.h"
#include <sstream>
#include <math.h>
#include <complex>

using namespace std;

const complex<double> I{0.0, 1.0};

/** string: Fourier Transformatua daukan .tau testu-fitxategia **/

FourTrans::FourTrans(string nom)
{
    string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
    string fitxategia;
    string lerroa; /** Irakurritako lerro osoa **/

    vector<complex<double> > cbek; /** Lerroak txertatzeko bitarteko bektorea **/

    string txartela;

    proiektua = nom;

    unsigned int maiz, kok; /** Irakurtzeko laguntzaileak **/
    double bik; /**Irakurtzeko double **/
    complex<double> cik;

    fitxategia = bidea + nom + ".tau";
    ifstream taufitx(fitxategia.c_str(), ios::in);

    getline(taufitx,lerroa); /** Maiztasun kopurua **/
    istringstream ls(lerroa);
    ls >> maizkop;

    getline(taufitx,lerroa); /** Lerro kopurua **/
    istringstream ls2(lerroa);
    ls2 >> lerrokop;

    espektrua.clear();
    taula.clear();

    for (maiz=0; maiz < maizkop; maiz++)
    {
        getline(taufitx,lerroa);
        istringstream ls(lerroa);
        ls >> txartela;
        ls >> bik;
        espektrua.push_back(bik);

        cbek.clear();
        for (kok=0; kok < lerrokop; kok++)
        {
            getline(taufitx,lerroa);
            istringstream ls(lerroa);
            ls >> cik;
            cbek.push_back(cik);
        }
        taula.push_back(cbek);
    }

    taufitx.close();

} /** eraikitzailearen amaiera **/

/** Desegilea **/
FourTrans::~FourTrans()
{

}

