#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <vector>
#include <complex>
#include "../erkideak/FourTrans.h"
#include "../erkideak/ObjektuPuntu.h"
#include "doiketa.h"
#include "kaminlu.h"
#include "siskalk.h"
#include "../scattering/neumann.h"

using namespace std;

int main()
{
    string proiektua_oro, taula;
    vector<float> honderl_oro; /** maiztasun bakoitzean egindako errore erlatiboa (%) **/

    int amkop;
    string hitza; /** laguntzailea **/
    vector <string> izenak;
    vector<vector<unsigned int> > indizeak;

    vector<unsigned int> itemp; /** Laguntzaileak **/
    vector<double> dtemp; /** Laguntzailea **/
    vector<complex<double> > ctemp; /** Laguntzailea **/

    vector<ObjektuPuntu> puntu;

    int proharmonikomax;

    unsigned int kont;
    unsigned int maizkont, objin, objout;
    int kok;

    cout << "Proiektuaren izena? " << endl;
    cin >> proiektua_oro;

    cout << "Fourier Transformatuak dauzkan .tau fitxategia? " << endl;
    cin >> taula;

    /************************* Fourier transformatua  objektua sortu *************************************/

    FourTrans ft(taula);

    /************************************************************************************/
    /************************* OBJEKTUEN IRAKURKETA *************************************/
    /************************************************************************************/

    /***** ObjektuPuntu klasekoak puntu bektorean gordeko dira *****/

    unsigned int objkop;
    cout << "Objektu kopurua guztira? " << endl;
    cin >> objkop;

    for (objin=0; objin<objkop; objin++)
    {
        cout << objin+1 << " Objektuaren izena? " << endl;
        cin >> hitza;
        izenak.push_back(hitza); /** Objektuaren izena **/

        if (objin==0)
            {
                kok=-1; /** objont=0 kasurako bakarrik **/
            }
            else
            {
                kok=indizeak[objin-1][1];
            }

        ObjektuPuntu alepuntu(hitza, objin, objkop, kok, 0);    /** Neuk sortutako .pun fitxategi batetik **/
        //ObjektuPuntu alepuntu(hitza, objin, objkop, kok);        /** GiDek sortutako .pun fitxategi batetik **/

        itemp.push_back(alepuntu.kok1);
        itemp.push_back(alepuntu.kok2);
        indizeak.push_back(itemp);
        itemp.clear();

        /**alepuntu objektua puntu bektoreari erantsi **/

        puntu.push_back(alepuntu);

        cout << puntu[objin].izena << " objektua irakurrita." << endl;

    } /**objin begiztaren (irakurketaren) amaiera **/

    /************************************************************************************/
    /*************************** AKOPLAMENDUEN DOIKETA **********************************/
    /***************************************************PI*********************************/

    vector<vector<unsigned int> > ekkop;
    ekkop.clear();
    ekkop.resize(objkop, vector<unsigned int>(objkop,0));

    doiketa(puntu, ekkop, objkop);

    for (objin=0; objin < objkop; objin++)
    {
           for (objout=0; objout < objkop; objout++)
           {
               cout << "akopla[" << objin << "],[" << objout << "]= " << puntu[objin].akopla[objout] << endl;
           }
    }

    for (objin=0; objin < objkop; objin++)
    {
        if (puntu[objin].mota=="ELAS")
        {
            amkop=puntu[objin].kok2-puntu[objin].kok1+1;
            puntu[objin].harmonikomax=floor(sqrt(amkop/8-1));
        } /** ELAS **/

        if (puntu[objin].mota=="OSKO")
        {
            amkop=puntu[objin].kok2-puntu[objin].kok1+1;
            puntu[objin].harmonikomax=floor(amkop/8-1);
        } /** OSKO **/

        if (puntu[objin].mota=="FFIN")
        {
            amkop=puntu[objin].kok2-puntu[objin].kok1+1;
            puntu[objin].harmonikomax=floor(sqrt(amkop/2-1));
        } /** FFIN **/

        if (puntu[objin].mota=="FINF")
        {
            amkop=puntu[objin].kok2-puntu[objin].kok1+1;
            puntu[objin].harmonikomax=floor(sqrt(amkop-1));
        } /** FINF **/

        for (objout=0; objout < objkop; objout++)
        {
            cout << "ekkop[" << objin << "],[" << objout << "]= " << ekkop[objin][objout] << endl;
        }
    }

    kok=0;
    for (objin=0; objin < objkop; objin++)
    {
        cout << puntu[objin].izena << " objektuaren harmoniko maximoa: " << puntu[objin].harmonikomax << endl;

        if (puntu[objin].harmonikomax > kok)
        {
            kok=puntu[objin].harmonikomax;
        }
        cout << "Obj " << objin << ": kok1= " << puntu[objin].kok1 << ": kok2= " << puntu[objin].kok2 << endl;
    }

    proharmonikomax=kok;

    cout << endl;
    cout << "Proiektuaren harmoniko maximoa: " << proharmonikomax << endl;
    cout << endl;

    /************************************************************************************/
    /*************************** KALKULUA **********************************************/
    /************************************************************************************/

    /************************* ALDAGAI KOMUNAK OBJEKTU GUZTIETARAKO ********************/

    unsigned int neurria = puntu[objkop-1].kok2+1;

    cout << "Neurria: " << neurria << endl;

    vector<vector<complex<double> > > amat; /** Koefiziente-matrizea **/
    vector<complex<double> > bbek; /** Gai aske bektorea **/
    vector<complex<double> > ebazpena; /** Maiztasun bateko ebazpena **/

    complex<double> hbal;
    vector<complex<double> > indar; /** Egiturako indarrPIa gordetzeko **/

    int harmonikoa;
    vector<vector<unsigned int> > zutabeak /** objektu bakoitzerako, hasierako eta amaerako amat osoko zutabeak mihiztatzeko **/;

    vector<unsigned int> shkop; /** Maiztasun eta objektu bakoitzeko konbergitu duen harmoniko kopurua **/
    vector<unsigned int> koefkopmax; /** Koefiziente kopuru maximoa (azken harmonikoraino) objektu bakoitzeko maiztasun guztietarako **/
    koefkopmax.clear();
    koefkopmax.resize(objkop,0);

    vector<vector<unsigned int> > koefkop; /** Maiztasun eta objektu bakoitzeko konbergitu duen koefiziente kopurua **/

    float errorea, erroremax; /** Sistemaren ebazpeneko hondarraren errorea eta errore maximoa (%) **/
    vector<float> erroreak; /** Urrats guztietako erroreak **/
    float abiadura, abiaduramax=0.03; /** azken urratseko errorearen aldaketaren abiadura **/
    unsigned int urraskop=2;
    float aurrekoerroreak[urraskop]; /** aurreko urratsetako erroreak **/
    unsigned int prokoefkop=0;
    float gbh, gbhmax=1.0; /** RAM memoria maximoa amat matrizea gordetzeko **/

    for (maizkont=0; maizkont<ft.maizkop; maizkont++)
    {
        erroreak.clear();
        erroremax=0.05;
        cout << ft.espektrua[maizkont] << " Hz maiztasuna hasita: " << endl;
        cout << "Errore minimoa: % " << fixed << setprecision(2) << erroremax << endl;
        cout << "Abiadura minimoa: " << fixed << setprecision(2) << abiaduramax << " %/urrats" << endl;
        cout << "RAM memoria maximoa: " << fixed << setprecision(2) << gbhmax << " Gb" << endl;

        for (kont=0; kont < urraskop; kont++)
        {
            aurrekoerroreak[kont]=50000.+kont*1000000.;
        }
        errorea=1000.;
        abiadura=(aurrekoerroreak[urraskop-1]-errorea)/urraskop;

        harmonikoa=-1;

        gbh=0.0;

        //while (abs(abiadura) > abiaduramax && gbh < gbhmax)
        while (errorea > erroremax && abiadura > abiaduramax && gbh < gbhmax && harmonikoa < proharmonikomax)
        //while (errorea > erroremax && abs(abiadura) > abiaduramax && gbh < gbhmax && harmonikoa < proharmonikomax)
        {
            harmonikoa++;
            //clock_t begin = clock();
            siskalk(puntu, objkop, ekkop, harmonikoa, maizkont, taula, amat, bbek, zutabeak, prokoefkop);
            //clock_t end = clock();
            //cout << "Sistema osatzeko " << double(end - begin) / CLOCKS_PER_SEC << " segundo." << endl;

            gbh=16*neurria*prokoefkop/1.e9; /** RAM memoria (Gb) = 16 byte/elementu * amat matrizeko elementu kopurua **/

            /** Proiektu osoaren maizkont maiztasuneko sistema osoa ebatzi **/

            ebazpena.clear();
            //begin = clock();

            kaminlu(amat, bbek, zutabeak, ebazpena, errorea); /** Ekuazio-sistemaren ebazpena (LU deskonposizioa) eta errorearen kalkulua **/
            //end = clock();
            //cout << "Sistema ebazteko " << double(end - begin) / CLOCKS_PER_SEC << " segundo." << endl;
            //cout << endl;

            if (isnan(errorea))
            {
                errorea=100;
            }
            erroreak.push_back(errorea);
            abiadura=(aurrekoerroreak[urraskop-1]-errorea)/urraskop;

            cout << harmonikoa << ".eraino: Errorea: % " << fixed << setprecision(3) << errorea << "; Abiadura: " << abiadura << " (%/u)" << "; Gb:  " << gbh << endl;

            for (kont=urraskop-1; kont > 0; kont--)
            {
                aurrekoerroreak[kont]=aurrekoerroreak[kont-1];
            }
            aurrekoerroreak[0]=errorea;

        } /** end of while errorea. **/

        /*** Errorerik txikiena egin duen harmoniko kopurua aukeratu eta birkalkulatu ***/

        kok=0;
        erroremax=1000.0;

        for (kont=0;kont<erroreak.size(); kont++)
        {
            if (isnan(erroreak[kont]))
            {
                erroreak[kont]=1000.0;
            }
            else
            {
                if (erroreak[kont] < erroremax)
                {
                        kok=kont;
                        erroremax=erroreak[kont];
                }
            }
        }

        harmonikoa=kok;

        siskalk(puntu, objkop, ekkop, harmonikoa, maizkont, taula, amat, bbek, zutabeak, prokoefkop);

        ebazpena.clear();
        kaminlu(amat, bbek, zutabeak, ebazpena, errorea); /** Ekuazio-sistemaren ebazpena (LU deskonposizioa) eta errorearen kalkulua **/

        /***********************************************************************/
        /*********************  MAIZTASUNEKO KALKULUAREN AMAIERA ***************/
        /***********************************************************************/


        cout << endl;
        cout << ft.espektrua[maizkont] << " Hz maiztasunaren amaiera. Harmonikoa: " << harmonikoa << ". Errorea: % " << errorea << endl;
        cout << endl;
        honderl_oro.push_back(errorea);

        /** Koefizienteak gorde **/

        for (objin=0; objin<objkop; objin++)
        {
            puntu[objin].shkop.push_back(harmonikoa);
            if (puntu[objin].mota=="ELAS")
            {
                kont=8*pow(harmonikoa+1,2);
                puntu[objin].koefkop.push_back(kont);
                if (kont > koefkopmax[objin])
                {
                    koefkopmax[objin]=kont;
                }
            }
            else if (puntu[objin].mota=="OSKO")
            {
                kont=8*(harmonikoa+1);
                puntu[objin].koefkop.push_back(kont);
                if (kont > koefkopmax[objin])
                {
                    koefkopmax[objin]=kont;
                }
            }
            else if (puntu[objin].mota=="FFIN")
            {
                kont=2*pow(harmonikoa+1,2);
                puntu[objin].koefkop.push_back(kont);
                if (kont > koefkopmax[objin])
                {
                    koefkopmax[objin]=kont;
                }
            }
            else if (puntu[objin].mota=="FINF")
            {
                kont=pow(harmonikoa+1,2);
                puntu[objin].koefkop.push_back(kont);
                if (kont > koefkopmax[objin])
                {
                koefkopmax[objin]=kont;
                }
            }
            ctemp.clear();
            for (kont=zutabeak[objin][0]; kont <= zutabeak[objin][1]; kont++)
            {
                ctemp.push_back(ebazpena[kont]);
            }
            puntu[objin].koef.push_back(ctemp);
        }

} /** end of maizkont. Proiektuaren espektruaren amaiera **/

/***********************************************************************/
/********************  ESPEKTRU OSOKO KALKULUAREN AMAIERA **************/
/***********************************************************************/

/***********************************************************************/
/********************  TESTU FITXATEGIAK IDATZI ************************/
/***********************************************************************/

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /* proiektuko fitxategien direktorioa */
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /* proiektuko fitxategien direktorioa */

string adifitx;

for (objin=0; objin<objkop; objin++)
{
    hitza = bidea + puntu[objin].izena + ".adi";
    ofstream adifitx(hitza.c_str());

    adifitx << proiektua_oro << endl;
    adifitx << puntu[objin].mota << endl;
    adifitx << ft.maizkop << endl;
    adifitx << koefkopmax[objin] << endl;
    adifitx << "ZENT  " << fixed << setprecision(12) << "\t" << puntu[objin].zentroa[0] << "\t" << puntu[objin].zentroa[1] << "\t" << puntu[objin].zentroa[2] << endl;
    adifitx << "MAT1  " << showpos << puntu[objin].mat[0] << "\t" << puntu[objin].mat[1] << "\t" << puntu[objin].mat[2]  << "\t" << puntu[objin].mat[3] << endl;

    for (maizkont=0; maizkont < ft.maizkop; maizkont++)
    {
        double f=ft.espektrua[maizkont];
        adifitx << "FREQ  " << maizkont+1 << "\t" << f << "\t" << puntu[objin].shkop[maizkont] << "\t" << puntu[objin].koefkop[maizkont] << "\t" << honderl_oro[maizkont] << endl;
        for (kont=0; kont<puntu[objin].koefkop[maizkont]; kont++)
        {
            adifitx << "KOEF  " << maizkont+1 << " " << kont+1 << scientific << setprecision(20) << "\t" << puntu[objin].koef[maizkont][kont] << endl;
        }
    }
    adifitx << "AMAI  ";
    cout << puntu[objin].izena << ".adi fitxategia idatzita." << endl;
    adifitx.close();
}

} /** AURRE aplikazioaren amaiera **/
