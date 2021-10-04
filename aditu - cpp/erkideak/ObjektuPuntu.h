#ifndef ObjektuPuntu_H_INCLUDED
#define ObjektuPuntu_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include<string>
#include<complex>

using namespace std;

class ObjektuPuntu
{
public:
    unsigned int ordinala; /** Ordinala proiektu osoko objektuen zerrendan **/
    unsigned int kok1, kok2;
    string izena;
    string mota; /** OSKO : OSKOl elastikoa / FFIN :  Fluido FINitua / FINF : Fluido INFinitua **/
    unsigned int puntukop; /** Nodo kopurua **/
    vector<double> zentroa; /** Erreferentzia-sistemaren jatorria **/
    vector<double> mat; /** Materialaren propietateak **/

    unsigned int luzera; /** Kopefiziente kopurua maiztasuneko **/

    vector<unsigned int> mbzenb; /** Nodoan ezarritako MB zenbakia (1,2,...) **/
    vector<unsigned int> mbmota; /** Nodoan ezarritako MB mota (0:desakoplatua / 1:akoplatua) **/
    vector<unsigned int> mbeku; /** Nodoan ezarritako MBak ezartzen dituen ekuazio kopurua **/
    vector<vector<int> > mbobj; /** Nodoan ezarritako MBak ezartzen duen taulako kokapena/objektua **/

    vector<vector<double> > puntuak;
    vector<vector<double> > normalak;

    vector<unsigned int> akopla;  /** Objektuaren akoplamanduak gainerako objektuekin **/

    int harmonikomax;

    vector<unsigned int > shkop;
    vector<unsigned int > koefkop;
    vector<vector<complex<double> > > koef;

    /** Eraikitzaileak **/

    ObjektuPuntu();  /** objektu hutsa sotrtzeko eraikitzailea (binarioa iraultzeko) **/

    ObjektuPuntu(string, unsigned int, unsigned int, int, int);  /** char: Neuk sortutako .mub testu-fitxategia; unsigned int : ordinala;  unsigned int : objkop ; unsigned int : kok1 **/

    ObjektuPuntu(string, unsigned int, unsigned int, int);  /** char: GiDek sortutako .mub testu-fitxategia; unsigned int : ordinala;  unsigned int : objkop ; unsigned int : kok1 **/

    /** Desegilea **/

    ~ObjektuPuntu();

};

#endif /** ObjektuPuntu_H_INCLUDED **/
