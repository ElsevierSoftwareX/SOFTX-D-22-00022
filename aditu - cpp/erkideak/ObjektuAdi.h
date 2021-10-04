#ifndef ObjektuAdi_H_INCLUDED
#define ObjektuAdi_H_INCLUDED
#include <iostream>
#include <fstream>
#include <vector>
#include<string>
#include<complex>

using namespace std;

class ObjektuAdi
{
public:
    unsigned int ordinala; /** Ordinala proiektu osoko objektuen zerrendan **/
    string proiektua; /** Proiektuaren izena **/
    vector<double> espektrua; /** Proiektuko maiztasun-espektrua **/
    vector<float> honderl; /** maiztasun bakoitzean egindako errore erlatiboa (%) **/
    string izena;
    string mota; /** ELAS : Egitura eLAStikoa / FFIN :  Fluido FINitua / FINF : Fluido INFinitua **/
    unsigned int nodokop; /** Nodo kopurua **/
    unsigned int maizkop; /** Maiztasun kopurua **/
    vector<double> zentroa; /** Erreferentzia-sistemaren jatorria **/
    vector<double> mat; /** Materialaren propietateak **/
    vector<vector<double> > nodoak;
    vector<vector<double> > normalak;

    vector<unsigned int> shkop;   /** Harmoniko eskeriko kopurua maiztasuneko **/
    vector<unsigned int> koefkop; /** Goordeko diren koefiziente kopurua maiztasuneko **/
    vector<vector<complex<double> > > koef;

    /** Eraikitzaileak **/

    ObjektuAdi(string);                /** adi fitxategi batetik eraikitzeko **/

    ObjektuAdi();                /** objektu hutsa sotrtzeko eraikitzailea (binarioa iraultzeko) **/

    /** Desegilea **/

    ~ObjektuAdi();

};

#endif /** ObjektuAdi_H_INCLUDED **/
