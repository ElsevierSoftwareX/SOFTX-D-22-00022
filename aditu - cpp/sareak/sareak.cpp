#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include<vector>
#include<math.h>
#include<random>

const double PI = 3.14159265358979323846;

using namespace std;

int main()
{
    string hitza;

    vector<vector<double> >  koord; /** Koordenatuak **/
    vector<vector<double> >  normalak; /** Normalak **/
    vector<vector<unsigned int> >  mubal; /** Koordenatuak **/

    koord.clear();
    normalak.clear();
    mubal.clear();

    vector<double> xpuntu(3, 0.0); /** puntuko koordenatuak sartzeko **/
    double mod;
    vector<double> npuntu(3, 0.0); /** puntuko normalaren osagaiak sartzeko **/
    vector<unsigned int> mbpuntu(4, 0.0); /** puntuko MBak sartzeko **/

    int zenbat, k;
    double erradioa, faktorea;
    double zentroa[3];
    int muba[4];

    unsigned int puntukop, n;

    unsigned int u, v, puntukont=0;
    double theta, fi, urraskopu, urraskopv;

    cout << "Fitxategiaren izena? " << endl;
    cin >> hitza;

    zentroa[0]=0.0;
    zentroa[1]=0.0;
    zentroa[2]=0.0;

    cout << "Zenbat esfera?" << endl;
    cin >> zenbat;

    for (k=0; k<zenbat; k++)
        {
        cout << k+1 << ". esfera. Puntu kopurua? " << endl;
        cin >> puntukop;

        cout << "Erradioa? " << endl;
        cin >> erradioa;

        cout << "Normalerako faktorea: (-1/+1) " << endl;
        cin >> faktorea;

        cout << k+1 << ". esferako MB zenbakia:" << endl;
        cin >> muba[0];

        cout << k+1 << ". esferako MB mota (0/1):" << endl;
        cin >> muba[1];

        cout << k+1 << ". esferako MBko am kopurua:" << endl;
        cin >> muba[2];

        cout << k+1 << ". esferako MBren taulako kokapena / objektua (1/obj):" << endl;
        cin >> muba[3];

        n=ceil(sqrt(puntukop/2));
        urraskopu=PI/n;
        urraskopv=(2*PI)/(2*n);

        for (u=0; u < n; u++)
        {
            theta=(0.0+urraskopu/2)+u*urraskopu;

            for (v=0; v < 2*n; v++)
            {
                fi=(0.0+urraskopv/2)+v*urraskopv;

                puntukont++;

                xpuntu[0]=zentroa[0]+erradioa*sin(theta)*cos(fi);
                xpuntu[1]=zentroa[1]+erradioa*sin(theta)*sin(fi);
                xpuntu[2]=zentroa[2]+erradioa*cos(theta);

                mod=sqrt(pow(xpuntu[0]-zentroa[0],2)+pow(xpuntu[1]-zentroa[1],2)+pow(xpuntu[2]-zentroa[2],2));
                npuntu[0]=faktorea*xpuntu[0]/mod;
                npuntu[1]=faktorea*xpuntu[1]/mod;
                npuntu[2]=faktorea*xpuntu[2]/mod;

                /** Mugalde baldintza: Neumann**/
                mbpuntu[0]=muba[0];
                mbpuntu[1]=muba[1];
                mbpuntu[2]=muba[2];
                mbpuntu[3]=muba[3];

                koord.push_back(xpuntu);
                normalak.push_back(npuntu);
                mubal.push_back(mbpuntu);
            }
        }
    }

/**** TESTU FITXATEGIA IDATZI ********/

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /* proiektuko fitxategien direktorioa */
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /* proiektuko fitxategien direktorioa */

hitza = bidea + hitza;
ofstream fitx(hitza.c_str());

fitx << "FINF" << endl;
fitx << puntukont << endl;
fitx << "ZENT  " << fixed << setprecision(12) << "\t" << zentroa[0] << "\t" << zentroa[1] << "\t" << zentroa[2] << endl;
fitx << "MAT1  " << showpos << +142355.29 << "\t" << +1.21 << "\t" << 343.0 << "\t" << +415.0 << endl;


for (u=0; u<puntukont; u++)
{
    fitx << "GRID  " << u+1 << fixed << setprecision(12) << "\t" << koord[u][0] << "\t" << koord[u][1] << "\t" << koord[u][2] << endl;
}

for (u=0; u<puntukont; u++)
{
    fitx << "NORM  " << u+1 << fixed << setprecision(12) << "\t" << normalak[u][0] << "\t" << normalak[u][1] << "\t" << normalak[u][2] << endl;
}

for (u=0; u<puntukont; u++)
{
    fitx << "MBAL  " << u+1 << fixed << "\t" << mubal[u][0] << "\t" << mubal[u][1] << "\t" << mubal[u][2] << "\t" << mubal[u][3]  << endl;
}

fitx << "AMAI" << endl;

cout << hitza << "Fitxategia idatzita." << endl;
fitx.close();

} /** SAREAK aplikazioaren amaiera **/
