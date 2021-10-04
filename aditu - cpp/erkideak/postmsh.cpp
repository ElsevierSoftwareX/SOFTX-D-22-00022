#include <iomanip>
#include <sstream>
#include <math.h>
#include <complex>
#include "../erkideak/ObjektuPuntu.h"
#include "../erkideak/postmsh.h"

using namespace std;

const double PI = 3.14159265358979323846;
const complex<double> I{0.0, 1.0};

template <typename T, typename U>
inline complex<T> operator*(complex<T> lhs, const U& rhs)
{
    return lhs *= rhs;
}

void postmsh(string punizena)
{

unsigned int kont;

string bidea= "C:/Users/mapgazug/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
//    string bidea= "C:/Users/aita/Documents/proiektuak/C++/datuak/"; /** proiektuko fitxategien direktorioa **/
string lerroa; /** Irakurritako lerro osoa **/

string txartela="HASI";
string pun, msh;

unsigned int etiketa, nodele; /** Irakurtzeko laguntzaileak **/
unsigned int osoa; /** Irakurtzeko int **/
double bik; /**Irakurtzeko double **/
vector<double> dtemp;
vector <unsigned int> itemp;

vector<vector<double>> puntuak;

vector<vector<unsigned int>> elementuak; /** Lerroak txertatzeko bitarteko bektorea **/

pun = bidea + punizena + ".pun";
msh = bidea + punizena + ".msh.post";

ifstream punfitx(pun.c_str(), ios::in);
ofstream mshfitx(msh.c_str());

getline(punfitx,lerroa); /** ObjektuAurre-mota **/

getline(punfitx,lerroa); /** Nodo kopurua **/
istringstream ls2(lerroa);
ls2 >> osoa;

mshfitx << "MESH " << punizena << "dimension " << 3 << " Elemtype " << "Quadrilateral " << "Nnode " << osoa << endl;

while(txartela.compare("AMAI") != 0)
{
getline(punfitx,lerroa);
istringstream ls(lerroa);
ls >> txartela;

if (txartela.compare("GRID")==0)
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
    elementuak.push_back(itemp);
    itemp.clear();
}

else if (txartela.compare("AMAI")==0)
{

}

} /** of while txartela **/

mshfitx << "coordinates" << endl;
for (kont=0; kont < puntuak.size(); kont++)
{
    mshfitx << kont+1 << " " << fixed << setprecision(10) << puntuak[kont][1] << " "  << puntuak[kont][2] << " " << puntuak[kont][0] << endl;
}
mshfitx << "end coordinates" << endl;

mshfitx << "elements" << endl;
for (kont=0; kont < elementuak.size(); kont++)
{
    mshfitx << kont+1 << " " << elementuak[kont][0] << " " << elementuak[kont][1] << " " << elementuak[kont][2] << " " << elementuak[kont][3] << endl;
}
mshfitx << "end elements" << endl;



}

