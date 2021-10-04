#include <iostream>
#include "doiketa.h"
#include "../erkideak/ObjektuPuntu.h"
#include <complex>
#include <vector>
const complex<double> I{0.0, 1.0};

using namespace std;

void doiketa(vector<ObjektuPuntu> &puntu, vector<vector<unsigned int> > &ekkop, unsigned int objkop)
{
    unsigned int obi, obj, konta, puntua;
    vector<vector<unsigned int> > dirichlet;
    dirichlet.clear();
    dirichlet.resize(objkop, vector<unsigned int>(objkop, 0));

    /** akoplamenduak murriztu **/

    for (obi=0; obi<objkop; obi++)
    {
        puntu[obi].akopla.clear();
        puntu[obi].akopla.resize(objkop, 0);
        for (puntua=0; puntua<puntu[obi].puntukop; puntua++)
        {
            switch (puntu[obi].mbmota[puntua])
            {
                case 0: /** Desakoplatua **/
                    puntu[obi].akopla[obi]+=puntu[obi].mbeku[puntua];
                    ekkop[obi][obi]+=puntu[obi].mbeku[puntua];
                break;

                case 1: /** Akoplatua Neumann **/
                    obj=puntu[obi].mbobj[puntua][0];
                    puntu[obi].akopla[obj]+=puntu[obi].mbeku[puntua];
                    ekkop[obj][obi]+=puntu[obi].mbeku[puntua];
                break;

                case 2: /** Akoplatua Dirichlet **/
                    obj=puntu[obi].mbobj[puntua][0];
                    dirichlet[obi][obj]+=puntu[obi].mbeku[puntua];
                break;
            } /** end of switch **/
        }
    }

    for (obi=0; obi<objkop; obi++)
    {
        for (obj=obi+1; obj<objkop; obj++)
        {
            if (dirichlet[obi][obj]>0)
            {
            if (dirichlet[obi][obj] >= dirichlet[obj][obi])
            {
                puntu[obi].akopla[obj]+=floor(dirichlet[obi][obj]/2);
                puntu[obj].akopla[obi]+=floor(dirichlet[obi][obj]/2);

                ekkop[obi][obj]+=2*floor(dirichlet[obi][obj]/2);
            }
            else
            {
                puntu[obi].akopla[obj]+=floor(dirichlet[obj][obi]/2);
                puntu[obj].akopla[obi]+=floor(dirichlet[obj][obi]/2);

                ekkop[obj][obi]+=2*floor(dirichlet[obj][obi]/2);
            }
            }
        }
    }

    /** kok1 eta kok2 kokapenak zehaztu **/

    konta=0;

    for (obi=0; obi<objkop; obi++)
    {
        if (obi==0)
        {
            puntu[obi].kok1=0;
        }
        else
        {
            puntu[obi].kok1=konta;
        }

        for (obj=0; obj<objkop; obj++)
        {
            konta=konta+puntu[obi].akopla[obj];
        }
        puntu[obi].kok2=konta-1;
    }

} /** end of doiketa **/
