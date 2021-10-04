#include <iomanip>
#include <sstream>
#include<vector>
#include <math.h>
const double PI = 3.14159265358979323846;

using namespace std;

void wsquare(vector<double> xnodo, double omega, double d, double h, double rho, double a, double b, double &wsq)
{

unsigned int m, n;

wsq=0.0;

for (m=1; m<101; m+=2)
{
    for (n=1; n<101; n+=2)
    {
        wsq=wsq+(sin(m*PI*xnodo[0]/a)*sin(n*PI*xnodo[1]/b))/(m*n*(pow(PI,4)*d*pow(pow(m/a,2)+pow(n/b,2),2)-rho*h*pow(omega,2)));
    }
}

wsq=wsq*16*1/pow(PI,2);

}

