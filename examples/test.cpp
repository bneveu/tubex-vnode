//created by bedouhene 09/12/2019

#include <iostream>
#include <vector>
#include "tubex.h"
#include "CtcVnodelp.h"
#include "tubex-solve.h"
//#include "CtcVnode.h"

using namespace std;
using namespace ibex;
using namespace tubex;
using namespace vnodelp;
template<typename var_type>

void test(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    yp[0] = -y[1];
    yp[1] = y[0];

}

int main(){    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution

    Interval domain(0.,2*M_PI);
    TubeVector x(domain, 2);
    IntervalVector v(2);
    v[0]=Interval(0.).inflate(1e-1);
    v[1]=Interval(0.).inflate(1e-1);
    x.set(v, 0.); // ini

    IntervalVector v1(2);
    v1[0]=Interval(0.).inflate(0);
    v1[1]=Interval(-1.).inflate(0);
    x.set(v1, 2*M_PI);


    int n=2;
    double t=0;
    double tend=2*M_PI;



    AD *ad=new FADBAD_AD(n,test,test);
    CtcVnodelp c;
    double t0=-1;
    bool incremental=false;
    c.Contract(ad,t,tend,n,x,t0,incremental);

    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;

}