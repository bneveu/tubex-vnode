//created by bedouhene 09/12/2019

#include <iostream>
#include <vector>
#include "tubex.h"
#include "CtcVnodelp.h"
#include "tubex-solve.h"

using namespace std;
using namespace ibex;
using namespace tubex;
using namespace vnodelp;
template<typename var_type>

void bvp09(int n, var_type*yp, const var_type*y, var_type t, void*param)
{   interval a=interval(0.7); interval b=interval(log(2./5.));
    yp[0] = -a*y[0];
    yp[1] = a*y[0] -b*y[1] ;
}


void contract(TubeVector& x, double& t0, bool incremental)
{
    // Differential equation
    int n=2;
    double t=0;
    double tend=6;

    AD *ad=new FADBAD_AD(n,bvp09,bvp09);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}

int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    
    cout << log(2./5.)<<endl;
    int n = 2;
    Interval domain(0.,6.);
    TubeVector x(domain,n,IntervalVector(n,Interval(-9999,9999)));


    // Boundary condition:
    IntervalVector init = x(x.domain().lb());
    init[0] = 1.25;
    init[1] = Interval(-9999.,9999.);
    x.set(init, x.domain().lb());

    // Additional restriction (maximum value):
    Interval domain_restriction(1.,3.);
    IntervalVector max_restriction(2);
    max_restriction[0] = Interval(-9999.,9999.);
    max_restriction[1] = Interval(1.1,1.3);

    /* =========== SOLVER =========== */

    Vector epsilon(n);
    epsilon[0] = 5.;
    epsilon[1] = 5.;


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(0);
    solver.set_propa_fxpt_ratio(0.99);
    solver.set_var3b_fxpt_ratio(0.9);

    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(1);
    //solver.set_refining_mode(3);
    solver.set_bisection_timept(0);
    // solver.set_contraction_mode(2);

    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;

    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
    
  
}
