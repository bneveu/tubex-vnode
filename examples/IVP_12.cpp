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

void ivp12(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    yp[0] = -sqr(y[0]);

}


void contract(TubeVector& x, double t0, bool incremental)
{

    // Differential equation


    int n=1;
    double t=0;
    double tend=5;



    AD *ad=new FADBAD_AD(n,ivp12,ivp12);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}

int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    TFunction f("x", "-x^2");
    Interval domain(0.,5.);
    //    TubeVector x(domain,0.1, 1);
    TubeVector x(domain, 1);
    IntervalVector v(1);
    v[0]=Interval(0.1,0.4);

    x.set(v, 0.); // ini

    double epsilon=0.2;
    Vector eps(1,epsilon);

    /* =========== SOLVER =========== */
    tubex::Solver solver(eps);

    solver.set_refining_fxpt_ratio(2);
    solver.set_propa_fxpt_ratio(0.);
    //    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_fxpt_ratio(-1);
    //solver.set_var3b_fxpt_ratio(0.99);
    solver.set_var3b_propa_fxpt_ratio(0.99);

    solver.set_var3b_timept(1);
    solver.set_trace(1);
    solver.set_max_slices(40000);
    //   solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(-2);
    solver.set_contraction_mode(2);
    solver.set_stopping_mode(0);

     list<TubeVector> l_solutions = solver.solve(x, f);
    //     list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << l_solutions.front()(5.) << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;



    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;

}
