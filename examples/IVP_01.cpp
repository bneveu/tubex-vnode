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

void ivp01(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    yp[0] = -y[0];

}


void contract(TubeVector& x, double t0, bool incremental)
{

    // Differential equation


    int n=1;
    double t=0;
    double tend=10;



    AD *ad=new FADBAD_AD(n,ivp01,ivp01);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}

int main()

{   

  TFunction f("x", "-x");
  float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    
    Interval domain(0.,10.);
    TubeVector x(domain, 1);
    IntervalVector v(1);
    double init = exp(-10);
    v[0]=Interval(init);

    x.set(v, 10.); // ini


    double eps=0.005;

    /* =========== SOLVER =========== */
    Vector epsilon(1, eps);


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.  );
    solver.set_var3b_fxpt_ratio(-1);

   // solver.set_var3b_timept(1);
    solver.set_trace(1);

    
    //    solver.set_stopping_mode(1);
    solver.set_refining_mode(0.);
    solver.set_bisection_timept(0);
    solver.set_contraction_mode(4);

    solver.set_max_slices(40000);
    //    list<TubeVector> l_solutions = solver.solve(x,f);


    // Parameters for vnode
    //    solver.set_max_slices(1);
    

    
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max gate diam : " <<l_solutions.front().max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

   
    return 0;
}
