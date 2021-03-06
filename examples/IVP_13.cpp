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

void bvp13(int n, var_type*yp, const var_type*y, var_type t, void*param)
{interval a=interval(3.); interval b=interval(2.);
    yp[0] =-y[0]-b*y[1];
    yp[1] = -a*y[0]-b*y[1];
}


void contract(TubeVector& x, double& t0, bool incremental)
{

    // Differential equation


    int n=2;
    double t=0;
    double tend=1;



    AD *ad=new FADBAD_AD(n,bvp13,bvp13);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}

int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    tubex::Function f("x1", "x2" ,"(-x1-2*x2;-3*x1-2*x2)");
    Interval domain(0.,1.);
    TubeVector x(domain,2);
    IntervalVector v(2);
    v[0]=Interval(5.9,6.1);
    v[1]=Interval(3.9,4.1);
    x.set(v, 0.); // ini

    double eps=1.0;

    /* =========== SOLVER =========== */
    Vector epsilon(2, eps);

    tubex::Solver solver(epsilon);

   solver.set_refining_fxpt_ratio(2);
   solver.set_propa_fxpt_ratio(0.999);
   solver.set_var3b_fxpt_ratio(-1);

  //  solver.set_var3b_timept(0);
    solver.set_trace(1);
    // solver.set_max_slices(2000);
    solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(-1);
    solver.set_contraction_mode(2);
    solver.set_stopping_mode(1);
    // list<TubeVector> l_solutions = solver.solve(x, f);
    // list<TubeVector> l_solutions = solver.solve(x, f,&contract);
   list<TubeVector> l_solutions = solver.solve(x,&contract);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<" , "<<l_solutions.front()[1].max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
    
    return 0;
}
