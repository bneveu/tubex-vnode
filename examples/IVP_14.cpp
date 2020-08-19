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

void ivp14(int n, var_type*yp, const var_type*y, var_type t, void*param)
{   
  interval p=interval(-10.);
  yp[0] = p*(y[0]- sin(y[1]))+cos(y[1]);
  yp[1]= interval(1.);
}


void contract(TubeVector& x, double t0, bool incremental)
{

    // Differential equation

    int n=2;
    double t=0;
    double tend=3;

    AD *ad=new FADBAD_AD(n,ivp14,ivp14);
    CtcVnodelp c;
  
    c.Contract(ad,t,tend,n,x,t0,incremental);
  
}

int main()

{    
   
    TFunction f("x1", "x2" ,"(-10*(x1-sin(x2))+cos(x2);1)");
    Interval domain(0.,3.);
    TubeVector x(domain, 2);
    IntervalVector v(2);
    v[0]=Interval(0.);
    v[1]=Interval(0.);
    x.set(v, 0.); // ini


    double eps=0.01;

    /* =========== SOLVER =========== */
    Vector epsilon(2, eps);


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.);
     solver.set_var3b_fxpt_ratio(-1);
    //solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(1000);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(0);
    solver.set_contraction_mode(2);
    solver.set_stopping_mode(0);
    //    cout << " before solve " << x << endl;

    // list<TubeVector> l_solutions = solver.solve(x, &contract);
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x, f);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;

    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

    return 0;
}
