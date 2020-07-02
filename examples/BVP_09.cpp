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
{   interval a=interval(0.7); interval b=interval(log(2.)/5.);
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

{   
  tubex::Function f("y1", "y2", "(-0.7*y1 ; 0.7*y1 - (ln(2)/5.)*y2)");
  float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    
    cout << log(2./5.)<<endl;
    int n = 2;
    Interval domain(0.,6.);
    TubeVector x(domain,n,IntervalVector(2,Interval(-1e5,1e5)));


    // Boundary condition:
    IntervalVector init = x(x.domain().lb());
    init[0] = Interval(1.25);
    // init[1] = Interval(1.1,1.3);
    //init[1]=Interval(-0.413732, -0.372333);
    x.set(init, 0.);

    // Additional restriction (maximum value):
    Interval domain_restriction(1.,3.);
    IntervalVector max_restriction(2);
    max_restriction[0] = Interval(-1e5,1e5);
    max_restriction[1] = Interval(1.1,1.3);

//     domain_restriction(1.);
//    IntervalVector max_restriction(2);
//    max_restriction[0] = Interval(-9999.,9999.);
//    max_restriction[1] = Interval(1.1,1.3);
    x.set(max_restriction,domain_restriction);
    /* =========== SOLVER =========== */

//    contract(x);

    Vector epsilon(n);
    epsilon[0] = 0.04;
    epsilon[1] = 0.04;
//
//
    tubex::Solver solver(epsilon);
//
    solver.set_refining_fxpt_ratio(2);
    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.99);

//
  //  solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(5000);
    //    solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(-2);
    solver.set_contraction_mode(2);
    solver.set_stopping_mode(0);
//
    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x,&contract);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<< ")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
    
  
}
