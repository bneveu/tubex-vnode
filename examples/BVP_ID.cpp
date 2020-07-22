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

void exp(int n, var_type*yp, const var_type*y, var_type t, void*param)
{   interval a=1,b=2,c=5;
    yp[0] = (a-b*y[0] - c*y[1]);
    yp[1] = (y[0]);
}


void contract(TubeVector& x, double t0, bool incremental)
{
    // Boundary constraints
    Variable vx0, vx1;
    SystemFactory fac;
    fac.add_var(vx0);
    fac.add_var(vx1);
    fac.add_ctr(sqr(vx0) + sqr(vx1) = 1);
    System sys(fac);
    ibex::CtcHC4 hc4(sys);
    IntervalVector bounds(2);
    bounds[0] = x[0](0.);
    bounds[1] = x[0](1.);
    hc4.contract(bounds);
    IntervalVector x0(2);
    IntervalVector x1(2);
    x0[0]=bounds[0];
    x0[1]=0.0;
    x1[0]=bounds[1];
    x1[1]=x[1](1.);

    x.set(x0, 0.);
    x.set(x1, 1.);
    // Differential equation


    int n=2;
    double t=0;
    double tend=1;



    AD *ad=new FADBAD_AD(n,exp,exp);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);

}


int main()
{   float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    
    /* =========== PARAMETERS =========== */

    int n = 2;
    TFunction f("x1", "x2", "(1-2*x1-5*x2;x1)");
    Vector epsilon(n, 0.02);
    Interval domain(0.,1.);
    //    TubeVector x(domain, n, Interval (-1.e100,1.e100));
    TubeVector x(domain, n);
    IntervalVector x0(2);
    IntervalVector x1(2);

    x0[0]=Interval(-10,10);
    x0[1]=Interval(0,0);

    x.set(x0,0.);
    x1[0]=Interval(-10,10);
    x1[1]=Interval(-10,10);
    x.set(x1,1.0);

    TrajectoryVector truth1(domain, TFunction("(exp(-t)*(-(cos(2*t)*(-1 + cos(4) + 2*sin(4) + 4*exp(1)*sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4))))) + sin(2*t)*(2 + 2*cos(4) - sin(4) + 2*exp(1)*(2*exp(1) + sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4)))))))/(5 + 3*cos(4) + 8*exp(2) - 4*sin(4))"));
    TrajectoryVector truth2(domain, TFunction("(exp(-t)*(-(cos(2*t)*(-1 + cos(4) + 2*sin(4) - 4*exp(1)*sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4))))) + sin(2*t)*(2 + 2*cos(4) + 4*exp(2) - sin(4) - 2*exp(1)*sqrt(2*(1 + cos(4) + 2*exp(2) - sin(4))))))/(5 + 3*cos(4) + 8*exp(2) - 4*sin(4))"));


    /* =========== SOLVER =========== */


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2);
    solver.set_propa_fxpt_ratio(0.99);
    //    solver.set_var3b_fxpt_ratio(0.99);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.99);
  //  solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(400);
    //solver.set_max_slices(1);
    solver.set_refining_mode(0);
    //solver.set_stopping_mode(1);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(2);
    solver.set_var3b_external_contraction(false);
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //      list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << l_solutions.front()(domain.ub()).diam() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<< ")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;
    cout << l_solutions.back()<<" ti-> " <<l_solutions.back()(domain.lb()) << " tf -> "<< l_solutions.back()(domain.ub()) <<" max diam : (" <<l_solutions.back()[0].max_gate_diam(t_max_diam)<<")"<< " volume :  "<< l_solutions.back().volume()<<" ti (diam) -> " <<l_solutions.back()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.back()(domain.ub()).diam() << endl;

    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
    
}
