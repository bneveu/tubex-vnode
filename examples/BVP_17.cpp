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

void bvp17(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    yp[0] = y[1];
    yp[1] = y[1]/0.2;
}

    AD *ad=new FADBAD_AD(2,bvp17,bvp17);

void contract(TubeVector& x, double t0, bool incremental)
{
  
    // Differential equation

    int n=2;
    double t=0;
    double tend=1;


    //    cout << " vnode contract " << x << endl;

    CtcVnodelp c;
    //c.preserve_slicing(false);
    // c.m_slicevnode=true;
   
    if (x.volume() < DBL_MAX) {c.preserve_slicing(true);
      c.set_ignoreslicing(true);
    }
    else {c.preserve_slicing(false);
       c.set_ignoreslicing(true);
    }
   
    c.set_vnode_hmin(5.e-4);
    c.Contract(ad,t,tend,n,x,t0,incremental);
}

int main()

{    
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    TFunction f("x1", "x2" ,"(x2;x2/0.2)");
    Interval domain(0.,1);
    TubeVector x(domain, 2);
    cout <<  " last slice " << *(x[0].last_slice()) << endl;
    IntervalVector v(2);
     v[0]=Interval(1.,1.);
     // v[1]=Interval(-1e300,1e300);
     v[1]=Interval(-10.,10.);
     //    v[1]=Interval(-1.e8,1.e8);

    x.set(v, 0.); // ini
    v[0]=Interval(0.,0.);
    v[1]=Interval(-10.,10.);
    //    v[1]=Interval(-1.e8,1.e8);
    // v[1]=Interval(-1e300,1e300);
    x.set(v,1.);

    
    double eps0=0.02;
    double eps1=0.02;
    
    /*
    double eps0=1.e-8;
    double eps1=1.e-8;
    */
    /* =========== SOLVER =========== */
    Vector epsilon(2);
	epsilon[0]=eps0;
	epsilon[1]=eps1;


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.);
    solver.set_propa_fxpt_ratio(0.);
    //solver.set_var3b_fxpt_ratio(0.99);
    //    solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    // solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(4000);
    //solver.set_max_slices(1);
    solver.set_refining_mode(0);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(4);
    solver.set_stopping_mode(2);
    solver.set_var3b_external_contraction(true);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    
    
    return 0;
}
