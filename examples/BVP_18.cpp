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

void bvp18(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    yp[0] = y[1];
    yp[1] = -exp(y[0]);
}


void contract(TubeVector& x, double t0, bool incremental)
{
  //    incremental=false;
    // Differential equation

    int n=2;
    double t=0;
    double tend=1;

    AD *ad=new FADBAD_AD(n,bvp18,bvp18);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}

int main() {
  TFunction f("x1", "x2" ,"(x2;-exp(x1))");
    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution

    /* =========== PARAMETERS =========== */

    Interval domain(0.,1.);
    TubeVector x(domain,2);
    IntervalVector v(2);
    v[0]=Interval(0.,0.);
    //    v[1]=Interval(-1.e8,1.e8);
    v[1]=Interval(-20.,20.);
    x.set(v, 0.); // ini
    v[0]=Interval(0.,0.);
    v[1]=Interval(-20.,20.);
    x.set(v,1.);
    
    double eps0=5.e-2;
    double eps1=5.e-2;
    /*
    double eps0=0.005;
    double eps1=0.005;
    */
    /* =========== SOLVER =========== */
    Vector epsilon(2);
    epsilon[0]=eps0;
    epsilon[1]=eps1;

    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.9);
    // solver.set_var3b_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_propa_fxpt_ratio(0.999);
    solver.set_var3b_external_contraction(true);
   // solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(2000);
    //solver.set_max_slices(1);
    solver.set_bisection_timept(3);
    solver.set_refining_mode(0);
    solver.set_stopping_mode(0);
    solver.set_contraction_mode(4);
    solver.set_var3b_external_contraction(true);
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
       cout << "nb sol " << l_solutions.size() << endl;
    if(l_solutions.size()>0){
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<< ")" << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;
    cout << l_solutions.back()<<" ti-> " <<l_solutions.back()(domain.lb()) << " tf -> "<< l_solutions.back()(domain.ub()) <<" max diam : (" <<l_solutions.back()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.back()[1].max_gate_diam(t_max_diam)<< ")"<< " volume :  "<< l_solutions.back().volume()<<" ti (diam) -> " <<l_solutions.back()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.back()(domain.ub()).diam() << endl;
    }

    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;
}
/*

solver.set_cid_timept(0);
solver.set_bisection_timept(-1);
solver.set_trace(1);
solver.set_max_slices(1);
solver.set_refining_mode(2);

  */
//
//TubeVector (dim 2) [0, 1]↦([-0.53469, 0.295158] ; [-1.20621, 0.77602]), 153 slices ti-> ([0, 0] ; [-0.0809922, 0.77602]) tf -> ([0, 0] ; [-1.08613, -0.319046]) max diam : (0.829847 ; 0.959983) volume :  1.44828 ti (diam) -> (0 ; 0.857012) tf (diam) -> (0 ; 0.767083)
//TubeVector (dim 2) [0, 1]↦([-0.415218, 4.12288] ; [-11.0195, 10.9699]), 173 slices ti-> ([0, 0] ; [10.82, 10.9699]) tf -> ([0, 0] ; [-10.982, -10.8134]) max diam : (0.569361 ; 1.10763) volume :  0.720347 ti (diam) -> (0 ; 0.149883) tf (diam) -> (0 ; 0.168621)

//TubeVector (dim 2) [0, 1]↦([-1.20465e-13, 0.179116] ; [-0.723309, 0.604904]), 42 slices ti-> ([0, 0] ; [0.44696, 0.604904]) tf -> ([0, 0] ; [-0.602787, -0.519447]) max diam : (0.0956414 ; 0.255698) volume :  0.317889 ti (diam) -> (0 ; 0.157943) tf (diam) -> (0 ; 0.0833401)
//TubeVector (dim 2) [0, 1]↦([-0.214989, 4.09685] ; [-10.8724, 10.8489]), 121 slices ti-> ([0, 0] ; [10.8466, 10.8489]) tf -> ([0, 0] ; [-10.849, -10.8466]) max diam : (0.423897 ; 0.999834) volume :  0.537764 ti (diam) -> (0 ; 0.00224738) tf (diam) -> (0 ; 0.00237805)
