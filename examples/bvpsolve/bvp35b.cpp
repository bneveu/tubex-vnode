//created by neveu sept 30 2020
// problem Bvpsolve35 


#include <iostream>

#include "tubex.h"
#include "CtcVnodelp.h"
#include "tubex-solve.h"

using namespace std;
using namespace ibex;
using namespace tubex;
using namespace vnodelp;
template<typename var_type>

void bvp35(int n, var_type*yp, const var_type*y, var_type t, void*param)
{

  interval ksi = 10;
  
  yp[0] = 1;
  yp[1] = y[2];
  yp[2] = ksi*(y[0]*y[2]-y[1]);


  
}
AD *ad=new FADBAD_AD(3,bvp35,bvp35);

void contract(TubeVector& x, double t0, bool incremental)
{

    int n=3;
    double t=-1;
    double tend=1;
    CtcVnodelp c;
   
    if (x.volume() < DBL_MAX && x.nb_slices()>=2  ) {   
      c.preserve_slicing(true);
      c.set_ignoreslicing(true);
    }
    else {
      c.preserve_slicing(false);
      c.set_ignoreslicing(true);
    }
     
    c.set_vnode_hmin(1.e-3);
    c.set_vnode_order(11);
    cout << " before vnode " << t0 << " " << incremental << endl;
    c.Contract(ad,t,tend,n,x,t0,incremental);
    cout << " after vnode " <<   x.nb_slices() << endl;
}

int main() {
  TFunction f("x1", "x2" ,"x3","(1;x3;10*(x1*x3-x2))");

   
    /* =========== PARAMETERS =========== */

    Interval domain(-1.,1.);
   
    TubeVector x (domain, 3);
   
    IntervalVector v(3);
    v[0]=Interval(-1.,-1.);
    v[1]=Interval(1);
    v[2]=Interval(-100,100);
    
    x.set(v, -1.); // ini
    v[0]=Interval(-1.,1.);
    v[1]=Interval(2);
    v[2]=Interval(-100,100);
    
    x.set(v,1.);
    
        
    double eps0=0.1; 
    double eps1=0.1;
    double eps2=0.1;
    
    
    /* =========== SOLVER =========== */
    Vector epsilon(3);
    epsilon[0]=eps0;
    epsilon[1]=eps1;
    epsilon[2]=eps2;


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2.0);
    //solver.set_propa_fxpt_ratio(0.999);
    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_fxpt_ratio(-1);
    //solver.set_var3b_fxpt_ratio(0.9);

    solver.set_var3b_propa_fxpt_ratio(0.9);

    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(5000);
    //    solver.set_max_slices(1);
    solver.set_bisection_timept(3);

    solver.set_refining_mode(2);
    solver.set_stopping_mode(0);
    solver.set_contraction_mode(4);
    solver.set_var3b_external_contraction(true);

    std::ofstream Out("err.txt");
    std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());
    list<TubeVector> l_solutions = solver.solve(x, f, &contract);
    //    list<TubeVector> l_solutions = solver.solve(x, &contract);
    //list<TubeVector> l_solutions = solver.solve(x, f);
    std::cerr.rdbuf(OldBuf);
    
    cout << "nb sol " << l_solutions.size() << endl;
    if(l_solutions.size()>0){
    double t_max_diam;

    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max gate diam : (" <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam)<<", "<< l_solutions.front()[2].max_gate_diam(t_max_diam)<< " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

    }

    return 0;
}
