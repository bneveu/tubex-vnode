//created by bneveu : 3rd variant of Bvp solve22

#include <iostream>
#include "tubex.h"
#include "CtcVnodelp.h"
#include "tubex-solve.h"

using namespace std;
using namespace ibex;
using namespace tubex;
using namespace vnodelp;
template<typename var_type>

void bvp22c(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    interval ksi=100;

    yp[0] = y[1];
    yp[1] = -(ksi*(y[1]+sqr(y[0])));
}

AD *ad=new FADBAD_AD(2,bvp22c,bvp22c);

void contract(TubeVector& x, double t0, bool incremental)
{

    // Differential equation

    int n=2;
    double t=0;
    double tend=1;


    CtcVnodelp c;
   
    if (x.volume() < DBL_MAX && x.nb_slices() >=2)
      {c.preserve_slicing(true);
      c.set_ignoreslicing(true);
      }
    else {c.preserve_slicing(false);
       c.set_ignoreslicing(true);
    }
    //    c.set_vnode_hmin(5.e-4);
    c.set_vnode_hmin(1.e-3);
    c.Contract(ad,t,tend,n,x,t0,incremental);


}

int main()

{    float temps;

    TFunction f("x1", "x2" ,"(x2;-100*(x2+x1^2))");
    Interval domain(0.,1.);

    TubeVector x(domain,2);

    IntervalVector v(2);
    v[0]=Interval(0.,0.);
    v[1]=Interval(-100,100);

    x.set(v, 0.); // ini
    v[0]=Interval(0.5,0.5);
    v[1]=Interval(-100,100);
    x.set(v,1.);

    double eps0=0.2;
    double eps1=0.2;

    /* =========== SOLVER =========== */
    Vector epsilon(2);
    epsilon[0]=eps0;
    epsilon[1]=eps1;


    tubex::Solver solver(epsilon);

    solver.set_refining_fxpt_ratio(2);

    solver.set_propa_fxpt_ratio(0.);

    //solver.set_var3b_fxpt_ratio(-1);
    solver.set_var3b_fxpt_ratio(0.9);
    solver.set_var3b_propa_fxpt_ratio(0.9);
    solver.set_var3b_timept(3);
    solver.set_trace(1);
    solver.set_max_slices(100000);
    //    solver.set_max_slices(1);
    solver.set_refining_mode(2);
    solver.set_bisection_timept(3);
    solver.set_contraction_mode(2);
    solver.set_var3b_external_contraction(true);
    solver.set_stopping_mode(0);;
    std::ofstream Out("err.txt");
    std::streambuf* OldBuf = std::cerr.rdbuf(Out.rdbuf());

    list<TubeVector> l_solutions = solver.solve(x,f, &contract);
    //list<TubeVector> l_solutions = solver.solve(x, &contract);
    //list<TubeVector> l_solutions = solver.solve(x,f);
    std::cerr.rdbuf(OldBuf);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    double t_max_diam;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front()[0].max_gate_diam(t_max_diam)<<", "<<l_solutions.front()[1].max_gate_diam(t_max_diam) << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;

    return 0;
}
