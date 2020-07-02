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


void contract(TubeVector& x, double& t0, bool incremental)
{
    //ibex-tubex eq diff
    tubex::Function f("x1", "x2" ,"(x2;-exp(x1))");

    // vnode Differential equation


    int n=2;
    double t=0;
    double tend=1;



    AD *ad=new FADBAD_AD(n,bvp18,bvp18);
    CtcVnodelp c;
    if (x.volume()>0.64)
        c.Contract(ad,t,tend,n,x,t0,incremental);
    else{
        TubeVector v = f.eval_vector(x);

        CtcDeriv ctc_deriv;
        ctc_deriv.set_fast_mode(true);
        ctc_deriv.contract(x, v, FORWARD | BACKWARD);
    }
}

int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution

    /* =========== PARAMETERS =========== */

    Interval domain(0.,1.);
    TubeVector x(domain, 2);
    IntervalVector v(2);
    v[0]=Interval(0.,0.);
    //    v[1]=Interval(-1.e8,1.e8);
    v[1]=Interval(-20.,20.);
    x.set(v, 0.); // ini
    v[0]=Interval(0.,0.);
    v[1]=Interval(-20.,20.);
    x.set(v,1.);

    double eps0=1;
    double eps1=1;
    /* =========== SOLVER =========== */
    Vector epsilon(2);
    epsilon[0]=eps0;
    epsilon[1]=eps1;


    tubex::Solver solver(epsilon);

   // solver.set_refining_fxpt_ratio(0.1);
    solver.set_propa_fxpt_ratio(0.);
    solver.set_var3b_fxpt_ratio(0.9);

    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(400000);
    solver.set_bisection_timept(-1);
    //solver.set_refining_mode(0);

     //solver.set_contraction_mode(2);
    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << "nb sol " << l_solutions.size() << endl;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_diam() << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;
    cout << l_solutions.back()<<" ti-> " <<l_solutions.back()(domain.lb()) << " tf -> "<< l_solutions.back()(domain.ub()) <<" max diam : " <<l_solutions.back().max_diam() << " volume :  "<< l_solutions.back().volume()<<" ti (diam) -> " <<l_solutions.back()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.back()(domain.ub()).diam() << endl;


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
