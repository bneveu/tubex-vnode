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

    // Differential equation


    int n=2;
    double t=0;
    double tend=1;



    AD *ad=new FADBAD_AD(n,bvp18,bvp18);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}
void contract_ctcderiv(TubeVector& x)

/* example 34  bvpsolve  with ksi=1
   2 tubes solutions
 */
{
    tubex::Function f("x1", "x2" ,"(x2;-exp(x1))");

    CtcPicard ctc_picard;
    ctc_picard.preserve_slicing(false);
    if (x.volume() > 50000.)
        ctc_picard.contract(f, x, FORWARD | BACKWARD);

    TubeVector v = f.eval_vector(x);

    CtcDeriv ctc_deriv;
    ctc_deriv.set_fast_mode(true);
    ctc_deriv.contract(x, v, FORWARD | BACKWARD);

    /*
    v=f.eval_vector(x);


    CtcDynCid* ctc_dyncid = new CtcDynCid(f);

    //CtcDynCidGuess*  ctc_dyncid = new CtcDynCidGuess(f);
    ctc_dyncid->set_fast_mode(true);
    CtcIntegration ctc_integration(f,ctc_dyncid);
    ctc_integration.contract(x,v,x[0].domain().lb(),FORWARD) ;
    ctc_integration.contract(x,v,x[0].domain().ub(),BACKWARD) ;
    delete ctc_dyncid;
    */
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
    solver.set_max_slices(1);
    solver.set_bisection_timept(-1);
    //solver.set_refining_mode(0);

    // solver.set_contraction_mode(2);
    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << "nb sol " << l_solutions.size() << endl;
    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_diam() << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;
    cout << l_solutions.back()<<" ti-> " <<l_solutions.back()(domain.lb()) << " tf -> "<< l_solutions.back()(domain.ub()) <<" max diam : " <<l_solutions.back().max_diam() << " volume :  "<< l_solutions.back().volume()<<" ti (diam) -> " <<l_solutions.back()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.back()(domain.ub()).diam() << endl;


    tubex::Function f("x1", "x2" ,"(x2;-exp(x1))");
    TubeVector tvy(domain, 2);
    tvy= l_solutions.front();

     eps0=0.05;
     eps1=0.05;
    /* =========== SOLVER =========== */
    epsilon[0]=eps0;
    epsilon[1]=eps1;


    tubex::Solver solver2(epsilon);

    //    solver.set_refining_fxpt_ratio(0.99999);
    solver2.set_refining_fxpt_ratio(2.0);
    //      solver.set_refining_fxpt_ratio(0.99);

    solver2.set_propa_fxpt_ratio(0.99999);

    solver2.set_var3b_fxpt_ratio(0.9999);
    //solver.set_var3b_fxpt_ratio(0.);
    solver2.set_var3b_timept(0);
    solver2.set_bisection_timept(3);
    solver2.set_trace(1);
    solver2.set_max_slices(20000);
    solver2.set_refining_mode(2);
    solver2.set_contraction_mode(2);
    list<TubeVector> l_solutions1 = solver2.solve(tvy, &contract_ctcderiv);
    cout << "nb sol " << l_solutions1.size() << endl;

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
