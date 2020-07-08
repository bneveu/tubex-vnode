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

void ivp12(int n, var_type*yp, const var_type*y, var_type t, void*param)
{
    yp[0] = -sqr(y[0]);

}


void contract(TubeVector& x, double& t0, bool incremental)
{

    // Differential equation


    int n=1;
    double t=0;
    double tend=5;



    AD *ad=new FADBAD_AD(n,ivp12,ivp12);
    CtcVnodelp c;

    c.Contract(ad,t,tend,n,x,t0,incremental);
}
void contract_ctcderiv(TubeVector& x, double& t0, bool incremental)
{
    TFunction f("x", "-x^2");

    CtcPicard ctc_picard;
    ctc_picard.preserve_slicing(true);

    if (x.volume() > 1.e300)
        ctc_picard.contract(f, x, TimePropag::FORWARD );

    if (x.volume() < 1.e300){
        /*
        CtcDeriv ctc_deriv;
        ctc_deriv.set_fast_mode(true);
        ctc_deriv.contract(x, f.eval_vector(x), FORWARD | BACKWARD);
        */

        TubeVector v = f.eval_vector(x);
        CtcDynCid* ctc_dyncid = new CtcDynCid(f);
        //    CtcDynCidGuess* ctc_dyncid = new CtcDynCidGuess(f);
        ctc_dyncid->set_fast_mode(true);
        CtcIntegration ctc_integration(f,ctc_dyncid);
        ctc_integration.contract(x,v,x[0].tdomain().lb(),TimePropag::FORWARD) ;
        v = f.eval_vector(x);
        ctc_integration.contract(x,v,x[0].tdomain().ub(),TimePropag::BACKWARD) ;
        delete ctc_dyncid;


    }

}
int main()

{    float temps;
    clock_t t1, t2;
    t1=clock();//sert à calculer le temps d'exécution
    
    Interval domain(0.,5.);
    TubeVector x(domain,0.1, 1);
    IntervalVector v(1);
    v[0]=Interval(0.1,0.4);

    x.set(v, 0.); // ini


    double epsilon=2.;
    Vector eps(1,epsilon);

    /* =========== SOLVER =========== */
    tubex::Solver solver(eps);

   // solver.set_refining_fxpt_ratio(0.5);
    solver.set_propa_fxpt_ratio(0.999);
    solver.set_var3b_fxpt_ratio(0.0);

    solver.set_var3b_timept(0);
    solver.set_trace(1);
    solver.set_max_slices(1);
   // solver.set_refining_mode(0.9);
    solver.set_bisection_timept(-1);
    // solver.set_contraction_mode(2);

    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << l_solutions.front()(5.) << endl;
    cout << "nb sol " << l_solutions.size() << endl;

    cout << l_solutions.front()<<" ti-> " <<l_solutions.front()(domain.lb()) << " tf -> "<< l_solutions.front()(domain.ub()) <<" max diam : " <<l_solutions.front().max_diam() << " volume :  "<< l_solutions.front().volume()<<" ti (diam) -> " <<l_solutions.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions.front()(domain.ub()).diam() << endl;


    TubeVector tvy(domain, 1);
    tvy= l_solutions.front();
    epsilon=0.2;
    eps[0]= epsilon;
    tubex::Solver solver2(eps);
    //      solver.set_refining_fxpt_ratio(0.999999);
    //    solver.set_refining_fxpt_ratio(0.99999);
    solver2.set_refining_fxpt_ratio(2.0);
    //    solver.set_refining_fxpt_ratio(0.98);
    //    solver.set_propa_fxpt_ratio(1.);
    //solver.set_propa_fxpt_ratio(0.999);
    //    solver.set_var3b_fxpt_ratio(0.9);
    solver2.set_propa_fxpt_ratio(0);
   // solver2.set_var3b_fxpt_ratio(-1);

    solver2.set_max_slices(40000);

   // solver2.set_var3b_timept(1);
    solver2.set_bisection_timept(-2);

    solver2.set_trace(1);
    solver2.set_refining_mode(0);
//    solver2.set_contraction_mode(2);
//
    list<TubeVector> l_solutions2 = solver2.solve(tvy, &contract_ctcderiv);
    cout <<  "nb sol " << l_solutions2.size() << endl;
    cout << l_solutions2.front()<<" ti-> " <<l_solutions2.front()(domain.lb()) << " tf -> "<< l_solutions2.front()(domain.ub()) <<" max diam : " <<l_solutions2.front().max_diam() << " volume :  "<< l_solutions2.front().volume()<<" ti (diam) -> " <<l_solutions2.front()(domain.lb()).diam() << " tf (diam) -> "<< l_solutions2.front()(domain.ub()).diam() << endl;


    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    cout << "temps ="<< temps << endl<<endl;
    return 0;

}
