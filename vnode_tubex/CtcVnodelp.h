#ifndef __CTCVNODELP_H__
#define __CTCVNODELP_H__

#include "ibex.h"
#include "tubex.h"
#include "tubex_DynCtc.h"
#include "tubex_Slice.h"
#include "vnode.h"
#include <vector>
#include "vnode.h"


using namespace ibex;

const int ord (11);
namespace tubex {
    typedef std::vector<std::pair<double, ibex::IntervalVector>> Vstate;

    class CtcVnodelp : public DynCtc {
    public:

        CtcVnodelp();

        void fill_state_vector(Tube &x,vector<double> &time_gate ,vector<ibex::IntervalVector>& si, const int& i);
        //ode type: 0 for ivp 1 for bvp
        void Contract(vnodelp::AD *ad,double t,double tend,int n, Tube &x, double t0, bool incremental);

        void Contract(vnodelp::AD *ad,double t,double tend,int n, TubeVector &x,double t0, bool incremental);


    };
}


#endif
