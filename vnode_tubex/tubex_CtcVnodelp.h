#ifndef __TUBEX_CTCVNODELP_H__
#define __TUBEX_CTCVNODELP_H__

#include "ibex.h"
#include "tubex.h"
#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include "vnode.h"
#include <vector>
#include "vnode.h"
const int ord (11);
namespace tubex {
    typedef std::vector<std::pair<double, ibex::IntervalVector>> Vstate;

    class CtcVnodelp : public Ctc {
    public:

        CtcVnodelp();

        void set_state(Vstate &s);

        void show_state_info();

        void OdeContractor(vnodelp::AD *ad, double t, double tend,int n, Tube &x, Vstate state);

        void OdeContractor(vnodelp::AD *ad, double t, double tend,int n, TubeVector &x, Vstate state);

        void Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x);

        void fill_state_vector(Tube &x,vector<double> &time_gate, vector<ibex::IntervalVector>& si,const int& i);

    private:
        Vstate state;


    };
}

#endif
