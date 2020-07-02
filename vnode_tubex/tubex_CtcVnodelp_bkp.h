
#ifndef __TUBEX_CTCVNODELP_H__
#define __TUBEX_CTCVNODELP_H__

#define FWD true
#define BWD false

#include "ibex.h"
#include "tubex.h"
#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include "vnode.h"
#include <vector>
#include "vnode.h"
const int ord (7);
namespace tubex {
    typedef std::vector<std::pair<double, ibex::IntervalVector>> Vstate;

    class CtcVnodelp : public Ctc {
    public:

        CtcVnodelp();

        void set_state(Vstate &s);

        void show_state_info();

        void FwdBwdIntegration(vnodelp::AD *ad, double t, double tend, int n,const TubeVector &x,TubeVector &transition_x, Vstate &state,bool direction);

        void OdeContractor(vnodelp::AD *ad, double t, double tend,int n, Tube &x, Vstate state);

        void OdeContractor(vnodelp::AD *ad, double t, double tend,int n, TubeVector &x, Vstate state);

        void GuaranteedIntegration(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, Vstate state);

        void Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x);

        void fill_state_vector( Tube &x,vector<double> &time_gate, vector<ibex::IntervalVector>& si,const int& i);


    private:
        Vstate state;


    };
}

#endif
