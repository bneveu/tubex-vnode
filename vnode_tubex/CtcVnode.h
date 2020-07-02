#ifndef __CTCVNODE_H__
#define __CTCVNODE_H__

#include "ibex.h"
#include "tubex.h"
#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include "vnode.h"
#include <vector>
#include "vnode.h"
const int ord (9);
namespace tubex {
    typedef std::vector<std::pair<double, ibex::IntervalVector>> Vstate;

    class CtcVnode : public Ctc {
    public:

        CtcVnode();

        void Contract(vnodelp::AD *ad,double t,double tend,int n, Tube &x,double initial_time=-1);

        void Contract(vnodelp::AD *ad,double t,double tend,int n, TubeVector &x, double initial_time=-1);

    };
}


#endif
