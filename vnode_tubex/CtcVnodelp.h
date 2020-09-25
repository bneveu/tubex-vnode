#ifndef __CTCVNODELP_H__
#define __CTCVNODELP_H__

#include "ibex.h"
#include "tubex.h"
#include "tubex_DynCtc.h"
#include "tubex_Slice.h"
#include "vnode.h"
#include <vector>

using namespace ibex;

namespace tubex {
    typedef std::vector<std::pair<double, ibex::IntervalVector>> Vstate;

    class CtcVnodelp : public DynCtc {
    public:

        CtcVnodelp();

        void set_vnode_order(unsigned int vorder);

        void set_vnode_hmin(double vhmin);

        void set_vnode_tol(double vatol, double vrtol);

        void set_ignoreslicing(bool ignorslicing);


        void Contract(vnodelp::AD *ad,double t,double tend,int n, Tube &x, double t0, bool incremental);

        void Contract(vnodelp::AD *ad,double t,double tend,int n, TubeVector &x,double t0, bool incremental);
    private:

        void fill_state_vector(Tube &x,vector<double> &time_gate ,vector<ibex::IntervalVector>& si, const int& i);
        void starting_condition(int n, TubeVector &x,double t, double tend, double t0,vector <double> &time_gate,vector<ibex::IntervalVector> &si,  Vstate &gates_vector, int &starting_index,bool incremental);
        bool  vnode_integration(vnodelp::AD *ad,const int n, TubeVector &x,TubeVector &transition_x, const double &t, const double &tend,const  int &starter_index, const Vstate &gates_vector, const bool &direction);

    double Vord;
    double Vhmin;
    bool ishminset;
    double Vatol;
    double Vrtol;
    bool IgnoreSlicing;


    };
}


#endif
