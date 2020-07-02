
    Interval domain(0.,10.);
    TubeVector x(domain, 1);
    IntervalVector v(1);
    double init = exp(-10);
    v[0]=Interval(init);

    x.set(v, 10.); // ini


    double eps=02;

    /* =========== SOLVER =========== */
    Vector epsilon(1, eps);


    tubex::Solver solver(epsilon);

    //      solver.set_refining_fxpt_ratio(0.9999);
//    solver.set_refining_fxpt_ratio(2.0);
    solver.set_propa_fxpt_ratio(0.999);
    solver.set_cid_fxpt_ratio(0.999);
    solver.set_cid_propa_fxpt_ratio(0.999);

    solver.set_cid_timept(2);
    solver.set_bisection_timept(2);
    solver.set_max_slices(20000);
    solver.set_refining_mode(0);
//    solver.set_trace(1);

    list<TubeVector> l_solutions = solver.solve(x, &contract);
    cout << l_solutions.front() << endl;
    cout << "nb sol " << l_solutions.size() << endl;
    return 0;#include <utility>

#include "tubex_CtcVnodelp.h"

using namespace std;
using namespace ibex;

namespace tubex {
    CtcVnodelp::CtcVnodelp() : Ctc() {

    }//ctcvnodelp

    void CtcVnodelp::set_state(Vstate &s){
    state = s;
    }//set

    void CtcVnodelp::show_state_info(){
        for(int i=0; i<state.size(); i++){
            cout <<"time state: " <<state[i].first;
            cout <<" - state: [ "<< state[i].second<<"]"<<endl;
        }//for
    }//show

    void CtcVnodelp::FwdBwdIntegration(vnodelp::AD *ad, double t, double tend, int n,const TubeVector &x,TubeVector &transition_x, Vstate &state,bool direction) {
        double final_time(tend);
        Vstate integration_states;
        if (direction==FWD) {
            //cout << "fwd" << endl;
                        for (int i = 0; i < state.size(); i++)
                integration_states.emplace_back(state[i]);
        }
        if (direction==BWD) {
           //cout << "bwd" << endl;
            final_time = t;
            for (int i = static_cast<int>(state.size() - 1); i >= 0; i--)
                integration_states.emplace_back(state[i]);
        }

        transition_x.set(integration_states[0].second, integration_states[0].first);
        int fc(0);
        //initial conditions
        vnodelp::iVector y(n);
        for (int i = 0; i < n; i++) {
            y[i] = interval(integration_states[fc].second[i].lb(), integration_states[fc].second[i].ub());
        }//for
        //solver
        vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
        vnode_solver->setOneStep(vnodelp::on);
        vnode_solver->setOrder(ord);

        while (integration_states[fc].first != final_time) {
            if (integration_states[fc].second.is_unbounded()){
                fc++;
                for (int i = 0; i < n; i++) {
                    y[i] = interval(integration_states[fc].second[i].lb(), integration_states[fc].second[i].ub());
		vnode_solver->setFirstEntry();
                }//for
                continue;
            }

            //time
            interval ti;
            interval tnext;
            ti = interval(integration_states[fc].first);
            if(fc!=integration_states.size()-1)
            tnext = interval(integration_states[fc + 1].first);
            else
            tnext = interval(final_time);

            if(ti==tnext)
                return;
           // cout << "time " << ti << " , " << tnext<<endl;

            //integration
            while (ti != tnext) {
                vnode_solver->integrate(ti, y, tnext);
                if (!vnode_solver->successful()) {
                    cout << "VNODE-LP could not reach t = " << tnext
                         << endl;//need to figure out what exception we can add to this case
                    transition_x.set_empty();
                    return ;
                }
//                for (int ith_sol=0; ith_sol<n; ith_sol++)
//                    cout <<tnext <<"solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                //global enclosure and time step
                iVector Y = vnode_solver->getAprioriEncl();
                interval Tj = vnode_solver->getT();

                if (direction==FWD) {
                    //forward tube
                    IntervalVector door(n);
                    IntervalVector boxe(n);
                    for (int k = 0; k < n; k++) {
                        door[k] = Interval(inf(y[k]), sup(y[k]));
                        boxe[k] = Interval(inf(Y[k]), sup(Y[k]));
                    }
                    door.intersects(x(sup(ti)));
                    if (sup(Tj) < sup(tnext))
                        transition_x.set(boxe, Interval(inf(Tj), sup(Tj)));
                    else
                        transition_x.set(boxe, Interval(inf(Tj), sup(tnext)));
                    transition_x.set(door, sup(ti));

                    //update of y
                    for (int k = 0; k < n; k++) {
                       //vérifier que suffisement bon?
			 y[k] = interval(door[k].lb(), door[k].ub());
                    }
                }
                if (direction==BWD) {
                    //bwd tube
                    IntervalVector door(n);
                    IntervalVector boxe(n);
                    for (int k = 0; k < n; k++) {
                        door[k] = Interval(inf(y[k]), sup(y[k]));
                        boxe[k] = Interval(inf(Y[k]), sup(Y[k]));
                    }//for
                    door &= x(inf(ti));
                    if (inf(Tj) > inf(tnext))
                        transition_x.set(boxe, Interval(inf(Tj), sup(Tj)));
                    else
                        transition_x.set(boxe, Interval(inf(tnext), sup(Tj)));
                    transition_x.set(door, inf(ti));

                    //update of y
                    for (int k = 0; k < n; k++)
                        y[k] = interval(door[k].lb(), door[k].ub());
                }
            }//while
            //update of the state vector
            for (int i = 0; i < n; i++) {
                integration_states[fc + 1].second[i] = integration_states[fc + 1].second[i] & Interval(inf(y[i]), sup(y[i]));

                if (integration_states[fc + 1].second[i].is_empty()) {
                    transition_x.set_empty();
                    return;
                }
            }
            fc++;
        }
        if(direction==FWD){
            for(int i=0; i<state.size(); i++){
                state[i]=integration_states[i];
            }
        }
        if(direction==BWD){
            for(int i=0; i<state.size(); i++){
                state[i]=integration_states[state.size()-1-i];
            }
        }

    }//integration garantie

    void CtcVnodelp::OdeContractor(vnodelp::AD *ad, double t, double tend,int n, Tube &x, Vstate state){
        Interval domain=Interval(t,tend);
        TubeVector y(domain);
        y[0]=x;
        OdeContractor(ad,t,tend,n,y, std::move(state));
        x=y[0];
    }//odecontractor


    void CtcVnodelp::OdeContractor(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, Vstate state){
        //update of the vector state
        for(int i=0; i<state.size(); i++){
            state[i].second&=x(state[i].first);
        }//for

        TubeVector fwd_x(x.domain(),n);

        FwdBwdIntegration(ad,t,tend,n,x,fwd_x,state,FWD);

        if(!fwd_x.is_empty()){
            x=x&fwd_x;

        }
        TubeVector bwd_x(x.domain(),n);
        FwdBwdIntegration(ad,t,tend,n,x,bwd_x,state,BWD);
        if(!bwd_x.is_empty())
        {
            x=x&bwd_x;

        }
    }//odecontract

//
//    void CtcVnodelp::GuaranteedIntegration(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, Vstate state) {
//        if(state.size()>1)
//            OdeContractor(ad,t,tend,n,x,state);
//        if(state[0].first>t)
//            if(state[0].first<tend)
//                OdeContractor(ad,t,tend,n,x,state);
//
//        if(state.size()==t) {
//            vnodelp::iVector y(n);
//            for (int i = 0; i < n; i++) {
//                y[i] = interval(state[0].second[i].lb(), state[0].second[i].ub());
//            }//for
//            //solver
//            vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
//            vnode_solver->setOneStep(vnodelp::on);
//            vnode_solver->setOrder(ord);
//            interval tinit(t);
//            interval tfinal(tend);
//            //integration
//            while (tinit != tfinal) {
//                vnode_solver->integrate(tinit, y, tfinal);
//                if (!vnode_solver->successful()) {
//                    cout << "VNODE-LP could not reach t = " << tend
//                         << endl;//need to figure out what exception we can add to this case
//                    x.set_empty();
//                    return ;
//                }
////                for (int ith_sol=0; ith_sol<n; ith_sol++)
////                    cout <<tnext <<"solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;
//
//                //global enclosure and time step
//                iVector Y = vnode_solver->getAprioriEncl();
//                interval Tj = vnode_solver->getT();
//
//                    //forward tube
//                    IntervalVector door(n);
//                    IntervalVector boxe(n);
//                    for (int k = 0; k < n; k++) {
//                        door[k] = Interval(inf(y[k]), sup(y[k]));
//                        boxe[k] = Interval(inf(Y[k]), sup(Y[k]));
//                    }
//                    door.intersects(x(sup(tinit)));
//                    if (sup(Tj) < sup(tend))
//                        x.set(boxe, Interval(inf(Tj), sup(Tj)));
//                    else
//                        x.set(boxe, Interval(inf(Tj), sup(tend)));
//                    x.set(door, sup(t));
//
//                    //update of y
//                    for (int k = 0; k < n; k++) {
//                        y[k] = interval(door[k].lb(), door[k].ub());
//                    }//for
//                }//while
//        }//if
//
//    }

    void CtcVnodelp::fill_state_vector( Tube &x,vector<double> &time_gate ,vector<ibex::IntervalVector>& si, const int& i){
        tubex::Slice *x_slice=x.first_slice();
        int k=0;
        while (x_slice!=NULL){
            Interval outgate=x_slice->output_gate();
            Interval dom=x_slice->domain();
            if(k==0) {
                Interval ingate = x_slice->input_gate();
                si[k][i] =ingate;
                //cout <<"[ "<<si[k][i].lb()<<" , "<< si[k][i].ub()<<" ]"<< endl;

                time_gate[k]=dom.lb();
                k++;
            }
            si[k][i] =outgate;
          //  cout <<"[ "<<si[k][i].lb()<<" , "<< si[k][i].ub()<<" ]"<< endl;
            time_gate[k]=dom.ub();
            k++;
            x_slice=x_slice->next_slice();
        }
    }//fill_state_vector

    void CtcVnodelp::Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x){

    Vstate s;
    vector<ibex::IntervalVector> si;
    vector <double>time_gate;
    IntervalVector empty(n);

    for(int i=0; i<=x[0].nb_slices();i++) {
        time_gate.emplace_back(std::nan("1"));
        si.emplace_back(empty);
    }

    for(int i=0; i<n; i++){
        fill_state_vector(x[i],time_gate,si,i);
    }
    //check
//    for (int i=0; i<n; i++) {
//        for (int j = 0; j < time_gate.size(); j++) {
//            // cout<<endl;
//            cout << "time: " << time_gate[j] << endl;//
//            cout << " -> state : "<<si[j][i]<<endl;
//        }
//    }
   // assert(time_gate.size()==si.size());

    for(int i=0; i<time_gate.size();i++) {
        s.emplace_back(time_gate[i],si[i]);
    }
    TubeVector y=x;
        if (x.is_empty())
            return;
    OdeContractor(ad,t,tend,n,y,s);
    if(x.is_empty())
        return;
    x=x&y;
    //cout <<" vol x = "<< x.volume() << endl;
    }


//contract
}//namespace
