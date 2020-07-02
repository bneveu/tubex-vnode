#include <utility>

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

    void CtcVnodelp::OdeContractor(vnodelp::AD *ad, double t, double tend,int n, Tube &x, Vstate state){
        Interval domain=Interval(t,tend);
        TubeVector y(domain);
        y[0]=x;
        OdeContractor(ad,t,tend,n,y, std::move(state));
        x=y[0];
    }//ctcvnodelp

    void CtcVnodelp::OdeContractor(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, Vstate state){
        //  assert(!state.empty());
        TubeVector oldx(x);

        for(int i=0; i<state.size(); i++){
            state[i].second&=x(state[i].first);
        }//for
        // cout <<"["<<t<<","<<tend<<"] " << "Forward step " <<endl;
        //forward step
        TubeVector fwd_x(x.domain(),n);
        fwd_x.set(state[0].second, state[0].first);
        int fwd(0);
        for( fwd; fwd<state.size()-1; fwd++){

            if(state[fwd].second.is_unbounded())
                continue;

            //  cout<<"fwd loop : "<<endl;
            //guaranteed integration
            //integration parameters
            //time
            interval ti=interval(state[fwd].first);
            interval tnext=interval(state[fwd+1].first);
            //initial conditions
            vnodelp::iVector y(n);
            for(int i=0; i<n; i++){
                y[i]=interval(state[fwd].second[i].lb(),state[fwd].second[i].ub());
            }//for
            //solver
            vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
            vnode_solver->setOneStep(vnodelp::on);
            vnode_solver->setOrder(ord);
            //vnode_solver->setHmin(1e-4);
            //integration

            while(ti!=tnext) {
            cout << "tnext : "<< tnext << endl;
                vnode_solver->integrate(ti, y, tnext);

                if (!vnode_solver->successful()){
                    cout << "VNODE-LP could not reach t = " << tnext << endl;//need to figure out what exception we can add to this case
                    x.set_empty();
                    x=oldx;

                    return;
                }
                for (int ith_sol=0; ith_sol<n; ith_sol++)
                    cout <<tnext <<"solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                //global enclosure and time step
                iVector Y = vnode_solver->getAprioriEncl();
                interval Tj = vnode_solver->getT();

                //forward tube
                IntervalVector door(n);
                IntervalVector boxe(n);
                for(int k=0; k<n; k++) {
                    door[k]=Interval(inf(y[k]), sup(y[k]));
                    boxe[k]=Interval(inf(Y[k]), sup(Y[k]));

                }
                door.intersects(x(sup(ti)));
                if(sup(Tj)<sup(tnext))
                    fwd_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                else
                    fwd_x.set(boxe,Interval(inf(Tj),sup(tnext)));
                fwd_x.set(door,sup(ti));

                //update of y
                for(int k=0; k<n; k++) {
                    y[k] = interval(door[k].lb(), door[k].ub());
                }
            }//while

            //update of the state vector
            for (int i=0; i<n; i++){
//            cout<<"solution at time t = "<< ti <<" is : " << y[i] <<endl;
//            cout << "from obs vector we have : "<<state[fwd+1].first<< " -> "<< state[fwd+1].second<<endl;
//            cout << "here state fwd+1"<< state[fwd+1].second[i] << endl;
//            cout << "here interval y"<< Interval(inf(y[i]),sup(y[i])) << endl;
                state[fwd+1].second[i]=state[fwd+1].second[i]&Interval(inf(y[i]),sup(y[i]));
//            cout << "here "<< state[fwd+1].second[i] << endl;
                //assert(!(state[fwd+1].second[i].is_empty()));
                if(state[fwd+1].second[i].is_empty()) {
                    x.set_empty();
                    return;
                }
            }//for
//            cout << "from intersection obs vector becomes : "<<state[fwd+1].first<< " -> "<< state[fwd+1].second<<endl;

            delete vnode_solver;
        }//for
        //
        if(state[fwd].first<tend){
            //    cout << "check tend : "<<endl;
            //guaranteed integration
            //integration parameters
            interval ti=interval(state[fwd].first);
            interval tnext=interval(tend);
            //initial conditions
            vnodelp::iVector y(n);
            for(int i=0; i<n; i++){
                y[i]=interval(state[fwd].second[i].lb(),state[fwd].second[i].ub());
            }
            //solver
            vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
            vnode_solver->setOneStep(vnodelp::on);
            vnode_solver->setOrder(ord);
            //integration
            while(ti!=tnext) {

                vnode_solver->integrate(ti, y, tnext);

                if (!vnode_solver->successful()){
                    cout << "VNODE-LP could not reach t = " << tnext << endl;//need to figure out what exception we can add to this case
                    x.set_empty();
                    x=oldx;
                    return;
                }
                for (int ith_sol=0; ith_sol<n; ith_sol++)
                    cout << "solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                //forward tube from last information to tend
                iVector Y = vnode_solver->getAprioriEncl();
                interval Tj = vnode_solver->getT();
                IntervalVector door(n);
                IntervalVector boxe(n);
                for(int k=0; k<n; k++) {
                    door[k]=Interval(inf(y[k]), sup(y[k]));
                    boxe[k]=Interval(inf(Y[k]), sup(Y[k]));

                }
                door&=x(sup(ti));
                if(sup(Tj)<sup(tnext))
                    fwd_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                else
                    fwd_x.set(boxe,Interval(inf(Tj),sup(tnext)));
                fwd_x.set(door,sup(ti));

                //update of y
                for(int k=0; k<n; k++)
                    y[k]=interval(door[k].lb(),door[k].ub());

            }//while
            for (int i=0; i<n; i++)
//                cout<<"solution at time t = "<< ti <<" is : " << y[i] <<endl;
                delete vnode_solver;
        }//if

        // cout << "Backward loop : "<<endl;
        //bwd
        TubeVector bwd_x(x.domain(),n);
        bwd_x.set(state[state.size()-1].second, state[state.size()-1].first);
        int bwd(state.size()-1);
        for( bwd; bwd>0; bwd--){

            if(state[bwd].second.is_unbounded())
                continue;

            // cout << "bwd loop : " << endl;

            //guaranteed integration
            //integration parameters
            //time
            interval ti=interval(state[bwd].first);
            interval tprevious=interval(state[bwd-1].first);
            //initial conditions
            vnodelp::iVector y(n);
            for(int i=0; i<n; i++){
                y[i]=interval(state[bwd].second[i].lb(),state[bwd].second[i].ub());
            }//for
            //solver
            vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
            vnode_solver->setOneStep(vnodelp::on);
            vnode_solver->setOrder(ord);
            //integration
            while(ti!=tprevious) {

                vnode_solver->integrate(ti, y, tprevious);

                if (!vnode_solver->successful()){
                    cout << "VNODE-LP could not reach t = " << ti << endl;//need to figure out what exception we can add to this case
                    x.set_empty();
                    x=oldx;
                    return;
                }
                for (int ith_sol=0; ith_sol<n; ith_sol++)
                    cout << "solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                //global enclosure and time step
                iVector Y = vnode_solver->getAprioriEncl();
                interval Tj = vnode_solver->getT();

                //bwd tube
                IntervalVector door(n);
                IntervalVector boxe(n);
                for(int k=0; k<n; k++) {
                    door[k]=Interval(inf(y[k]), sup(y[k]));
                    boxe[k]=Interval(inf(Y[k]), sup(Y[k]));
//                    cout << "taille boxe "<<boxe.diam()<<endl;


                }//for
                door&=x(inf(ti));
                if(inf(Tj)>inf(tprevious))
                    fwd_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                else
                    fwd_x.set(boxe,Interval(inf(tprevious),sup(Tj)));
                fwd_x.set(door,inf(ti));

                //update of y
                for(int k=0; k<n; k++)
                    y[k]=interval(door[k].lb(),door[k].ub());

            }//while

            //update of state vector
            for (int i=0; i<n; i++){
//                cout<<"solution at time t = "<< ti <<" is : " << y[i] <<endl;
//                cout << "from obs vector we have : "<<state[bwd-1].first<< " -> "<< state[bwd-1].second<<endl;
                state[bwd-1].second[i]=state[bwd-1].second[i]&Interval(inf(y[i]),sup(y[i]));
                //assert(!(state[bwd-1].second[i].is_empty()));
                if((state[bwd-1].second[i].is_empty())) {
                    x.set_empty();
                    return;
                }
            }//for
//            cout << "from intersection obs vector becomes : "<<state[bwd-1].first<< " -> "<< state[bwd-1].second<<endl;

            delete vnode_solver;
        }//for
        //
        if(state[bwd].first>t){
            cout<< "state bwd first "<<state[bwd].first << " ? t " << t<< endl;
            //   cout << "check t0 : "<<endl;
            //guaranteed integration
            //integration parameters
            interval ti=interval(state[bwd].first);
            interval tprevious=interval(t);
            //initial condition
            vnodelp::iVector y(n);
            for(int i=0; i<n; i++){
                y[i]=interval(state[bwd].second[i].lb(),state[bwd].second[i].ub());
            }//for
            //solver
            vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
            vnode_solver->setOneStep(vnodelp::on);
            vnode_solver->setOrder(ord);
            //integration
            while(ti!=tprevious) {

                vnode_solver->integrate(ti, y, tprevious);

                if (!vnode_solver->successful()){
                    cout << "VNODE-LP could not reach t = " << ti << endl;//need to figure out what exception we can add to this case
                    x.set_empty();
                    x=oldx;
                    return;
                }
                for (int ith_sol=0; ith_sol<n; ith_sol++)
                    cout << "solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                //global enclosure and time step
                iVector Y = vnode_solver->getAprioriEncl();
                interval Tj = vnode_solver->getT();

                //bwd tube from last information to t0
                IntervalVector door(n);
                IntervalVector boxe(n);
                for(int k=0; k<n; k++) {
                    door[k]=Interval(inf(y[k]), sup(y[k]));
                    boxe[k]=Interval(inf(Y[k]), sup(Y[k]));
//                    cout << "taille boxe "<<boxe.diam()<<endl;

                }//for
                door&=x(inf(ti));
                if(inf(Tj)>inf(tprevious))
                    bwd_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                else
                    bwd_x.set(boxe,Interval(inf(tprevious),sup(Tj)));
                bwd_x.set(door,inf(ti));

                //update of y
                for(int k=0; k<n; k++)
                    y[k]=interval(door[k].lb(),door[k].ub());
            }//while
            for (int i=0; i<n; i++)
//                cout<<"solution at time t = "<< ti <<" is : " << y[i] <<endl;

                delete vnode_solver;
        }//if
        //fusion of initial tube x, and transition tubes fwd_x and bwd_x
        x=x&fwd_x;
        x=x&bwd_x;
    }//odecontract

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
    }

    void CtcVnodelp::Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x){

        cout << "Vnode" << endl;

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
    for (int i=0; i<n; i++) {
        for (int j = 0; j < time_gate.size(); j++) {
            // cout<<endl;
            cout << "time: " << time_gate[j] << endl;//
            cout << " -> state : "<<si[j][i]<<endl;
        }
    }
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
    }//contract
}//namespace
