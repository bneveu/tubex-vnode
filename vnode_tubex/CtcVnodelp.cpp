#include "CtcVnodelp.h"

namespace codac {

    CtcVnodelp::CtcVnodelp() : DynCtc() {
        //default Vnode order is 20, we set it here to 11
        Vord = 11;
	Vatol = 1e-12;
        Vrtol = 1e-12;
        IgnoreSlicing = false;
	Vhmin = 1e-3;
	ishminset = true;
    }//ctcvnode

    void CtcVnodelp::set_vnode_order(unsigned int vorder) {
        Vord = vorder;
    }

    void CtcVnodelp::set_vnode_hmin(double vhmin) {
        Vhmin = vhmin;
        ishminset = true;
    }

    void CtcVnodelp::set_vnode_tol(double vatol, double vrtol) {
        Vatol = vatol;
        Vrtol = vrtol;
    }

    void CtcVnodelp::set_ignoreslicing(bool ignoreslicing) {
        IgnoreSlicing = ignoreslicing;
    }
     void CtcVnodelp::disable_hmin() {
        ishminset=false;
    }


    void CtcVnodelp::starting_condition(int n, TubeVector &x, double t, double tend, double t0,vector<double> &time_gate,vector<ibex::IntervalVector> &si,Vstate &gates_vector, int &starter_index, bool incremental){

        IntervalVector empty(n);


        bool starter_index_exists = false;

        ////////////////////// incrémental = 1 //////////////////////////////
        if (incremental) {
            time_gate.emplace_back(std::nan("1"));
            si.emplace_back(empty);
            if (IgnoreSlicing) {
                starter_index=0;
                for (int i = 0; i < n; i++) {
                    Slice *x_slice = x[i].slice(t0);
                    if (x_slice->tdomain().lb() == t0) {
                        Interval ingate = x_slice->input_gate();
                        si[0][i] = ingate;
                        time_gate[0] = t0;

                    } else if (x_slice->tdomain().ub() == t0) {
                        Interval outgate=x_slice->output_gate();
                        si[0][i] = outgate;
                        time_gate[0] = t0;
                    } else {
                        cout << "Can not find a gate for time bisection" << endl;
                        incremental=false;
                        break;
                    }

                }
                starter_index_exists = true;
                gates_vector.emplace_back(time_gate[0], si[0]);
            }
            else if (!IgnoreSlicing) {
                for (int i = 0; i < x[0].nb_slices() + 1; i++) {
                    time_gate.emplace_back(std::nan("1"));
                    si.emplace_back(empty);
                }
                for (int i = 0; i < n; i++){
                    fill_state_vector(x[i], time_gate, si, i);
                }
                for (int i = 0; i < time_gate.size(); i++) {
                    gates_vector.emplace_back(time_gate[i], si[i]);
                    if (gates_vector[i].first == t0) {
                        starter_index = i;
                        starter_index_exists = true;
                    }
                }
                if(!starter_index_exists)
                    incremental = false;
            }
        }
            ////////////////////// incrémental = 0 //////////////////////////////
        else if(!incremental) {
            if (IgnoreSlicing) {

                double tstart=t0;
                IntervalVector istart(n);

                tstart= x[0].first_slice()->tdomain().lb();
                for (int j=0; j<n; j++){
                    istart[j]=x[j].first_slice()->input_gate();
                }
                double air=0; double minair=0;
                for (int j=0; j<n; j++){
                    air+=x[j].first_slice()->input_gate().diam();
                }
                minair=air;
                Slice* s[n];
                for (int j=0; j<n; j++)
                    s[j]=x[j].first_slice();
                for(const Slice *si = s[0] ; si != NULL ; si = si->next_slice()){
                    air=0;
                    for (int j=0; j<n; j++){
                        air+=s[j]->output_gate().diam();
                    }
                    if (air < minair){
                        minair=air; tstart=si->tdomain().ub();
                        for (int j=0; j<n; j++) {
                            istart[j]=s[j]->output_gate();
                        }
                    }
                    for (int j=0; j<n; j++)
                        s[j]=s[j]->next_slice();
                }
                gates_vector.emplace_back(tstart, istart);
            }
            else if(!IgnoreSlicing) {
                for (int i = 0; i < x[0].nb_slices() + 1; i++) {
                    time_gate.emplace_back(std::nan("1"));
                    si.emplace_back(empty);
                }
                for (int i = 0; i < n; i++) {
                    fill_state_vector(x[i], time_gate, si, i);
                }
                for (int i = 0; i < time_gate.size(); i++) {
                    gates_vector.emplace_back(time_gate[i], si[i]);
                    if (i == 0)
                        starter_index = i;
                    else {
                        double air(0);
                        double min_air(0);
                        for (int j = 0; j < n; j++) {
                            air += gates_vector[i].second[j].diam();
                            min_air += gates_vector[starter_index].second[j].diam();
                        }

                        if (air < min_air) {
                            starter_index = i;
                        }
                    }
                }
            }
        }

    }

    void CtcVnodelp::fill_state_vector(Tube &x, vector<double> &time_gate, vector<ibex::IntervalVector> &si, const int &i) {
        codac::Slice *x_slice = x.first_slice();
        int k = 0;
        while (x_slice != NULL) {
            Interval outgate = x_slice->output_gate();
            Interval dom = x_slice->tdomain();
            if (k == 0) {
                Interval ingate = x_slice->input_gate();
                si[k][i] = ingate;

                time_gate[k] = dom.lb();
                k++;
            }
            si[k][i] = outgate;
            time_gate[k] = dom.ub();
            k++;
            x_slice = x_slice->next_slice();
        }
    }

    bool CtcVnodelp::vnode_integration(vnodelp::AD *ad, const int n, TubeVector &x, TubeVector &transition_x,
                                       const double &t, const double &tend, const int &starter_index,
                                       const Vstate &gates_vector, const bool &direction) {
        //ODE Vnodelp
        vnodelp::VNODE *vnode_solver = new vnodelp::VNODE(ad);
        vnode_solver->setOneStep(vnodelp::on);
        vnode_solver->setOrder(Vord);
        if (ishminset)
            vnode_solver->setHmin(Vhmin);


        //main integration loop
        interval final_time;
        int k = 1;
        iVector y(n);
        interval ti;
        interval tnext;
        if (direction == true) {
            final_time = interval(tend);
            //Time and initial condition
            ti = interval(gates_vector[starter_index].first);
            if (IgnoreSlicing) {
                tnext = interval(tend);
            } else {
                tnext = interval(gates_vector[starter_index + k].first);
            }
            //initial condition
            for (int i = 0; i < n; i++) {
                y[i] = interval(gates_vector[starter_index].second[i].lb(),
                                gates_vector[starter_index].second[i].ub());
            }
        } else if (direction == false) {
            final_time = interval(t);
            //Time and initial condition
            ti = interval(gates_vector[starter_index].first);
            if (IgnoreSlicing) {
                tnext = interval(t);
            } else {
                tnext = interval(gates_vector[starter_index - k].first);
            }
            //initial condition

            for (int i = 0; i < n; i++) {
                y[i] = interval(gates_vector[starter_index].second[i].lb(),
                                gates_vector[starter_index].second[i].ub());
            }
        }
        while (ti != final_time) {
            //"slice" integration loop
            while (ti != tnext) {

                vnode_solver->integrate(ti, y, tnext);
	       
                if (!vnode_solver->successful()) {
		  //                    cout << "VNODE-LP could not reach t = " << tnext << endl;
                    transition_x.set_empty();
                    delete vnode_solver;
                    return false;
                }
//                for (int ith_sol = 0; ith_sol < n; ith_sol++)
//                    cout << "solution enclosure at t = " << ti << " [" << inf(y[ith_sol]) << ", " << sup(y[ith_sol])
//                         << "]" << endl;

                //global enclosure and time step
                iVector Y = vnode_solver->getAprioriEncl();
                interval Tj = vnode_solver->getT();

                //putting results in fwd_x/bwd_x

                //parameters to save data
                IntervalVector door(n);
                IntervalVector boxe(n);

                //results conversion vnode interval -> ibex interval
                for (int i = 0; i < n; i++) {
                    door[i] = Interval(inf(y[i]), sup(y[i]));
                    boxe[i] = Interval(inf(Y[i]), sup(Y[i]));
                }

                if (direction == true) {
                    if (!door.intersects(x(sup(ti)))) {
                        x.set_empty();
                        delete vnode_solver;
                        return false;
                    }

                    //saving data in tubes
                    if (sup(Tj) < sup(tnext))
                        transition_x.set(boxe, Interval(inf(Tj), sup(Tj)));
                    else
                        transition_x.set(boxe, Interval(inf(Tj), sup(tnext)));

                    transition_x.set(door, sup(ti));
                } else if (direction == false) {
                    if (!door.intersects(x(inf(ti)))) {
                        x.set_empty();
                        delete vnode_solver;
                        return false;
                    }

                    //saving data in tubes
                    if (inf(Tj) > inf(tnext))
                        transition_x.set(boxe, Interval(inf(Tj), sup(Tj)));
                    else
                        transition_x.set(boxe, Interval(inf(tnext), sup(Tj)));

                    transition_x.set(door, inf(ti));
                }
            }
            k++;
            if (direction == true) {
                if (k < gates_vector.size()) {
                    tnext = gates_vector[starter_index + k].first;

                }
            } else if (direction == false) {
                if (starter_index - k >= 0) {
                    tnext = gates_vector[starter_index - k].first;

                }
            }
        }

        delete vnode_solver;
        return true;
    }

    void CtcVnodelp::Contract(vnodelp::AD *ad, double t, double tend, int n, Tube &x, double t0, bool incremental) {

        ibex::Interval domain = ibex::Interval(t, tend);
        TubeVector y(domain, n);
        y[0] = x;
        Contract(ad, t, tend, n, y, t0, incremental);
        x = y[0];
    }//

    void    CtcVnodelp::Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, double t0, bool incremental) {

        Interval domain(t, tend);
	//        cout << "vnodelp..." << domain << endl;

        //initialisation of vstate (get time discretization of x and gates)
        Vstate gates_vector;
        vector<ibex::IntervalVector> si;
        vector<double> time_gate;
        int starter_index(0);

        starting_condition(n,x,t,tend,t0,time_gate,si,gates_vector,starter_index,incremental);
        //////////////////////////////////////////////
	//	cout << " gate " << gates_vector[starter_index].second << " t " << gates_vector[starter_index].first << endl;
        TubeVector fwd_x(domain,n);
        fwd_x.set(gates_vector[starter_index].second, gates_vector[starter_index].first);
        TubeVector bwd_x(domain,n);
        bwd_x.set(gates_vector[starter_index].second, gates_vector[starter_index].first);
        bool successfull_integration=true;

        ////////////////////////////////fwd integration///////////////////////////////////////

        //fwd tubevector
        if (gates_vector[starter_index].first<tend) {
            // integration parameters
            bool FWD=true;
            bool direction=FWD;

            successfull_integration=vnode_integration(ad,n,x, fwd_x,t,tend,starter_index,gates_vector,direction);

            if(!successfull_integration)
                return;
        }
        ////////////////////////////////bwd integration///////////////////////////////////////
        //bwd tubevector
        if (gates_vector[starter_index].first>t) {
            bool BWD=false;
            bool direction=BWD;

            successfull_integration=vnode_integration(ad,n,x, bwd_x,t,tend,starter_index,gates_vector,direction);

            if(!successfull_integration)
                return;
        }

        if(m_preserve_slicing)
            x&=fwd_x&bwd_x;
        else
            x = (fwd_x&bwd_x) & x;

    }

  void  CtcVnodelp::contract(std::vector<Domain*>& v_domains) {;}
}
