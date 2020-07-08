#include "CtcVnodelp.h"
//#include <vector>
namespace tubex {
    CtcVnodelp::CtcVnodelp() : DynCtc() {
    }//ctcvnode

	    void CtcVnodelp::fill_state_vector( Tube &x,vector<double> &time_gate ,vector<ibex::IntervalVector>& si, const int& i){
		tubex::Slice *x_slice=x.first_slice();
        int k=0;
        while (x_slice!=NULL){
            Interval outgate=x_slice->output_gate();
            Interval dom=x_slice->tdomain();
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

    void CtcVnodelp::Contract(vnodelp::AD *ad, double t, double tend, int n, Tube &x, double t0,bool incremental) {

        ibex::Interval domain = ibex::Interval(t, tend);
        TubeVector y(domain,n);
        y[0] = x;
        Contract( ad, t, tend, n, y,t0, incremental);
        x=y[0];
    }//

    void CtcVnodelp::Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, double t0,bool incremental) {
        //check if initial time is given

//        assert(!x.is_empty());
//	cout << "vnodelp..."<<endl;
        Interval domain(t, tend);
        TubeVector tempo_x(domain,n);


//        TubeVector fwd_x(domain, n);
//        TubeVector bwd_x(domain, n);


        //initialisation of vstate (get time discretization of x and gates)
        Vstate gates_vector;
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

        int starter_index(0);
        bool starter_index_exists=false;
        //bisection time starter index for integration
        if(incremental==true){
            for (int i = 0; i < time_gate.size(); i++) {
                gates_vector.emplace_back(time_gate[i], si[i]);
                if(gates_vector[i].first==t0) {
                    starter_index = i;
                    starter_index_exists = true;
                   // cout << "incrémental"<<endl;
                }

            }

        }
        if(starter_index_exists==false){
            incremental=false;
        }
        //no bisection time -> smallest gate
        if(incremental==false){
          //  cout << "non incremental"<< endl;

            for (int i = 0; i < time_gate.size(); i++) {
                gates_vector.emplace_back(time_gate[i], si[i]);
//                cout << "time : " << gates_vector[i].first << " - > [ " << gates_vector[i].second << " ] --> vol : " << gates_vector[i].second.volume()
                  //   << endl;
                if (i == 0)
                    starter_index = i;
                else {
			double air(0);
			double min_air(0);
			for (int j=0; j<n; j++){
			air+=gates_vector[i].second[j].diam();
			min_air+=gates_vector[starter_index].second[j].diam();
			}	
                   // if (gates_vector[i].second.volume() < gates_vector[min_starter_index].second.volume())
                     if(air<min_air) {
                         starter_index = i;
                     }
                }
            }

//            cout << "min vol index : " << min_starter_index << "- > " << gates_vector[min_starter_index].first<<" : " <<gates_vector[min_starter_index].second << endl;

        }

        tempo_x.set(gates_vector[starter_index].second, gates_vector[starter_index].first);





        ////////////////////////////////fwd integration///////////////////////////////////////

        //fwd tubevector
	//    cout << "# "<< min_starter_index << " ## "<<gates_vector.size()<< endl;
        if (starter_index<gates_vector.size()-1) {
            //fwd_x.set(gates_vector[min_starter_index].second, gates_vector[min_starter_index].first);
            // integration parameters

            //time
            interval ti = interval(gates_vector[starter_index].first);
            int k=1;
            interval tnext = interval(gates_vector[starter_index + k].first);

            //initial condition
            iVector y(n);
            for (int i = 0; i < n; i++) {
                y[i] = interval(gates_vector[starter_index].second[i].lb(),
                                gates_vector[starter_index].second[i].ub());
            }

            //integration object
            vnodelp::VNODE *vnode_fwd_solver = new vnodelp::VNODE(ad);
            vnode_fwd_solver->setOneStep(vnodelp::on);
            vnode_fwd_solver->setOrder(ord);
	    vnode_fwd_solver->setHmin(0.0005);
            //main integration loop
            while(ti!=interval(tend)){
                //"slice" integration loop
                while(ti!=tnext){
                    vnode_fwd_solver->integrate(ti,y,tnext);
                    if (!vnode_fwd_solver->successful()) {
                        cout << "VNODE-LP could not reach t = " << tnext
                             << endl;//need to figure out what exception we can add to this case
                        //tempo_x.set_empty();
                        return ;
                    }
                    //printVector(y);
		    //                for (int ith_sol=0; ith_sol<n; ith_sol++)
		  //                    cout <<"solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                    //global enclosure and time step
                    iVector Y = vnode_fwd_solver->getAprioriEncl();
                    interval Tj = vnode_fwd_solver->getT();

                    //putting results in fwd_x

                    //parameters to save data
                    IntervalVector door(n);
                    IntervalVector boxe(n);
                    for(int k=0; k<n; k++) {
                        door[k]=Interval(inf(y[k]), sup(y[k]));
                        boxe[k]=Interval(inf(Y[k]), sup(Y[k]));

                    }
            //        cout <<sup(ti)<< "door : "<<door << " -> "<< door.intersects(x(sup(ti)))<<endl;
//                    door.intersects(x(sup(ti)));
//                    cout <<sup(ti)<< "door : "<<door <<endl;

                    //test that the computed gate is still in the tube
                  //  for(int e=0; e<n; e++) {
                        if (!door.intersects(x(sup(ti)))) {
                            x.set_empty();
                          //  cout << "EMPTY"<<endl;
                            return;
                        }
                //    }
                    //
                    if(sup(Tj)<sup(tnext))
                        //fwd_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                        tempo_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                    else
                        //fwd_x.set(boxe,Interval(inf(Tj),sup(tnext)));
                        tempo_x.set(boxe,Interval(inf(Tj),sup(tnext)));
                    //fwd_x.set(door,sup(ti));
                    tempo_x.set(door,sup(ti));

                }



                k++;
//                cout << "k "<< k << endl;
                if(k<gates_vector.size())
                tnext=gates_vector[starter_index + k].first;

            }
        }
        ////////////////////////////////bwd integration///////////////////////////////////////
        //bwd tubevector
       // cout << "# "<< min_starter_index << " ## "<<gates_vector.size()<< endl;
        if (starter_index>0) {
           // cout << "bwd"<<endl;
           // bwd_x.set(gates_vector[min_starter_index].second, gates_vector[min_starter_index].first);
            // integration parameters

            //time
            interval ti = interval(gates_vector[starter_index].first);
            int k=-1;
            interval tprev = interval(gates_vector[starter_index + k].first);

            //initial condition
            iVector y(n);
            for (int i = 0; i < n; i++) {
                y[i] = interval(gates_vector[starter_index].second[i].lb(),
                                gates_vector[starter_index].second[i].ub());
            }

            //integration object
            vnodelp::VNODE *vnode_bwd_solver = new vnodelp::VNODE(ad);
            vnode_bwd_solver->setOneStep(vnodelp::on);
            vnode_bwd_solver->setOrder(ord);
		vnode_bwd_solver->setHmin(0.0005);
            //main integration loop
            while(ti!=interval(t)){
                //"slice" integration loop
                while(ti!=tprev){
                    vnode_bwd_solver->integrate(ti,y,tprev);
                    if (!vnode_bwd_solver->successful()) {
                        cout << "VNODE-LP could not reach t = " << tprev
                             << endl;//need to figure out what exception we can add to this case
                        tempo_x.set_empty();
                        return ;
                    }
  //                  for (int ith_sol=0; ith_sol<n; ith_sol++)
    //                    cout <<"solution enclosure at t = " << inf(ti) << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;

                    //global enclosure and time step
                    iVector Y = vnode_bwd_solver->getAprioriEncl();
                    interval Tj = vnode_bwd_solver->getT();


                    //putting results in bwd_x

                    //parameters to save data
                    IntervalVector door(n);
                    IntervalVector boxe(n);
                    for(int k=0; k<n; k++) {
                        door[k]=Interval(inf(y[k]), sup(y[k]));
                        boxe[k]=Interval(inf(Y[k]), sup(Y[k]));

                    }

                        if (!door.intersects(x(inf(ti)))) {
                            x.set_empty();
                           // cout << "EMPTY"<<endl;

                            return;
                        }
                    if(inf(Tj)>inf(tprev))
                        tempo_x.set(boxe,Interval(inf(Tj),sup(Tj)));
                    else
                        tempo_x.set(boxe,Interval(inf(tprev),sup(Tj)));
                    tempo_x.set(door,inf(ti));

                }



                k--;
      //          cout << "k "<< k << endl;
                if(starter_index + k>=0)
                    tprev=gates_vector[starter_index + k].first;

            }

        }
  // cout << tempo_x << endl;
//cout << "x" << endl;
//        tubex::Slice *x_slice=x[0].first_slice();
//        while (x_slice!=NULL){
//            cout << *x_slice<<endl;
//
//            x_slice=x_slice->next_slice();
//        }
//        cout << "///////////////////////" << endl;
//        tubex::Slice *y_slice=x[1].first_slice();
//        while (y_slice!=NULL){
//            cout << *y_slice<<endl;
//
//            y_slice=y_slice->next_slice();
//        }
//        cout << "vnode x" << endl;
//
//        tubex::Slice *tx_slice=tempo_x[0].first_slice();
//        while (tx_slice!=NULL){
//            cout << *tx_slice<<endl;
//
//            tx_slice=tx_slice->next_slice();
//        }
//        cout << "///////////////////////" << endl;
//
//        tubex::Slice *ty_slice=tempo_x[1].first_slice();
//        while (ty_slice!=NULL){
//            cout << *ty_slice<<endl;
//
//            ty_slice=ty_slice->next_slice();
//        }
    x=tempo_x&x;
//cout << "############################# x&vnodex "<< endl;
//
//        tubex::Slice *xv_slice=x[0].first_slice();
//        while (xv_slice!=NULL){
//            cout << *xv_slice<<endl;
//
//            xv_slice=xv_slice->next_slice();
//        }
//        tubex::Slice *yv_slice=x[1].first_slice();
//        while (yv_slice!=NULL){
//            cout << *yv_slice<<endl;
//
//            yv_slice=yv_slice->next_slice();
//        }
//    cout << " t : "<< domain.lb()<< "-> " << x(domain.lb()) << " t : " <<domain.ub()<<"-> " << x(domain.ub()) << " tube : " << x <<endl;

    }
}

//x: 0.25,0.28
