#include "CtcVnode.h"
//#include <vector>
namespace tubex {
    CtcVnode::CtcVnode() : Ctc() {
    }//ctcvnode


    void CtcVnode::Contract(vnodelp::AD *ad, double t, double tend, int n, Tube &x, double initial_time) {

        ibex::Interval domain = ibex::Interval(t, tend);
        TubeVector y(domain);
        y[0] = x;
        Contract(ad, t, tend, n, y, initial_time);
        x=y[0];
    }//

    void CtcVnode::Contract(vnodelp::AD *ad, double t, double tend, int n, TubeVector &x, double initial_time) {
        //check if initial time is given
        if(initial_time==-1){
            initial_time=t;
            cout << "no initial time given starts at t =  "<< t << endl;
        }
        //check if x is empty at the initial time
        if (x(initial_time).is_empty()){
            cout << x(initial_time)<<endl;
            return;

        }
        //parameters for fwd transition tubevectors
        ibex::Interval domain = ibex::Interval(t, tend);
        TubeVector fwd_tempo_x(domain, n);

        //getting tubevector time discretizition

        int initial_time_index(-1); //parameter to check if the initial time given is in the tubevector
        int count_index=0; //parameter to save its position in a vactor

        Slice *s_x = x[0].first_slice();
        std::vector<double> time_dis;
        time_dis.emplace_back(t);

        if (initial_time == t)
            initial_time_index=count_index;

        while(s_x!=NULL){
            time_dis.emplace_back(s_x->tdomain().ub());
            count_index++;
            if(s_x->tdomain().ub()==initial_time)
                initial_time_index=count_index;;
            s_x = s_x->next_slice();
        }

        for(int i=0; i<time_dis.size(); i++){
            std::cout << time_dis[i] << " | ";
        }
        std::cout << endl <<"here : " <<initial_time_index << endl;
        if (initial_time_index==-1)
            return;
        //fwd resolution
        interval ti=interval(time_dis[initial_time_index]);
        interval tfinal=(time_dis[time_dis.size()-1]);
        int time_index=initial_time_index +1;
        iVector y(n);

        for (int i = 0; i<n; i++){
            y[i]=interval(x[i](inf(ti)).lb(), x[i](inf(ti)).ub());
        }
        //setting first gate for transition tubevector
        fwd_tempo_x.set(x(inf(ti)),inf(ti));

        vnodelp::VNODE *fwd_vnode_solver = new vnodelp::VNODE(ad);
        interval step =interval(time_dis[time_index]);


        while(ti!=tfinal){
            cout << "ti "<< ti << " step "<< step << endl;
            fwd_vnode_solver->setOneStep(vnodelp::on);
            fwd_vnode_solver->integrate(ti,y,step);

                for (int ith_sol=0; ith_sol<n; ith_sol++)
                    cout <<step <<"solution enclosure at t = " << ti << " ["<<inf(y[ith_sol]) << ", " <<sup(y[ith_sol])<<"]"<<endl;



            //global enclosure and time step
            iVector Y = fwd_vnode_solver->getAprioriEncl();
            interval Tj = fwd_vnode_solver->getT();

            //forward tube
            ibex::IntervalVector door(n);
            ibex::IntervalVector boxe(n);
            for(int k=0; k<n; k++) {
                door[k]=ibex::Interval(inf(y[k]), sup(y[k]));
                boxe[k]=ibex::Interval(inf(Y[k]), sup(Y[k]));

            }


            //next time step
            if(inf(ti)>=sup(step)) {
                time_index++;
                step = interval(time_dis[time_index]);
            }
        }


    }
}
