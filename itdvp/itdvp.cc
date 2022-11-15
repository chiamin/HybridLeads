#include <typeinfo>
#include <iomanip>
#include "itensor/all.h"
#include "IUtility.h"
#include "uGauge.h"
#include "GenMPO.h"
#include "FixedPointTensor.h"
#include "Solver.h"
#include "iTDVP.h"
using namespace itensor;
using namespace std;



int main(int argc, char* argv[])
{
    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto t           = input.getReal("t");
    auto mu          = input.getReal("mu");
    auto V           = input.getReal("V");

    auto time_steps  = input.getInt("time_steps");
    auto dt_str      = input.getString("dt");
    auto dt = ((dt_str == "inf" || dt_str == "Inf" || dt_str == "INF") ? INFINITY : stod(dt_str));
    auto D           = input.getInt("max_dim");
    auto ErrGoal     = input.getReal("ErrGoal");
    auto MaxIter     = input.getInt("MaxIter");
    auto ErrGoalInit = input.getReal("ErrGoalInit",1e-12);
    auto MaxIterInit = input.getInt("MaxIterInit",100);
    auto SeedInit    = input.getInt("SeedInit",0);
    auto write       = input.getYesNo("write");
    auto out_dir     = input.getString("out_dir",".");

    //------------------------------------------
    cout << setprecision(18) << endl;
    // Make MPO
    auto sites = Fermion (3, {"ConserveQNs",false});
    // ************************************************************************
    // Change the Hamiltonian to whatever you want
    auto H = single_empurity_mpo (sites, 2, t, t, t, mu, mu, mu, V, V, V);
    // ************************************************************************
    auto [W, is, iwl, iwr] = get_W (H);
    // W: MPO tensor
    // is: physical index
    // iwl: left index
    // iwr: right index

    // Initialize MPS
    auto A = ITensor(); // ill-defined tensor
    auto [AL, AR, AC, C, La0, Ra0] = itdvp_initial (W, is, iwl, iwr, A, D, ErrGoalInit, MaxIterInit, SeedInit);
    // If A is ill-defined, A will be random generated in itdvp_initial
    // This is just for tensors to be used in iTDVP

    // iTDVP
    Args args = {"ErrGoal=",1e-4,"MaxIter",MaxIter};
    ITensor LW, RW;
    Real en, err;
    for(int i = 1; i <= time_steps; i++)
    {
        cout << "time step " << i << endl;
        // Run iTDVP
        // If dt is real,       do imaginary time evolution
        //          imaginary,     real
        tie (en, err, LW, RW) = itdvp (W, AL, AR, AC, C, La0, Ra0, dt, args);
        cout << "energy, error = " << en << " " << err << endl;

        // Decrease the ErrGoal dynamically
        if (args.getReal("ErrGoal") > ErrGoal)
            args.add("ErrGoal=",err*0.1);
    }
    //------------------------------------------

    if (write)
    {
        writeToFile (out_dir+"/AL.itensor",AL);
        writeToFile (out_dir+"/AR.itensor",AR);
        writeToFile (out_dir+"/AC.itensor",AC);
        writeToFile (out_dir+"/C.itensor",C);
        writeToFile (out_dir+"/LW.itensor",LW);
        writeToFile (out_dir+"/RW.itensor",RW);
        writeToFile (out_dir+"/W.itensor",W);
        writeToFile (out_dir+"/La.itensor",La0);
        writeToFile (out_dir+"/Ra.itensor",Ra0);
        global::IS.write (out_dir+"/global.inds");
    }
   return 0;
}
