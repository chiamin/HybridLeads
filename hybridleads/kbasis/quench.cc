#include <iomanip>

#include "Timer.h"
#include "itensor/all.h"
Timers timer;
#include "BdGBasis.h"
#include "ContainerUtility.h"
#include "Hamiltonian.h"
#include "IUtility.h"
#include "InitState.h"
#include "MPSUtility.h"
#include "MixedBasis.h"
#include "MyObserver.h"
#include "OneParticleBasis.h"
#include "ReadInput.h"
#include "ReadWriteFile.h"
#include "TDVPObserver.h"
#include "basisextension.h"
#include "tdvp.h"
using namespace itensor;
using namespace std;

// Define some parameters
struct Para {
  Real tcL = 0., tcR = 0.;

  void write(ostream& s) const {
    iut::write(s, tcL);
    iut::write(s, tcR);
  }

  void read(istream& s) {
    iut::read(s, tcL);
    iut::read(s, tcR);
  }
};

void writeAll(const string& filename, const MPS& psi, const MPO& H,
              const Para& para, const Args& args_basis, int step,
              const ToGlobDict& to_glob, const ToLocDict& to_loc) {
  ofstream ofs(filename);
  itensor::write(ofs, psi);
  itensor::write(ofs, H);
  itensor::write(ofs, args_basis);
  itensor::write(ofs, step);
  para.write(ofs);
  iut::write(ofs, to_glob);
  iut::write(ofs, to_loc);
}

void readAll(const string& filename, MPS& psi, MPO& H, Para& para,
             Args& args_basis, int& step, ToGlobDict& to_glob,
             ToLocDict& to_loc) {
  ifstream ifs = open_file(filename);
  itensor::read(ifs, psi);
  itensor::read(ifs, H);
  itensor::read(ifs, args_basis);
  itensor::read(ifs, step);
  para.read(ifs);
  iut::read(ifs, to_glob);
  iut::read(ifs, to_loc);
}

void print_orbs(const vector<SortInfo>& orbs) {
  cout << "Orbitals: name, ki, energy" << endl;
  for (int i = 1; i <= orbs.size(); i++) {
    auto [name, ki, en] = orbs.at(i - 1);
    cout << i << ": " << name << " " << ki << " " << en << endl;
  }
}

// Make MPO for the current operator
template <typename Basis1, typename Basis2, typename SiteType>
MPO get_current_mpo(const SiteType& sites, const Basis1& basis1,
                    const Basis2& basis2, int i1, int i2,
                    const ToGlobDict& to_glob) {
  AutoMPO ampo(sites);
  add_CdagC(ampo, basis1, basis2, i1, i2, 1., to_glob);
  auto mpo = toMPO(ampo);
  return mpo;
}

inline Real get_current(const MPO& JMPO, const MPS& psi) {
  auto J = innerC(psi, JMPO, psi);
  return -2. * imag(J);
}

int main(int argc, char* argv[]) {
  string infile = argv[1];
  InputGroup input(infile, "basic");

  auto L_lead = input.getInt("L_lead");
  auto L_device = input.getInt("L_device");
  auto t_lead = input.getReal("t_lead");
  auto t_device = input.getReal("t_device");
  auto t_contactL = input.getReal("t_contactL");
  auto t_contactR = input.getReal("t_contactR");
  auto mu_leadL = input.getReal("mu_leadL");
  auto mu_leadR = input.getReal("mu_leadR");
  auto mu_device = input.getReal("mu_device");
  auto mu_biasL = input.getReal("mu_biasL");
  auto mu_biasS = input.getReal("mu_biasS");
  auto mu_biasR = input.getReal("mu_biasR");
  auto damp_decay_length = input.getInt("damp_decay_length", 0);

  auto dt = input.getReal("dt");
  auto time_steps = input.getInt("time_steps");
  auto NumCenter = input.getInt("NumCenter");
  auto Truncate = input.getYesNo("Truncate");
  auto mixNumCenter = input.getYesNo("mixNumCenter", false);
  auto globExpanNStr = input.getString("globExpanN", "inf");
  int globExpanN;
  if (globExpanNStr == "inf" or globExpanNStr == "Inf" or
      globExpanNStr == "INF")
    globExpanN = std::numeric_limits<int>::max();
  else
    globExpanN = std::stoi(globExpanNStr);
  auto globExpanItv = input.getInt("globExpanItv", 1);
  auto globExpanCutoff = input.getReal("globExpanCutoff", 1e-8);
  auto globExpanKrylovDim = input.getInt("globExpanKrylovDim", 3);
  auto globExpanHpsiCutoff = input.getReal("globExpanHpsiCutoff", 1e-8);
  auto globExpanHpsiMaxDim = input.getInt("globExpanHpsiMaxDim", 300);
  auto globExpanMethod = input.getString("globExpanMethod", "DensityMatrix");

  auto UseSVD = input.getYesNo("UseSVD", true);
  auto SVDmethod =
      input.getString("SVDMethod", "gesdd");  // can be also "ITensor"
  auto WriteDim = input.getInt("WriteDim");

  auto write = input.getYesNo("write", false);
  auto write_dir = input.getString("write_dir", ".");
  auto write_file = input.getString("write_file", "");
  auto read = input.getYesNo("read", false);
  auto read_dir = input.getString("read_dir", ".");
  auto read_file = input.getString("read_file", "");

  auto sweeps = iut::Read_sweeps(infile, "sweeps");

  cout << setprecision(14) << endl;

  MPS psi;
  MPO H;
  // Define
  int step = 1;
  auto sites = Fermion();
  Para para;
  Args args_basis;

  ToGlobDict to_glob;
  ToLocDict to_loc;
  OneParticleBasis leadL, leadR, charge;
  OneParticleBasis scatterer;

  // -- Initialization --
  if (!read) {
    // Factor for exponentially decaying hoppings
    Real damp_fac =
        (damp_decay_length == 0 ? 1. : exp(-1. / damp_decay_length));
    // Create bases for the leads
    cout << "H left lead" << endl;
    leadL =
        OneParticleBasis("L", L_lead, t_lead, mu_leadL, damp_fac, true, true);
    cout << "H right lead" << endl;
    leadR =
        OneParticleBasis("R", L_lead, t_lead, mu_leadR, damp_fac, false, true);
    // Create basis for scatterer
    cout << "H dev" << endl;
    scatterer = OneParticleBasis("S", L_device, t_device, mu_device);

    // Combine and sort all the basis states
    auto info = sort_by_energy(leadL, leadR, scatterer);
    tie(to_glob, to_loc) = make_orb_dicts(info);
    print_orbs(info);

    // SiteSet
    int N = to_glob.size();
    sites = Fermion(N);

    // Make Hamiltonian MPO for time evolution
    para.tcL = t_contactL;
    para.tcR = t_contactR;
    auto ampo =
        get_ampo_tight_binding(leadL, leadR, scatterer, sites, para, to_glob);
    H = toMPO(ampo);
    cout << "MPO dim = " << maxLinkDim(H) << endl;

    // Initialize MPS
    // *********************************
    {
      /*auto para0 = para;
      para0.tcL = 0.;
      para0.tcR = 0.;
      auto leadL0 = OneParticleBasis ("L", L_lead, t_lead, mu_leadL+mu_biasL,
      damp_fac, true, true); auto leadR0 = OneParticleBasis ("R", L_lead,
      t_lead, mu_leadR+mu_biasR, damp_fac, false, true);*/
      psi = get_non_inter_ground_state(
          leadL, leadR, scatterer, sites, mu_leadL + mu_biasL,
          mu_device + mu_biasS, mu_leadR + mu_biasR, to_glob);
      /*auto ampo0 = get_ampo_tight_binding (leadL0, leadR0, scatterer, sites,
      para, to_glob); auto H0 = toMPO (ampo0); auto sweeps0 = iut::Read_sweeps
      (infile, "DMRG_sweeps"); dmrg (psi, H0, sweeps0, {"WriteDim",WriteDim});*/
    }
    // *********************************
    psi.position(1);

    // Check initial energy
    cout << "Initial energy = " << inner(psi, H, psi) << endl;
  } else {
    readAll(read_dir + "/" + read_file, psi, H, para, args_basis, step, to_glob,
            to_loc);
    sites = Fermion(siteInds(psi));
  }
  // -- End of initialization --

  // -- Observer --
  auto obs = TDVPObserver(sites, psi);
  // Current MPO
  auto jmpoL = get_current_mpo(sites, leadL, leadL, -2, -1, to_glob);
  auto jmpoR = get_current_mpo(sites, leadR, leadR, 1, 2, to_glob);

  // -- Time evolution --
  cout << "Start time evolution" << endl;
  cout << sweeps << endl;
  psi.position(1);
  Real en, err;

  Args args_tdvp_expansion = {"Cutoff",      globExpanCutoff,
                              "Method",      "DensityMatrix",
                              "KrylovOrd",   globExpanKrylovDim,
                              "DoNormalize", true,
                              "Quiet",       true};
  Args args_tdvp = {"Quiet",       true,     "NumCenter",    NumCenter,
                    "DoNormalize", true,     "Truncate",     Truncate,
                    "UseSVD",      UseSVD,   "SVDmethod",    SVDmethod,
                    "WriteDim",    WriteDim, "mixNumCenter", mixNumCenter};

  LocalMPO PH(H, args_tdvp);
  while (step <= time_steps) {
    cout << "step = " << step << endl;

    // Subspace expansion
    if (maxLinkDim(psi) < sweeps.mindim(1) or
        (step < globExpanN and (step - 1) % globExpanItv == 0)) {
      timer["glob expan"].start();
      addBasis(psi, H, globExpanHpsiCutoff, globExpanHpsiMaxDim,
               args_tdvp_expansion);
      PH.reset();
      timer["glob expan"].stop();
    }

    // Time evolution
    timer["tdvp"].start();
    TDVPWorker(psi, PH, 1_i * dt, sweeps, obs, args_tdvp);
    timer["tdvp"].stop();
    auto d1 = maxLinkDim(psi);

    // Measure currents by MPO
    timer["current mps"].start();
    auto jL = get_current(jmpoL, psi);
    auto jR = get_current(jmpoR, psi);
    cout << "\tI L/R = " << jL << " " << jR << endl;
    timer["current mps"].stop();

    step++;
    if (write) {
      timer["write"].start();
      writeAll(write_dir + "/" + write_file, psi, H, para, args_basis, step,
               to_glob, to_loc);
      timer["write"].stop();
    }
  }
  timer.print();
  return 0;
}
