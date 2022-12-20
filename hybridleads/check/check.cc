#include "Corr.h"
#include "IUtility.h"
#include "MPSUtility.h"
#include "MyObserver.h"
#include "ReadInput.h"
#include "TDVPObserver.h"
#include "itensor/all.h"
#include "tdvp.h"
// #include
// "/home/chiamin/Project/2021/topological_qubit/code/mixedbasis/git/SingleParticle.h"
using namespace itensor;
using namespace std;

Matrix Hamilt_k(int L, Real t, Real mu, Real damp_fac = 1.,
                bool damp_left = true, bool verbose = false) {
  cout << "L = " << L << endl;
  Matrix H(L, L);
  for (int i = 0; i < L; i++) {
    H(i, i) = -mu;
    if (i != L - 1) {
      int damp_dist = (damp_left ? L - 2 - i : i);
      Real ti = t * pow(damp_fac, damp_dist);
      H(i, i + 1) = -ti;
      H(i + 1, i) = -ti;
      if (verbose) cout << "Hk, t " << i << " = " << ti << endl;
    }
  }
  return H;
}

int get_Npar(const Matrix& Hk, Real mu) {
  Matrix U;
  Vector ens;
  diagHermitian(Hk, U, ens);

  int N = 0;
  for (int i = 0; i < ens.size(); i++) {
    auto en = ens(i);
    if (en < mu) N += 1;
  }
  return N;
}

vector<string> n_product_state(int L, int Np) {
  vector<string> nstr(L);
  for (int i = 0; i < L; i++) {
    string state;
    if (i % 2 == 0 and Np > 0) {
      nstr.at(i) = "Occ";
      Np--;
    } else
      nstr.at(i) = "Emp";
  }
  if (Np > 0) {
    for (int i = 0; i < L; i++)
      if (Np > 0 and nstr.at(i) == "Emp") {
        nstr.at(i) = "Occ";
        Np--;
      }
  }
  return nstr;
}

void set_initstate(InitState& init, int Np, int i1, int i2) {
  int L = i2 - i1 + 1;
  if (Np < 0 or Np > L) {
    cout << "Error: not valid Np: " << Np << endl;
    cout << "L = " << L << endl;
    throw;
  }
  auto nstr = n_product_state(L, Np);
  for (int i = i1; i <= i2; i++) init.set(i, nstr.at(i - i1));
}

MPS make_initstate(const Fermion& sites, int Np) {
  int N = length(sites);
  int Npar = Np;
  InitState init(sites);
  for (int i = 1; i <= N; i++) {
    string state;
    if (i % 2 == 0 && Np-- > 0)
      state = "Occ";
    else
      state = "Emp";
    init.set(i, state);
    cout << i << ": " << state << endl;
  }
  if (Np > 0) {
    for (int i = 1; i <= N; i += 2)
      if (Np-- > 0) init.set(i, "Occ");
  }
  auto psi = MPS(init);

  auto Nmpo = Make_NMPO(sites);
  auto Ntot = inner(psi, Nmpo, psi);
  if (Ntot != Npar) {
    cout << "particle number not match:" << Ntot << " " << Npar << endl;
    throw;
  }
  cout << "Ntot = " << Ntot << endl;
  return psi;
}

AutoMPO set_H(const SiteSet& sites, int L_lead, int L_device, Real t, Real mu,
              Real V, bool periodic = false) {
  int N = length(sites);
  AutoMPO ampo(sites);
  for (int i = 1; i <= N; ++i) {
    if (i != N) {
      ampo += -t, "Cdag", i, "C", i + 1;
      ampo += -t, "Cdag", i + 1, "C", i;
      cout << "H t " << i << " " << t << endl;
    }
    if (mu != 0.) ampo += -mu, "N", i;
  }
  if (periodic) {
    ampo += -t, "Cdag", N, "C", 1;
    ampo += -t, "Cdag", 1, "C", N;
    cout << "H t " << N << " " << t << endl;
  }

  int Ri = L_lead + L_device;
  for (int i = L_lead + 1; i <= Ri; i++)
    for (int j = L_lead + 1; j < i; j++) {
      ampo += V, "N", i, "N", j;
    }
  return ampo;
}

AutoMPO set_H(const SiteSet& sites, int L_lead, int L_device, Real t_lead,
              Real t_device, Real t_contact, Real muL, Real muS, Real muR,
              Real V, bool periodic = false) {
  int N = length(sites);
  AutoMPO ampo(sites);
  for (int i = 1; i <= N; ++i) {
    if (i != N) {
      Real t;
      if (i == L_lead or i == L_lead + L_device)
        t = t_contact;
      else if (i < L_lead or i > L_lead + L_device)
        t = t_lead;
      else
        t = t_device;
      ampo += -t, "Cdag", i, "C", i + 1;
      ampo += -t, "Cdag", i + 1, "C", i;
      cout << "H t " << i << " " << t << endl;
    }
    {
      Real mu;
      if (i <= L_lead)
        mu = muL;
      else if (i <= L_lead + L_device)
        mu = muS;
      else
        mu = muR;
      ampo += -mu, "N", i;
      cout << "H mu " << i << " " << mu << endl;
    }
  }
  if (periodic) {
    ampo += -t_lead, "Cdag", N, "C", 1;
    ampo += -t_lead, "Cdag", 1, "C", N;
    cout << "H t " << N << " " << t_lead << endl;
  }

  int Ri = L_lead + L_device;
  for (int i = L_lead + 1; i <= Ri; i++)
    for (int j = L_lead + 1; j < i; j++) {
      ampo += V, "N", i, "N", j;
      cout << "H V " << V << " " << i << " " << j << endl;
    }
  return ampo;
}

AutoMPO set_H_SC(const SiteSet& sites, int L_lead, int L_device, Real t_lead,
                 Real t_device, Real t_contact, Real muL, Real muS, Real muR,
                 Real V, Real Delta) {
  int N = length(sites);
  AutoMPO ampo(sites);
  for (int i = 1; i <= N; ++i) {
    if (i != N) {
      Real t;
      if (i == L_lead or i == L_lead + L_device)
        t = t_contact;
      else if (i < L_lead or i > L_lead + L_device)
        t = t_lead;
      else
        t = t_device;
      ampo += -t, "Cdag", i, "C", i + 1;
      ampo += -t, "Cdag", i + 1, "C", i;
      cout << "H t " << i << " " << t << endl;
    }
    {
      Real mu;
      if (i <= L_lead)
        mu = muL;
      else if (i <= L_lead + L_device)
        mu = muS;
      else
        mu = muR;
      ampo += -mu, "N", i;
      cout << "H mu " << i << " " << mu << endl;
    }
  }

  int Ri = L_lead + L_device;
  if (V != 0.) {
    for (int i = L_lead + 1; i <= Ri; i++)
      for (int j = L_lead + 1; j < i; j++) {
        ampo += V, "N", i, "N", j;
        cout << "H V " << V << " " << i << " " << j << endl;
      }
  }
  // SC
  if (Delta != 0.) {
    for (int i = L_lead + 1; i < Ri; i++) {
      ampo += -Delta, "C", i, "C", i + 1;
      ampo += -Delta, "Cdag", i + 1, "Cdag", i;
    }
  }
  return ampo;
}

template <typename MPSType>
ITensor print_wf(const MPSType& psi) {
  ITensor pp(1.);
  vector<Index> iis;
  for (int i = 1; i <= length(psi); i++) {
    pp *= psi(i);
    auto is = findIndex(psi(i), "Site,0");
    iis.push_back(is);
    if constexpr (is_same_v<MPO, MPSType>) iis.push_back(prime(is));
  }
  pp.permute(iis);
  PrintData(pp);
  return pp;
}

Real den(const SiteSet& sites, const MPS& psi, int i) {
  AutoMPO ampo(sites);
  ampo += 1., "N", i;
  auto n = toMPO(ampo);
  return inner(psi, n, psi);
}

Real den(const SiteSet& sites, const MPS& psi, int i1, int i2) {
  Real n = 0.;
  for (int i = i1; i <= i2; i++) n += den(sites, psi, i);
  return n;
}

int main(int argc, char* argv[]) {
  string infile = argv[1];
  InputGroup input(infile, "basic");

  auto L_lead = input.getInt("L_lead");
  auto L_device = input.getInt("L_device");
  int L = 2 * L_lead + L_device;
  auto t_lead = input.getReal("t_lead");
  auto t_device = input.getReal("t_device");
  auto t_contact = input.getReal("t_contact");
  auto mu_leadL = input.getReal("mu_leadL");
  auto mu_leadR = input.getReal("mu_leadR");
  auto mu_device = input.getReal("mu_device");
  auto V = input.getReal("V");
  auto Delta = input.getReal("Delta");

  auto do_write = input.getYesNo("write_to_file");
  auto out_dir = input.getString("outdir", ".");
  auto out_minm = input.getInt("out_minm", 0);
  auto ConserveQNs = input.getYesNo("ConserveQNs", true);
  auto ConserveNf = input.getYesNo("ConserveNf", true);
  auto WriteDim = input.getInt("WriteDim", -1);
  auto sweeps = iutility::Read_sweeps(infile, "DMRG");

  auto mu_biasL = input.getReal("mu_biasL");
  auto mu_biasS = input.getReal("mu_biasS");
  auto mu_biasR = input.getReal("mu_biasR");
  auto dt = input.getReal("dt");
  auto time_steps = input.getInt("time_steps");
  auto NumCenter = input.getInt("NumCenter");
  auto Truncate = input.getYesNo("Truncate");
  auto sweepst = iutility::Read_sweeps(infile, "TDVP");

  cout << "device site = " << L_lead + 1 << " " << (L_lead + L_device) << endl;
  cout << setprecision(14);
  // Site set
  using SitesType = Fermion;
  auto sites =
      SitesType(L, {"ConserveQNs", ConserveQNs, "ConserveNf", ConserveNf});

  // ----------- Getting initial state ----------------
  // Compute particle numbers
  auto Hk_L = Hamilt_k(L_lead, t_lead, mu_leadL);
  auto Hk_R = Hamilt_k(L_lead, t_lead, mu_leadR);
  auto Hk_S = Hamilt_k(L_device, t_device, mu_device);

  int Np_L = get_Npar(Hk_L, mu_leadL + mu_biasL);
  int Np_R = get_Npar(Hk_R, mu_leadR + mu_biasR);
  int Np_S = get_Npar(Hk_S, mu_device + mu_biasS);
  cout << "Np L,R,S = " << Np_L << " " << Np_R << " " << Np_S << endl;

  // Initialze MPS for DMRG
  InitState init(sites);
  set_initstate(init, Np_L, 1, L_lead);
  set_initstate(init, Np_S, L_lead + 1, L_lead + L_device);
  set_initstate(init, Np_R, L_lead + L_device + 1, L);
  auto psi = MPS(init);
  psi.position(1);

  cout << "Np_L = " << den(sites, psi, 1, L_lead) << endl;
  cout << "Np_S = " << den(sites, psi, L_lead + 1, L_lead + L_device) << endl;
  cout << "Np_R = " << den(sites, psi, L_lead + L_device + 1, L) << endl;
  cout << "Ntot = " << den(sites, psi, 1, L) << endl;

  // Make initial-state Hamiltonian MPO
  auto ampo =
      set_H(sites, L_lead, L_device, t_lead, t_device, 0., mu_leadL + mu_biasL,
            mu_device + mu_biasS, mu_leadR + mu_biasR, V);
  auto H = toMPO(ampo);
  cout << "MPO dim = " << maxLinkDim(H) << endl;

  // DMRG to find the initial state
  Real en0 = dmrg(psi, H, sweeps);  // Get ground state
  cout << "Initial energy = " << en0 << endl;

  cout << "Np_L = " << den(sites, psi, 1, L_lead) << endl;
  cout << "Np_S = " << den(sites, psi, L_lead + 1, L_lead + L_device) << endl;
  cout << "Np_R = " << den(sites, psi, L_lead + L_device + 1, L) << endl;
  cout << "Ntot = " << den(sites, psi, 1, L) << endl;
  // --------------------------------------------------------------------------

  // ------------------ Time evolution ---------------------
  // Define time-evolution Hamiltonian MPO
  ampo = set_H(sites, L_lead, L_device, t_lead, t_device, t_contact, mu_leadL,
               mu_device, mu_leadR, V);
  auto Ht = toMPO(ampo);

  // Define current operator MPO
  int N = length(psi);
  vector<int> spec_links = {L_lead, L_lead + L_device};
  vector<MPO> JMPOs(N);
  for (int i = 1; i < L; i++) {
    AutoMPO ampoj(sites);
    ampoj += 1., "Cdag", i, "C", i + 1;
    JMPOs.at(i) = toMPO(ampoj);
  }

  // Time evolution by using TDVP
  Args args_tdvp = {"Quiet",       true, "NumCenter", NumCenter,
                    "DoNormalize", true, "Truncate",  Truncate};
  auto obs = TDVPObserver(sites, psi);
  for (int step = 1; step <= time_steps; step++) {
    cout << "step = " << step << endl;
    tdvp(psi, Ht, 1_i * dt, sweepst, obs, args_tdvp);

    // Measure currents by MPO
    for (int j = 1; j < L; j++) {
      auto Jtmp = innerC(psi, JMPOs.at(j), psi);
      auto J = -2. * imag(Jtmp);
      cout << "\t*current spec " << j << " " << j + 1 << " = " << J << endl;
    }
  }
  return 0;
}
