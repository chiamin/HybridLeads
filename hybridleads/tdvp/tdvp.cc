#include "Entanglement.h"
#include "FixedPointTensor.h"
#include "GenMPO.h"
#include "GeneralUtility.h"
#include "GlobalIndices.h"
#include "IUtility.h"
#include "MPSUtility.h"
#include "MyLocalmpo.h"
#include "MyObserver.h"
#include "ReadInput.h"
#include "TDVPWorker.h"
#include "iTDVP.h"
#include "itensor/all.h"
#include "uGauge.h"
using namespace itensor;
using namespace std;

vector<Index> get_site_inds(const MPS& psi) {
  int N = length(psi);
  vector<Index> sites(N);
  for (int i = 1; i <= N; i++) {
    auto ii = findIndex(psi(i), "Site");
    sites.at(i - 1) = ii;
  }
  return sites;
}

Fermion get_SiteSet(const MPS& mps) {
  auto inds = get_site_inds(mps);
  return Fermion(inds);
}

ITensor get_W(const MPO& H, int i) {
  auto W = H(i);
  auto iWl = leftLinkIndex(H, i);
  auto iWr = rightLinkIndex(H, i);
  auto is = findIndex(W, "Site,0");

  auto is2 = global::IS.is();
  auto iwl2 = global::IS.iwl();
  auto iwr2 = global::IS.iwr();
  auto is2pr = prime(is2);
  W.replaceInds({is, prime(is), iWl, iWr}, {is2, is2pr, iwl2, iwr2});
  return W;
}

inline ITensor to_itdvp_AR(const ITensor& AR, const Index& iAl,
                           const Index& iAr) {
  auto ir = global::IS.ir();
  auto is = global::IS.is();
  auto iAs = findIndex(AR, "Site");
  auto ARre = replaceInds(AR, {iAl, iAr, iAs}, {prime(ir, 2), ir, is});
  global::IS.check("AR", ARre);
  return ARre;
}

inline ITensor to_itdvp_AL(const ITensor& AL, const Index& iAl,
                           const Index& iAr) {
  auto il = global::IS.il();
  auto is = global::IS.is();
  auto iAs = findIndex(AL, "Site");
  auto ALre = replaceInds(AL, {iAl, iAr, iAs}, {il, prime(il, 2), is});
  global::IS.check("AL", ALre);
  return ALre;
}

ITensor get_L(const ITensor& AR, const Index& iAl, const Index& iAr,
              Real crit = 1e-15) {
  auto il = global::IS.il();
  auto ir = global::IS.ir();
  auto is = global::IS.is();
  auto iAs = findIndex(AR, "Site");
  auto ARt = replaceInds(AR, {iAl, iAr, iAs}, {il, prime(il, 2), is});

  auto [AL, C] = Orthogonalize(ARt, crit);
  C.replaceInds({il}, {ir});

  auto Cdag = dag(C);
  Cdag.prime(ir);

  auto L = C * Cdag;
  global::IS.check("L", L);
  assert(check_leading_eigen(AR, AR, L));
  return L;
}

ITensor get_R(const ITensor& AL, const Index& iAl, const Index& iAr,
              Real crit = 1e-15) {
  auto il = global::IS.il();
  auto is = global::IS.is();
  auto iAs = findIndex(AL, "Site");
  auto ALt = replaceInds(AL, {iAl, iAr, iAs}, {prime(il, 2), il, is});

  auto [AR, C] = Orthogonalize(ALt, crit);
  C.mapPrime(1, 2);

  auto Cdag = dag(C);
  Cdag.mapPrime(0, 2);

  auto R = C * Cdag;
  global::IS.check("R", R);
  assert(check_leading_eigen(AL, AL, R));
  return R;
}

void expandL(MPS& psi, MPO& H, unique_ptr<MyLocalMPO>& PH, int n, int NumCenter,
             const ITensor& AL, const Index& iALl, const Index& iALr) {
  int oc = orthoCenter(psi);
  PH->position(oc, psi);

  int N = length(psi);
  int N2 = N + n;

  auto WL = H(1);

  auto iALs = findIndex(AL, "Site");
  auto iWLl = iut::leftIndex(H, 1);
  auto iWLr = iut::rightIndex(H, 1);
  auto iWLs = findIndex(H(1), "Site,0");

  // Set MPS and MPO
  MPS psi2(N2);
  MPO H2(N2);
  // Insert original tensors
  for (int i = 1; i <= N; i++) {
    psi2.ref(n + i) = psi(i);
    H2.ref(n + i) = H(i);
  }
  // Insert new tensor to the left
  auto il0 = iut::leftIndex(psi, 1);
  auto iHl0 = iut::leftIndex(H, 1);
  auto il = sim(il0);
  auto iHl = sim(iHl0);
  psi2.ref(n + 1).replaceInds({il0}, {il});
  H2.ref(n + 1).replaceInds({iHl0}, {iHl});
  for (int i = n; i >= 1; i--) {
    auto is2 = sim(iALs);
    // Replace the site and the right indices
    psi2.ref(i) = replaceInds(AL, {iALs, iALr}, {is2, dag(il)});
    H2.ref(i) =
        replaceInds(WL, {iWLs, prime(iWLs), iWLr}, {is2, prime(is2), dag(iHl)});
    il = sim(iALl);
    iHl = sim(iWLl);
    if (i != 1) {
      psi2.ref(i).replaceInds({iALl}, {il});
      H2.ref(i).replaceInds({iWLl}, {iHl});
    }
  }
  psi2.rightLim(n + 2);
  psi2.position(n + 1);
  psi = psi2;
  H = H2;

  // Set PH
  Args args = {"NumCenter", NumCenter};
  auto PH2 = make_unique<MyLocalMPO>(H, PH->L(1), PH->R(N), args);
  for (int i = N - 1; i >= 1; i--) PH2->setR(i + n, PH->R(i));
  PH2->setRHlim(n + 1 + NumCenter);
  PH2->position(n + 1, psi);
  PH = move(PH2);
}

void expandR(MPS& psi, MPO& H, unique_ptr<MyLocalMPO>& PH, int n, int NumCenter,
             const ITensor& AR, const Index& iARl, const Index& iARr) {
  int oc = orthoCenter(psi);
  PH->position(oc, psi);

  int N = length(psi);
  int N2 = N + n;

  auto WR = H(N);

  auto iARs = findIndex(AR, "Site");
  auto iWRl = iut::leftIndex(H, N);
  auto iWRr = iut::rightIndex(H, N);
  auto iWRs = findIndex(H(N), "Site,0");

  // Set MPS and MPO
  MPS psi2(N2);
  MPO H2(N2);
  // Insert original tensors
  for (int i = 1; i <= N; i++) {
    psi2.ref(i) = psi(i);
    H2.ref(i) = H(i);
  }
  // Insert new tensors to the right
  auto ir0 = iut::rightIndex(psi, N);
  auto iHr0 = iut::rightIndex(H, N);
  auto ir = sim(ir0);
  auto iHr = sim(iHr0);
  psi2.ref(N).replaceInds({ir0}, {ir});
  H2.ref(N).replaceInds({iHr0}, {iHr});
  for (int i = N + 1; i <= N2; i++) {
    auto is2 = sim(iARs);
    // Replace the site and the left indices
    psi2.ref(i) = replaceInds(AR, {iARs, iARl}, {is2, dag(ir)});
    H2.ref(i) =
        replaceInds(WR, {iWRs, prime(iWRs), iWRl}, {is2, prime(is2), dag(iHr)});
    ir = sim(iARr);
    iHr = sim(iWRr);
    if (i != N2) {
      psi2.ref(i).replaceInds({iARr}, {ir});
      H2.ref(i).replaceInds({iWRr}, {iHr});
    }
  }
  psi2.leftLim(N - 1);
  psi2.position(N);
  psi = psi2;
  H = H2;

  // Set PH
  Args args = {"NumCenter", NumCenter};
  auto PH2 = make_unique<MyLocalMPO>(H, PH->L(1), PH->R(N), args);
  for (int i = 2; i <= N; i++) PH2->setL(i, PH->L(i));
  PH2->setLHlim(N - NumCenter);
  PH2->position(N, psi2);
  PH = move(PH2);
}

void get_init(const string& infile, MPS& psi, MPO& H, Fermion& sites, Index& il,
              Index& ir, ITensor& WL, ITensor& WR, ITensor& AL_left,
              ITensor& AR_left, ITensor& AC_left, ITensor& C_left,
              ITensor& La_left, ITensor& Ra_left, ITensor& LW_left,
              ITensor& RW_left, ITensor& AL_right, ITensor& AR_right,
              ITensor& AC_right, ITensor& C_right, ITensor& La_right,
              ITensor& Ra_right, ITensor& LW_right, ITensor& RW_right,
              int& idev_first, int& idev_last) {
  InputGroup input(infile, "basic");

  auto L_window = input.getInt("L_window");
  auto L_device = input.getInt("L_device");
  auto t_lead = input.getReal("t_lead");
  auto t_device = input.getReal("t_device");
  auto t_contact = input.getReal("t_contact");
  auto mu_leadL = input.getReal("mu_leadL");
  auto mu_leadR = input.getReal("mu_leadR");
  auto mu_device = input.getReal("mu_device");
  auto V_lead = input.getReal("V_lead");
  auto V_device = input.getReal("V_device");
  auto V_contact = input.getReal("V_contact");

  auto psi_dir = input.getString("psi_dir");
  auto psi_file = input.getString("psi_file");
  auto itdvp_dir = input.getString("itdvp_dir");
  auto ErrGoal = input.getReal("ErrGoal");
  auto MaxIter_LRW = input.getInt("MaxIter_LRW");

  // Read the iTDVP tensors
  ITensor AL, AR, AC, C, LW, RW, La, Ra;
  readFromFile(itdvp_dir + "/AL.itensor", AL);
  readFromFile(itdvp_dir + "/AR.itensor", AR);
  readFromFile(itdvp_dir + "/AC.itensor", AC);
  readFromFile(itdvp_dir + "/C.itensor", C);
  readFromFile(itdvp_dir + "/LW.itensor", LW);
  readFromFile(itdvp_dir + "/RW.itensor", RW);
  readFromFile(itdvp_dir + "/La.itensor", La);
  readFromFile(itdvp_dir + "/Ra.itensor", Ra);
  global::IS.read(itdvp_dir + "/global.inds");
  il = global::IS.il();
  ir = global::IS.ir();
  int dim = ir.dim();
  // Read MPS in the window
  readFromFile(psi_dir + "/" + psi_file, psi);
  int N = length(psi);
  // Site indices
  auto site_inds = get_site_inds(psi);
  sites = Fermion(site_inds);
  // Hamiltonian
  AutoMPO ampo;
  tie(ampo, idev_first, idev_last) = t_mu_V_ampo(
      sites, L_device, t_lead, t_lead, t_device, t_contact, t_contact, mu_leadL,
      mu_leadR, mu_device, V_lead, V_lead, V_device, V_contact, V_contact);
  auto locamu = read_bracket_values<int, Real>(infile, "local_mu", 1);
  auto localtV =
      read_bracket_values<int, int, Real, Real>(infile, "local_tV", 1);
  ampo_add_mu(ampo, locamu, true);
  ampo_add_tV(ampo, localtV, true);
  H = toMPO(ampo);
  to_inf_mpo(H, global::IS.iwl(), global::IS.iwr());
  // Boundary tensors
  cout << "Compute boundary tensors" << endl;
  WL = get_W(H, 2);
  WR = get_W(H, N - 1);
  Args args_itdvp = {"ErrGoal=", ErrGoal, "MaxIter", MaxIter_LRW};
  auto L = get_LR<LEFT>(C);
  auto R = get_LR<RIGHT>(C);
  Real enL, enR;
  tie(LW, enL) = get_LRW<LEFT>(AL, WL, R, La, args_itdvp);
  tie(RW, enR) = get_LRW<RIGHT>(AR, WR, L, Ra, args_itdvp);
  // Check
  mycheck(commonIndex(LW, H(1)), "LW and H(1) has no common Index");
  mycheck(commonIndex(RW, H(N)), "RW and H(N) has no common Index");
  mycheck(commonIndex(LW, psi(1)), "LW and psi(1) has no common Index");
  mycheck(commonIndex(RW, psi(N)), "RW and psi(N) has no common Index");
  // Left boundaries
  AL_left = AL, AR_left = AR, AC_left = AC, C_left = C, La_left = La,
  Ra_left = Ra, LW_left = LW, RW_left = RW;
  // Right boundaries
  AL_right = AL, AR_right = AR, AC_right = AC, C_right = C, La_right = La,
  Ra_right = Ra, LW_right = LW, RW_right = RW;
}

void writeAll(const string& filename, const MPS& psi, const MPO& H,
              const Index& il, const Index& ir, const ITensor& WL,
              const ITensor& WR, const ITensor& AL_left, const ITensor& AR_left,
              const ITensor& AC_left, const ITensor& C_left,
              const ITensor& La_left, const ITensor& Ra_left,
              const ITensor& LW_left, const ITensor& RW_left,
              const ITensor& AL_right, const ITensor& AR_right,
              const ITensor& AC_right, const ITensor& C_right,
              const ITensor& La_right, const ITensor& Ra_right,
              const ITensor& LW_right, const ITensor& RW_right, int step,
              int idev_first, int idev_last, bool expand, bool expand_next) {
  ofstream ofs(filename);
  write(ofs, psi);
  write(ofs, H);
  write(ofs, il);
  write(ofs, ir);
  write(ofs, WL);
  write(ofs, WR);
  write(ofs, AL_left);
  write(ofs, AR_left);
  write(ofs, AC_left);
  write(ofs, C_left);
  write(ofs, La_left);
  write(ofs, Ra_left);
  write(ofs, LW_left);
  write(ofs, RW_left);
  write(ofs, AL_right);
  write(ofs, AR_right);
  write(ofs, AC_right);
  write(ofs, C_right);
  write(ofs, La_right);
  write(ofs, Ra_right);
  write(ofs, LW_right);
  write(ofs, RW_right);
  write(ofs, step);
  write(ofs, idev_first);
  write(ofs, idev_last);
  write(ofs, expand);
  write(ofs, expand_next);
  global::IS.write(ofs);
}

void readAll(const string& filename, MPS& psi, MPO& H, Index& il, Index& ir,
             ITensor& WL, ITensor& WR, ITensor& AL_left, ITensor& AR_left,
             ITensor& AC_left, ITensor& C_left, ITensor& La_left,
             ITensor& Ra_left, ITensor& LW_left, ITensor& RW_left,
             ITensor& AL_right, ITensor& AR_right, ITensor& AC_right,
             ITensor& C_right, ITensor& La_right, ITensor& Ra_right,
             ITensor& LW_right, ITensor& RW_right, int& step, int& idev_first,
             int& idev_last, bool& expand, bool& expand_next) {
  ifstream ifs = open_file(filename);
  read(ifs, psi);
  read(ifs, H);
  read(ifs, il);
  read(ifs, ir);
  read(ifs, WL);
  read(ifs, WR);
  read(ifs, AL_left);
  read(ifs, AR_left);
  read(ifs, AC_left);
  read(ifs, C_left);
  read(ifs, La_left);
  read(ifs, Ra_left);
  read(ifs, LW_left);
  read(ifs, RW_left);
  read(ifs, AL_right);
  read(ifs, AR_right);
  read(ifs, AC_right);
  read(ifs, C_right);
  read(ifs, La_right);
  read(ifs, Ra_right);
  read(ifs, LW_right);
  read(ifs, RW_right);
  read(ifs, step);
  read(ifs, idev_first);
  read(ifs, idev_last);
  read(ifs, expand);
  read(ifs, expand_next);
  global::IS.read(ifs);
}

int main(int argc, char* argv[]) {
  string infile = argv[1];
  InputGroup input(infile, "basic");

  auto dt = input.getReal("dt");
  auto time_steps = input.getInt("time_steps");
  auto NumCenter = input.getInt("NumCenter");
  auto ConserveQNs = input.getYesNo("ConserveQNs", false);
  auto expandN = input.getInt("expandN", 0);
  auto expand_checkN = input.getInt("expand_checkN", 5);
  auto expandS_crit = input.getReal("expandS_crit", 1e10);
  auto max_window = input.getInt("max_window", 10000);
  auto sweeps = iut::Read_sweeps(infile);

  auto ErrGoal = input.getReal("ErrGoal");
  auto MaxIter_LRW = input.getInt("MaxIter_LRW");
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

  auto out_dir = input.getString("outdir", ".");
  if (write_dir == "." && out_dir != ".") write_dir = out_dir;

  // Declare variables
  MPS psi;
  MPO H;
  Fermion sites;
  Index il, ir;
  ITensor WL, WR, AL_left, AR_left, AC_left, C_left, La_left, Ra_left, LW_left,
      RW_left, AL_right, AR_right, AC_right, C_right, La_right, Ra_right,
      LW_right, RW_right;
  bool expand = false, expand_next = false;
  int step = 1;
  int idev_first, idev_last;

  // Initialize variables
  if (!read) {
    get_init(infile, psi, H, sites, il, ir, WL, WR, AL_left, AR_left, AC_left,
             C_left, La_left, Ra_left, LW_left, RW_left, AL_right, AR_right,
             AC_right, C_right, La_right, Ra_right, LW_right, RW_right,
             idev_first, idev_last);
  }
  // Read variables
  else {
    readAll(read_dir + "/" + read_file, psi, H, il, ir, WL, WR, AL_left,
            AR_left, AC_left, C_left, La_left, Ra_left, LW_left, RW_left,
            AL_right, AR_right, AC_right, C_right, La_right, Ra_right, LW_right,
            RW_right, step, idev_first, idev_last, expand, expand_next);
    auto site_inds = get_site_inds(psi);
    sites = Fermion(site_inds);
  }

  // Args parameters
  Args args_itdvp = {"ErrGoal=", ErrGoal, "MaxIter", MaxIter_LRW};
  Args args_obs = {"ConserveQNs", ConserveQNs};
  Args args_tdvp = {"Quiet",       true,      "NumCenter", NumCenter,
                    "DoNormalize", true,      "UseSVD",    UseSVD,
                    "SVDmethod",   SVDmethod, "WriteDim",  WriteDim};

  // Observer
  auto obs = make_unique<MyObserver>(sites, psi, args_obs);

  // Effective Hamiltonian
  auto PH = make_unique<MyLocalMPO>(H, LW_left, RW_right, args_tdvp);

  // Time evolution
  cout << "Start time evolution" << endl;
  cout << sweeps << endl;
  psi.position(1);
  Real en, err;
  int N = length(psi);

  for (int i = 0; i < time_steps; i++) {
    cout << "step = " << step++ << endl;

    // Extend left edge
    if (expand) {
      // Expand
      expandL(psi, H, PH, expandN, NumCenter, AL_left, il, prime(il, 2));
      idev_first += expandN;
      idev_last += expandN;
      cout << "expand left " << expandN << endl;
      // Observer
      sites = get_SiteSet(psi);
      obs = make_unique<MyObserver>(sites, psi, args_obs);
      N = length(psi);
      psi.position(1);
    }
    cout << "device site = " << idev_first << " " << idev_last << endl;

    // Evolve left edge
    tie(en, err, LW_left, RW_left) =
        itdvp(WL, AL_left, AR_left, AC_left, C_left, La_left, Ra_left, 1_i * dt,
              args_itdvp);
    PH->L(1, LW_left);

    // From left to right
    TDVPWorker<Fromleft>(psi, *PH, 1_i * dt, sweeps, *obs, args_tdvp);

    if (expandN != 0) {
      auto const& specR = obs->spec(N - expand_checkN);
      auto SR_sys = iut::EntangEntropy(specR);
      auto SR = iut::EntangEntropy_singular(C_right);
      expand_next = (abs(SR_sys - SR) > expandS_crit);
    }

    // Extend right boundary
    if (expand) {
      // Expand
      expandR(psi, H, PH, expandN, NumCenter, AR_right, prime(ir, 2), ir);
      // Observer
      sites = get_SiteSet(psi);
      obs = make_unique<MyObserver>(sites, psi, args_obs);
      N = length(psi);
      psi.position(N);
    }

    // Evolve right boundary
    tie(en, err, LW_right, RW_right) =
        itdvp(WR, AL_right, AR_right, AC_right, C_right, La_right, Ra_right,
              1_i * dt, args_itdvp);
    PH->R(N, RW_right);

    // From right to left
    TDVPWorker<Fromright>(psi, *PH, 1_i * dt, sweeps, *obs, args_tdvp);

    if (expandN != 0) {
      auto const& specL = obs->spec(expand_checkN);
      auto SL_sys = iut::EntangEntropy(specL);
      auto SL = iut::EntangEntropy_singular(C_left);
      if (abs(SL_sys - SL) > expandS_crit) expand_next = true;
    }

    // if (maxLinkDim(psi) >= sweeps.maxdim(1))
    // args_tdvp.add ("NumCenter",1);
    expand = expand_next;
    if (N >= max_window) {
      expand = false;
      expandN = 0;
    }

    if (write) {
      writeAll(write_dir + "/" + write_file, psi, H, il, ir, WL, WR, AL_left,
               AR_left, AC_left, C_left, La_left, Ra_left, LW_left, RW_left,
               AL_right, AR_right, AC_right, C_right, La_right, Ra_right,
               LW_right, RW_right, step, idev_first, idev_last, expand,
               expand_next);
    }
  }

  return 0;
}
