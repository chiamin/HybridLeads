#ifndef __INITSTATE_H_CMC__
#define __INITSTATE_H_CMC__
#include "SortBasis.h"

// The charging energy is Ec * (N - Ng)^2
// Find out N that lowest the charging energy in even and odd number of particle
// sectors
tuple<Real, Real, int, int> en_charging_energy(int maxOcc, Real Ec, Real Ng) {
  Real en_even = std::numeric_limits<double>::max(),
       en_odd = std::numeric_limits<double>::max();
  vector<int> ns(1, 0);
  for (int n = 1; n <= maxOcc; n++) {
    ns.push_back(n);
    ns.push_back(-n);
  }
  int n_even, n_odd;
  for (int n : ns) {
    // Compute charging energy
    Real nn = n - Ng;
    Real enC = Ec * nn * nn;
    //
    if (n % 2 == 0 and enC < en_even) {
      en_even = enC;
      n_even = n;
    }
    if (n % 2 == 1 and enC < en_odd) {
      en_odd = enC;
      n_odd = n;
    }
  }
  return {en_even, en_odd, n_even, n_odd};
}

template <typename BasisL, typename BasisR, typename BasisS, typename SiteType,
          typename Para>
MPS get_ground_state_BdG_scatter(const BasisL& leadL, const BasisR& leadR,
                                 const BasisS& scatterer, const SiteType& sites,
                                 Real muL, Real muR, const Para& para,
                                 int maxOcc, const ToGlobDict& to_glob) {
  int N = to_glob.size();
  mycheck(length(sites) == N, "size not match");

  // Get the ground state for SC (even particle number)
  vector<string> state(N + 1);

  // Leads
  auto occ_negative_en_states = [&to_glob, &state](const auto& basis, Real mu) {
    string p = basis.name();
    for (int k = 1; k <= basis.size(); k++) {
      int i = to_glob.at({p, k});
      auto en = basis.en(k);
      if (en < mu)
        state.at(i) = "Occ";
      else
        state.at(i) = "Emp";
    }
  };
  occ_negative_en_states(leadL, muL);
  occ_negative_en_states(leadR, muR);

  // Scatterer
  string sname = scatterer.name();
  for (int k = 1; k <= scatterer.size(); k++) {
    int i = to_glob.at({sname, k});
    state.at(i) = "Emp";
  }

  // Superconducting gap
  Real SC_gap = scatterer.en(1);
  cout << "SC gap = " << SC_gap << endl;

  // Capacity site
  // Find out the charge numbers, which are integers, that lowest the charging
  // energy, for even and odd parities
  auto [enC0, enC1, n_even, n_odd] =
      en_charging_energy(maxOcc, para.Ec, para.Ng);
  // Combine the scatter energies to decide to choose even or odd sector
  Real en0 = enC0, en1 = enC1 + SC_gap;
  cout << "n (even,odd) = " << n_even << ", " << n_odd << endl;
  cout << "E (even,odd) = " << en0 << ", " << en1 << endl;
  if (en1 < en0)  // First excited state in superconductor
  {
    int is1 = to_glob.at({"S", 1});
    state.at(is1) = "Occ";
    cout << "Ground state has odd parity" << endl;
  } else {
    cout << "Ground state has even parity" << endl;
  }
  int n = (en0 <= en1 ? n_even : n_odd);
  cout << "Initial charge = " << n << endl;
  // Set the charge site
  int ic = to_glob.at({"C", 1});
  state.at(ic) = str(n);
  // Set state
  InitState init(sites);
  for (int i = 1; i <= N; i++) init.set(i, state.at(i));

  auto psi = MPS(init);
  /*
      // If Majorana zero mode exists, use the equal-weight superposition of
  occupied and unoccupied state auto const& en = visit (basis::en(1),
  sys.parts().at("S")); if (abs(en) < 1e-14)
      {
          int iglob = sys.to_glob ("S",1);
          auto A = psi(iglob);
          auto links = findInds (A, "Link");
          auto il1 = links(1);
          auto il2 = links(2);
          auto is  = findIndex (A, "Site");
          mycheck (dim(il1) == 1 and dim(il2) == 1, "Link indices dimension
  error"); Real ele = 1./sqrt(2); A.set (il1=1, il2=1, is=1, ele); A.set (il1=1,
  il2=1, is=2, ele); psi.ref(iglob) = A; PrintData(psi(iglob));
      }
  exit(0);*/
  return psi;
}

template <typename BasisL, typename BasisR, typename BasisS, typename SiteType>
MPS get_non_inter_ground_state(const BasisL& leadL, const BasisR& leadR,
                               const BasisS& scatterer, const SiteType& sites,
                               Real muL, Real muS, Real muR,
                               const ToGlobDict& to_glob) {
  int N = to_glob.size();
  mycheck(length(sites) == N, "size not match");

  int Ns = 0, Np = 0;
  Real E = 0.;
  vector<string> state(N + 1, "Emp");

  // Leads and scatterer
  auto occ_negative_en_states = [&to_glob, &state, &E, &Np, &Ns](
                                    const auto& basis, Real mu) {
    string p = basis.name();
    for (int k = 1; k <= basis.size(); k++) {
      int i = to_glob.at({p, k});
      auto en = basis.en(k);
      if (en < mu) {
        state.at(i) = "Occ";
        E += en;
        Np++;
        if (p == "S") Ns++;
      } else {
        state.at(i) = "Emp";
      }
    }
  };
  occ_negative_en_states(leadL, muL);
  occ_negative_en_states(leadR, muR);
  occ_negative_en_states(scatterer, muS);

  InitState init(sites);
  for (int i = 1; i <= N; i++) init.set(i, state.at(i));

  // Print information
  cout << "initial energy = " << E << endl;
  cout << "initial particle number = " << Np << endl;
  return MPS(init);
}

template <typename BasisS>
tuple<MPS, MPS, Real, Real, Real, Real, int> get_scatter_ground_state_SC(
    const BasisS& scatterer, Real mu, Real Delta, const Sweeps& sweeps,
    const ToGlobDict& to_glob, const Args& args) {
  int Ns = scatterer.size();
  SpecialFermion sites(Ns, {args, "in_scatter", true});

  // Find the first site of the scatter
  string sname = scatterer.name();
  int imin = std::numeric_limits<int>::max();
  for (int i = 1; i <= Ns; i++) {
    int j = to_glob.at({sname, i});
    if (imin > j) imin = j;
  }

  int L_offset = imin - 1;
  AutoMPO ampo(sites);
  // Diagonal terms
  for (int i = 1; i <= Ns; i++) {
    int j = to_glob.at({sname, i}) - L_offset;
    auto en = scatterer.en(i);
    ampo += en - mu, "N", j;
  }
  // Superconducting
  for (int i = 1; i < Ns; i++) {
    auto terms =
        quadratic_operator(scatterer, scatterer, i, i + 1, false, false);
    for (auto [c12, k1, dag1, k2, dag2] : terms) {
      int j1 = to_glob.at({sname, k1}) - L_offset;
      int j2 = to_glob.at({sname, k2}) - L_offset;
      if (j1 != j2) {
        auto c = Delta * c12;
        auto cc = iut::conj(c);
        string op1 = (dag1 ? "Cdag" : "C");
        string op2 = (dag2 ? "Cdag" : "C");
        string op1dag = (dag1 ? "C" : "Cdag");
        string op2dag = (dag2 ? "C" : "Cdag");
        ampo += -c, op1, j1, op2, j2;
        ampo += -cc, op1dag, j2, op2dag, j1;
      }
    }
  }
  auto H0 = toMPO(ampo);

  // Solve the ground states in even and odd parities by DMRG
  InitState init(sites);
  auto psi0 = MPS(init);
  init.set(1, "Occ");
  auto psi1 = MPS(init);

  auto en0 = dmrg(psi0, H0, sweeps, {"Quiet", true});
  auto en1 = dmrg(psi1, H0, sweeps, {"Quiet", true});

  cout << setprecision(14);
  AutoMPO Nampo(sites);
  for (int i = 1; i <= length(sites); i++) Nampo += 1.0, "N", i;
  auto Nmpo = toMPO(Nampo);
  Real Np0 = inner(psi0, Nmpo, psi0), Np1 = inner(psi1, Nmpo, psi1);

  return {psi0, psi1, en0, en1, Np0, Np1, L_offset};
}

template <typename BasisL, typename BasisR, typename BasisS, typename BasisC,
          typename SiteType, typename Para>
MPS get_ground_state_SC(const BasisL& leadL, const BasisR& leadR,
                        const BasisS& scatterer, const BasisC& charge,
                        const SiteType& sites, Real muL, Real muS, Real muR,
                        const Para& para, const Sweeps& sweeps,
                        const ToGlobDict& to_glob, const Args& args) {
  int N = to_glob.size();
  mycheck(length(sites) == N, "size not match");

  Real E_lead = 0.;
  int Np_lead = 0;
  vector<string> state(N + 1, "Emp");

  // Leads
  auto occ_neg_en_levels = [&E_lead, &Np_lead, &state, &to_glob](
                               const auto& basis, Real mu) {
    string p = basis.name();
    for (int i = 1; i <= basis.size(); i++) {
      auto en = basis.en(i);
      int j = to_glob.at({p, i});
      if (en < mu) {
        state.at(j) = "Occ";
        E_lead += en - mu;
        Np_lead++;
      } else {
        state.at(j) = "Emp";
      }
    }
  };
  occ_neg_en_levels(leadL, muL);
  occ_neg_en_levels(leadR, muR);
  cout << "lead E = " << E_lead << endl;
  cout << "lead Np = " << Np_lead << endl;

  // Get ground state of the scatterer
  auto [psi0, psi1, enSC0, enSC1, Np0, Np1, L_offset] =
      get_scatter_ground_state_SC(scatterer, muS, para.Delta, sweeps, to_glob,
                                  args);

  // Capacity site
  // Find out the charge numbers, which are integers, that lowest the charging
  // energy, for even and odd parities
  int maxOcc = args.getInt("MaxOcc");
  auto [enC0, enC1, n_even, n_odd] =
      en_charging_energy(maxOcc, para.Ec, para.Ng);
  // Combine the scatter energies to decide to choose even or odd sector
  Real en0 = enC0 + enSC0, en1 = enC1 + enSC1;
  Real en = (en0 < en1 ? en0 : en1);
  int n = (en0 < en1 ? n_even : n_odd);
  Real Np = (en0 < en1 ? Np0 : Np1);
  auto psiS = (en0 < en1 ? psi0 : psi1);
  // Set state
  int ic = to_glob.at({"C", 1});
  state.at(ic) = str(n);
  // Print
  cout << "Init scatter Np (even,odd) = (" << Np0 << "," << Np1 << ") -> " << Np
       << endl;
  cout << "Init scatter E (even,odd) = (" << en0 << "," << en1 << ") -> " << en
       << endl
       << "\tE_C = (" << enC0 << "," << enC1 << ")" << endl
       << "\tE_SC, gap = (" << enSC0 << "," << enSC1 << "), "
       << abs(enSC0 - enSC1) << endl;
  cout << "Init scatter total charge = (" << n_even << "," << n_odd << ") -> "
       << n << endl;

  // Initialize the leads and the charge site
  InitState init(sites);
  for (int i = 1; i <= N; i++) init.set(i, state.at(i));
  auto psi = MPS(init);

  // Replace the tensors in the scatter
  for (int i = 1; i <= length(psiS); i++) {
    int i0 = i + L_offset;
    auto iis = findIndex(psi(i0), "Site");
    auto iis2 = findIndex(psiS(i), "Site");
    Index iil;
    if (i == 1)
      iil = leftLinkIndex(psi, i0);
    else if (i == length(psiS))
      iil = rightLinkIndex(psi, i0);
    psiS.ref(i).replaceInds({iis2}, {iis});
    psi.ref(i0) = psiS(i);
    if (iil) {
      mycheck(dim(iil) == 1, "not dummy index");
      psi.ref(i0) *= setElt(iil = 1);
    }
    state.at(i0) = "*";
  }
  psi.position(1);
  psi.normalize();

  return psi;
}

#endif
