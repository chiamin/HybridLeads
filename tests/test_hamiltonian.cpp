#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../kbasis/Hamiltonian.h"
#include "../kbasis/OneParticleBasis.h"
#include "../kbasis/SortBasis.h"
#include "itensor/all.h"

using namespace itensor;
using namespace Catch;

struct Para {
  Real tcL = 0., tcR = 0.;
};

TEST_CASE("Check adding CdagC into AutoMPO", "[add_CdagC]") {}

TEST_CASE("Check AutoMPO of tight-binding model", "[get_ampo_tight_binding]") {
  int N = 12;
  auto t = 0.5;
  auto mu = 0.1;
  auto sites = Fermion(N);
  ToGlobDict to_glob;
  ToLocDict to_loc;
  Para para;
  para.tcL = t;
  para.tcR = t;

  // Construct AutoMPO in a hybrid basis jointed by 3 bases
  OneParticleBasis left_basis("left_sys", N / 3, t, mu);
  OneParticleBasis middle_basis("middle_sys", N / 3, t, mu);
  OneParticleBasis right_basis("right_sys", N / 3, t, mu);
  auto info = sort_by_energy(left_basis, right_basis, middle_basis);
  std::tie(to_glob, to_loc) = make_orb_dicts(info);
  auto ampo = get_ampo_tight_binding(left_basis, right_basis, middle_basis,
                                     sites, para, to_glob);
  auto H = toMPO(ampo);

  // Construct AutoMPO in real-space basis
  auto expected_ampo = AutoMPO(sites);
  for (int i = 1; i <= N; ++i) {
    expected_ampo += -mu, "N", i;
  }
  for (int i = 1; i < N; ++i) {
    expected_ampo += -t, "Cdag", i, "C", i + 1;
    expected_ampo += -t, "Cdag", i + 1, "C", i;
  }
  auto expected_H = toMPO(expected_ampo);

  // Create a random starting MPS
  auto state = InitState(sites);
  for (auto i : range1(N)) {
    if (i % 2 == 1)
      state.set(i, "1");
    else
      state.set(i, "0");
  }
  auto psi0 = randomMPS(state);

  // Run DMRG in 2 bases respectively
  auto sweeps = Sweeps(8);
  sweeps.maxdim() = 10, 20, 100, 200, 200;
  sweeps.cutoff() = 1E-10;
  auto [energy, psi] = dmrg(H, psi0, sweeps, {"Silent", true});
  auto [expected_energy, expected_psi] =
      dmrg(expected_H, psi0, sweeps, {"Silent", true});
  CHECK(energy == Approx(expected_energy).epsilon(1e-8));
}
