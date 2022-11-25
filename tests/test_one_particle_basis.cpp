#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../kbasis/OneParticleBasis.h"
#include "itensor/all.h"

using namespace itensor;
using namespace Catch;

TEST_CASE("Check matrix elements of tight-binding Hamiltonian",
          "[tight_binding_Hamilt]") {
  int mat_dim = 4;
  auto mat = tight_binding_Hamilt(mat_dim, 0.5, 1.);
  float expected_mat[4][4] = {{-1.0, -0.5, 0.0, 0.0},
                              {-0.5, -1.0, -0.5, 0.0},
                              {0.0, -0.5, -1.0, -0.5},
                              {0.0, 0.0, -0.5, -1.0}};

  for (int i = 0; i < mat_dim; i++) {
    for (int j = 0; j < mat_dim; j++) {
      CHECK(mat(i, j) == Approx(expected_mat[i][j]).epsilon(1e-12));
    }
  }
}

TEST_CASE("Check one particle basis", "[OneParticleBasis]") {
  int mat_dim = 4;
  OneParticleBasis basis("test_hamlt", mat_dim, 0.5, 1.);

  auto k_coef = basis.C_op(3, false);

  for (auto& tuple : k_coef) {
    std::cout << get<0>(tuple) << " " << get<1>(tuple) << " " << get<2>(tuple)
              << std::endl;
  }

  k_coef = basis.C_op(4, false);

  for (auto& tuple : k_coef) {
    std::cout << get<0>(tuple) << " " << get<1>(tuple) << " " << get<2>(tuple)
              << std::endl;
  }
}

TEST_CASE("Check auto MPO in real space basis", "[RealSpaceBasis]") {
  int N = 3;
  auto t = 0.5;
  auto mu = 0.1;
  auto sites = Fermion(N);
  auto ampo = AutoMPO(sites);

  for (int i = 1; i <= N; ++i) {
    ampo += -mu, "N", i;
  }

  for (int i = 1; i < N; ++i) {
    ampo += -t, "Cdag", i, "C", i + 1;
    ampo += -t, "Cdag", i + 1, "C", i;
  }

  // Construct full Hamiltonian matrix from MPO
  auto H = toMPO(ampo, {"Exact=", true});
  auto T = H(1) * H(2) * H(3);
  auto idxs = inds(T);

  auto [Comb, c] = combiner(idxs[0], idxs[2], idxs[4]);
  auto [Combp, cp] = combiner(idxs[1], idxs[3], idxs[5]);

  auto mat = Comb * T * Combp;
  auto M = Matrix(8, 8);

  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      M(i, j) = elt(mat, i, j);
    }
  }

  // Create a random starting MPS
  auto state = InitState(sites);
  for (auto i : range1(N)) {
    if (i % 2 == 1)
      state.set(i, "1");
    else
      state.set(i, "0");
  }
  auto psi0 = randomMPS(state);

  // Run DMRG
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 10, 20, 100, 200, 200;
  sweeps.cutoff() = 1E-8;
  auto [energy, psi] = dmrg(H, psi0, sweeps, {"Quiet", true});

  // Run exact diagonalization
  Matrix U;
  Vector ens;
  diagHermitian(M, U, ens);
  double min_en = *min_element(ens.begin(), ens.end());

  // Compare ED ground state energy with DMRG
  CHECK(min_en == Approx(energy).epsilon(1e-12));
}

TEST_CASE("Check auto MPO in hybrid basis", "[HybridBasis]") {
  int N = 4;
  auto t = 0.5;
  auto mu = 0.1;
  auto sites = Fermion(N);
  auto ampo = AutoMPO(sites);
  OneParticleBasis basis("test_hamlt", N, t, mu);

  for (int i = 1; i <= N; ++i) {
    ampo += -mu, "N", i;
  }

  for (int i = 1; i < N; ++i) {
    if (i < 2) {
      ampo += -t, "Cdag", i, "C", i + 1;
      ampo += -t, "Cdag", i + 1, "C", i;
    } else if (i == 2) {
      auto k_coef = basis.C_op(i, false);
      auto dag_k_coef = basis.C_op(i, true);
      for (int k; k < N; ++k) {
        ampo += -t * get<1>(k_coef[k]), "Cdag", i, "C", i + 1;
        ampo += -t * get<1>(dag_k_coef[k]), "Cdag", i + 1, "C", i;
      }
    } else {
      auto k_coef = basis.C_op(i, false);
      auto dag_k_coef = basis.C_op(i, true);
      for (int k; k < N; ++k) {
        for (int kp; kp < N; ++kp) {
          ampo += -t * get<1>(dag_k_coef[k]) * get<1>(k_coef[kp]), "Cdag", i,
              "C", i + 1;
          ampo += -t * get<1>(dag_k_coef[kp]) * get<1>(k_coef[k]), "Cdag",
              i + 1, "C", i;
        }
      }
    }
  }
  auto H = toMPO(ampo);

  auto control_ampo = AutoMPO(sites);
  for (int i = 1; i <= N; ++i) {
    control_ampo += -mu, "N", i;
  }

  for (int i = 1; i < N; ++i) {
    control_ampo += -t, "Cdag", i, "C", i + 1;
    control_ampo += -t, "Cdag", i + 1, "C", i;
  }
  auto control_H = toMPO(control_ampo);

  // Create a random starting MPS
  auto state = InitState(sites);
  for (auto i : range1(N)) {
    if (i % 2 == 1)
      state.set(i, "1");
    else
      state.set(i, "0");
  }
  auto psi0 = randomMPS(state);

  // Run DMRG
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 10, 20, 100, 200, 200;
  sweeps.cutoff() = 1E-8;
  // auto [energy, psi] = dmrg(H, psi0, sweeps, {"Quiet", true});
  auto [control_energy, control_psi] =
      dmrg(control_H, psi0, sweeps, {"Quiet", true});

  // printfln("Ground state energy = %.20f", energy);
  printfln("Ground state energy = %.20f", control_energy);
}
