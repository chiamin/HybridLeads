#include <armadillo>
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
  auto t = GENERATE(0.5, 1.0);
  auto mu = GENERATE(0.0, 0.1);
  auto ham_elems = tight_binding_Hamilt(mat_dim, t, mu);
  double expected_elems[mat_dim][mat_dim] = {{-mu, -t, 0.0, 0.0},
                                             {-t, -mu, -t, 0.0},
                                             {0.0, -t, -mu, -t},
                                             {0.0, 0.0, -t, -mu}};
  arma::mat ham_mat(&ham_elems(0, 0), mat_dim, mat_dim);
  arma::mat expected_mat(&expected_elems[0][0], mat_dim, mat_dim);
  CHECK(approx_equal(ham_mat, expected_mat, "absdiff", 1e-12));
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

TEST_CASE("Check AutoMPO in real space basis", "[RealSpaceBasis]") {
  int N = 3;
  auto t = 0.5;
  auto mu = GENERATE(0.0, 0.1);
  auto sites = Fermion(N);
  auto ampo = AutoMPO(sites);

  // Construct AutoMPO
  for (int i = 1; i <= N; ++i) {
    ampo += -mu, "N", i;
  }
  for (int i = 1; i < N; ++i) {
    ampo += -t, "Cdag", i, "C", i + 1;
    ampo += -t, "Cdag", i + 1, "C", i;
  }
  auto H = toMPO(ampo, {"Exact=", true});

  // Construct full Hamiltonian matrix from MPO
  auto T = H(1) * H(2) * H(3);
  auto idxs = inds(T);

  // Note that the order of indices is not taken care by ITensor
  // See: https://github.com/chiamin/HybridLeads/issues/4
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
  auto [energy, psi] = dmrg(H, psi0, sweeps, {"Silent", true});

  // Run exact diagonalization
  Matrix U;
  Vector ens;
  diagHermitian(M, U, ens);
  double min_en = *min_element(ens.begin(), ens.end());

  // Compare ED ground state energy with DMRG
  CHECK(min_en == Approx(energy).epsilon(1e-8));
}

TEST_CASE("Check AutoMPO in hybrid basis", "[HybridBasis]") {
  int N = GENERATE(8, 16, 20);
  auto t = 0.5;
  auto mu = GENERATE(0.0, 0.1);
  auto sites = Fermion(N);
  auto ampo = AutoMPO(sites);

  // Construct AutoMPO in hybrid basis
  // I. 1st approach
  SECTION("Use basis transformation") {
    auto elems = tight_binding_Hamilt(N, t, mu);
    arma::mat ham_mat(&elems(0, 0), N, N);
    arma::mat Uik = arma::eye(N, N);
    arma::vec evals;
    arma::mat evecs;
    arma::mat sub_ham_mat = ham_mat.submat(0, 0, N / 2 - 1, N / 2 - 1);
    eig_sym(evals, evecs, sub_ham_mat);
    Uik.submat(0, 0, N / 2 - 1, N / 2 - 1) = evecs;
    arma::mat hybrid_ham_mat = Uik.t() * ham_mat * Uik;

    // Note that AutoMPO uses 1-index
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        double coef = hybrid_ham_mat(i, j);
        if (i == j) {
          ampo += coef, "N", i + 1;
        } else {
          ampo += coef, "Cdag", i + 1, "C", j + 1;
        }
      }
    }
  }

  // II. 2nd approach
  SECTION("Direct construction of AutoMPO with OneParticleBasis::C_op") {
    OneParticleBasis basis("sub_sys", N / 2, t, mu);
    // 1. The k-space part
    for (int k = 1; k <= N / 2; ++k) {
      auto coef = basis.en(k);
      ampo += coef, "N", k;
    }
    // 2. The mixing part
    auto coef = basis.C_op(N / 2, false);
    for (int k = 1; k <= N / 2; ++k) {
      ampo += -t * get<1>(coef[k - 1]), "Cdag", k, "C", N / 2 + 1;
      ampo += -t * get<1>(coef[k - 1]), "Cdag", N / 2 + 1, "C", k;
    }
    // 3. The real-space part
    for (int i = N / 2 + 1; i <= N; ++i) {
      ampo += -mu, "N", i;
    }
    for (int i = N / 2 + 1; i < N; ++i) {
      ampo += -t, "Cdag", i, "C", i + 1;
      ampo += -t, "Cdag", i + 1, "C", i;
    }
  }

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
