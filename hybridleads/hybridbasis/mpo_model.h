#ifndef __MPOMODEL_H__
#define __MPOMODEL_H__

#include <armadillo>

#include "itensor/all.h"
#include "kbasis/OneParticleBasis.h"

using namespace itensor;

class MPOModel {
 public:
  /**
   * @brief Base class of MPO model in the hybrid basis jointed by left lead
   * (real), system (k), and right lead (real).
   *
   * @param left_size
   * @param system_size
   * @param right_size
   * @throws `invalid_argument` when lead size is less than 4.
   * @note User should implement the virtual member functions in child class.
   */
  MPOModel(int left_size, int system_size, int right_size) {
    if ((left_size < 4) || (right_size < 4)) {
      throw std::invalid_argument(
          "Lead size should be equal or greater than 4.");
    }
    n_left = left_size;
    n_sys = system_size;
    n_right = right_size;
    n_tot = left_size + system_size + right_size;
  }

  /**
   * @brief Return the size of each block in joint system.
   * @returns std::tuple<int, int, int> (left_size, system_size, right_size).
   */
  std::tuple<int, int, int> sizes() {
    return std::make_tuple(n_left, n_sys, n_right);
  }

  /**
   * @brief Return `SiteSet` of this model.
   * @returns SiteSet
   */
  SiteSet sites() { return _sites; }

  /**
   * @brief Return the `AutoMPO` instance in hybrid basis.
   * @returns MPO
   */
  MPO mpo() { return toMPO(ampo); }

  /**
   * @brief Return the single particle Hamiltonian in "standard" basis.
   * @returns arma::mat - sp_ham.
   */
  arma::mat single_particle_ham() { return sp_ham; }

  /**
   * @brief Return the single particle Hamiltonian in hybrid basis.
   * @returns arma::mat - hybrid_ham.
   */
  arma::mat hybrid_basis_ham() { return hybrid_ham; }

 protected:
  int n_left, n_sys, n_right, n_tot;
  SiteSet _sites;
  AutoMPO ampo;
  arma::mat sp_ham;
  arma::mat hybrid_ham;

  /**
   * @brief Rotate the basis of single particle Hamiltonian into the hybrid
   * basis.
   */
  void basis_transformer() {
    arma::vec evals;
    arma::mat evecs;
    arma::mat sub_ham_mat =
        sp_ham.submat(n_left, n_left, n_left + n_sys - 1, n_left + n_sys - 1);
    arma::eig_sym(evals, evecs, sub_ham_mat);
    arma::mat unitary_mat = arma::eye(n_tot, n_tot);
    unitary_mat.submat(n_left, n_left, n_left + n_sys - 1, n_left + n_sys - 1) =
        evecs;
    hybrid_ham = unitary_mat.t() * sp_ham * unitary_mat;
  }

  virtual void gen_single_particle_ham() {}

  virtual void gen_auto_mpo() {}
};

class TightBinding : public MPOModel {
 public:
  /**
   * @brief Construct the tight-binding model in hybridbasis.
   *
   * @param left_size
   * @param system_size
   * @param right_size
   * @param args Arguments containing the model parameters, allowed keywords:
   * (1) t_left, t_left_sys, t_sys, t_right_sys, t_right
   * (2) mu_left, mu_sys, mu_right
   */
  TightBinding(int left_size, int system_size, int right_size, Args const& args)
      : MPOModel(left_size, system_size, right_size) {
    t_left = args.getReal("t_left", 0.0);
    t_left_sys = args.getReal("t_left_sys", 0.0);
    t_sys = args.getReal("t_sys", 0.0);
    t_right_sys = args.getReal("t_right_sys", 0.0);
    t_right = args.getReal("t_right", 0.0);
    mu_left = args.getReal("mu_left", 0.0);
    mu_sys = args.getReal("mu_sys", 0.0);
    mu_right = args.getReal("mu_right", 0.0);
    _sites = Fermion(n_tot);
    ampo = AutoMPO(_sites);
    gen_auto_mpo();
  }

 protected:
  Real t_left, t_left_sys, t_sys, t_right_sys, t_right;
  Real mu_left, mu_sys, mu_right;

  arma::mat block_tight_binding_ham(int n, Real t, Real mu) {
    auto elems = tight_binding_Hamilt(n, t, mu);
    arma::mat block_ham(&elems(0, 0), n, n);
    return block_ham;
  }

  void gen_single_particle_ham() {
    sp_ham = arma::zeros(n_tot, n_tot);
    arma::mat block_11 = block_tight_binding_ham(n_left, t_left, mu_left);
    arma::mat block_22 = block_tight_binding_ham(n_sys, t_sys, mu_sys);
    arma::mat block_33 = block_tight_binding_ham(n_right, t_right, mu_right);
    sp_ham.submat(0, 0, n_left - 1, n_left - 1) = block_11;
    sp_ham.submat(n_left, n_left, n_left + n_sys - 1, n_left + n_sys - 1) =
        block_22;
    sp_ham.submat(n_left + n_sys, n_left + n_sys, n_tot - 1, n_tot - 1) =
        block_33;
    sp_ham(n_left - 1, n_left) = -t_left_sys;
    sp_ham(n_left, n_left - 1) = -t_left_sys;
    sp_ham(n_left + n_sys - 1, n_left + n_sys) = -t_right_sys;
    sp_ham(n_left + n_sys, n_left + n_sys - 1) = -t_right_sys;
  }

  void gen_auto_mpo() {
    gen_single_particle_ham();
    basis_transformer();
    for (int i = 0; i < n_tot; ++i) {
      for (int j = 0; j < n_tot; ++j) {
        double coef = hybrid_ham(i, j);
        if (i == j) {
          ampo += coef, "N", i + 1;
        } else {
          ampo += coef, "Cdag", i + 1, "C", j + 1;
        }
      }
    }
  }
};

class AndersonImpurity : public MPOModel {
 public:
  AndersonImpurity(int left_size, int system_size, int right_size,
                   Args const& args)
      : MPOModel(left_size, system_size, right_size) {}
};

class KondoImpurity : public MPOModel {
 public:
  KondoImpurity(int left_size, int system_size, int right_size,
                Args const& args)
      : MPOModel(left_size, system_size, right_size) {}
};

#endif
