#ifndef KBASIS_ONEPARTICLEBASIS_H_
#define KBASIS_ONEPARTICLEBASIS_H_

#include "GeneralUtility.h"
#include "itensor/all.h"

using Real = itensor::Real;
using tuple = std::tuple;

/**
 * @brief Generate the tight-binding Hamiltonian in single particle basis.
 *
 * @param L The size of matrix, L x L.
 * @param t The hopping strength.
 * @param mu The on-site energy.
 * @param damp_fac The base of power law hopping. Default to 1 as taking no
 * effect.
 * @param damp_from_right
 * @param verbose Print out the matrix elements.
 * @returns Matrix - The Hamiltonian matrix, which is tridiagonal.
 */
itensor::Matrix tight_binding_Hamilt(
    int L, Real t, Real mu, Real damp_fac = 1., bool damp_from_right = true,
    bool verbose = false
) {
  itensor::Matrix H(L, L);
  for (int i = 0; i < L; i++) {
    H(i, i) = -mu;
    if (i != L - 1) {
      int damp_dist = (damp_from_right ? L - 2 - i : i);
      Real ti = t * std::pow(damp_fac, damp_dist);
      H(i, i + 1) = -ti;
      H(i + 1, i) = -ti;
      if (verbose) std::cout << "Hk, t " << i << " = " << ti << std::endl;
    }
  }
  return H;
}

class OneParticleBasis {
 public:
  OneParticleBasis() {}
  OneParticleBasis(const string& name, const itensor::Matrix& H) : _name(name), _H(H) {
    itensor::diagHermitian(H, _Uik, _ens);
  }
  OneParticleBasis(
      const string& name, int L, Real t, Real mu, Real damp_fac = 1.,
      bool damp_from_right = true, bool verbose = false
  )
      : _name(name) {
    _H = tight_binding_Hamilt(L, t, mu, damp_fac, damp_from_right, verbose);
    itensor::diagHermitian(_H, _Uik, _ens);
  }
  OneParticleBasis(const string& name, int L) : _name(name) {
    _H = itensor::Matrix(L, L);
    itensor::diagHermitian(_H, _Uik, _ens);
  }

  // Functions that every basis class must have
  const string& name() const { return _name; }
  std::vector<tuple<int, auto, bool>> C_op(int i, bool dag) const;
  Real en(int k) const {
    mycheck(k > 0 and k <= _ens.size(), "out of range");
    return _ens(k - 1);
  }
  Real mu(int k) const {
    mycheck(k > 0 and k <= _ens.size(), "out of range");
    return -_H(k - 1, k - 1);
  }
  int size() const { return _ens.size(); }

  void write(ostream& s) const {
    itensor::write(s, _name);
    itensor::write(s, _Uik);
    itensor::write(s, _ens);
  }
  void read(istream& s) {
    itensor::read(s, _name);
    itensor::read(s, _Uik);
    itensor::read(s, _ens);
  }

 private:
  string _name;
  itensor::Matrix _Uik, _H;
  itensor::Vector _ens;
};

/**
 * @brief Get the operator information in this basis for the operator Cdag_i,
 * @f{eqnarray*}
 *    C_i = \sum_k U_ik C_k \\
 *    C_i^\dagger = \sum_k U_ik^* C_k^\dagger
 * @f}
 *
 * @param i The real-space site index, 1-index.
 * @param dag Whether the operator has a dagger or not.
 * @returns vector<tuple<int,auto,bool>> - The basis index (k, 0-index) in
 * ascending order, coefficient U_ik, and whether the operator has a dagger or
 * not.
 */
std::vector<tuple<int, auto, bool>> OneParticleBasis ::C_op(int i, bool dag) const {
  mycheck(i > 0 and i <= nrows(_Uik), "out of range");

  auto tmp = _Uik(0, 0);
  std::vector<tuple<int, decltype(tmp), bool>> k_coef_dag;

  for (int k = 0; k < this->size(); k++)  // Here k is zero-index
  {
    if (dag) {
      k_coef_dag.emplace_back(k + 1, iut::conj(_Uik(i - 1, k)), dag);
    } else {
      k_coef_dag.emplace_back(k + 1, _Uik(i - 1, k), dag);
    }
  }
  return k_coef_dag;
}

auto write(ostream& s, const OneParticleBasis& t) { t.write(s); }
auto read(istream& s, OneParticleBasis& t) { t.read(s); }

#endif
