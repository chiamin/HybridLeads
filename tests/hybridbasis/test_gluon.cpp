#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/trompeloeil.hpp>
#include <list>

#include "hybridbasis/fixed_point_tensor.h"
#include "hybridbasis/gluon.h"
#include "hybridbasis/mpo_model.h"
#include "hybridbasis/utils.h"

using namespace itensor;
using namespace Catch;

TEST_CASE("Check gluon indices", "[TestCommonInds]") {
  int n_left = 4;
  int n_sys = GENERATE(2, 4, 6);
  int n_right = n_left;
  int n_tot = n_left + n_sys + n_right;
  Args model_args = {"t_left",      0.5,  "t_left_sys", 0.01, "t_sys",       0.5,
                     "t_right_sys", 0.01, "t_right",    0.5,  "mu_left",     0.0,
                     "mu_sys",      0.0,  "mu_right",   0.0,  "ConserveQNs", false};
  TightBinding model(n_left, n_sys, n_right, model_args);
  auto mpo = model.mpo();
  Args itdvp_args = {"max_bond_dim", 6};
  auto sites = Fermion(n_sys + 2, {"ConserveQNs", false});

  Gluon gluon(mpo, sites, n_left - 1, n_right - 1, itdvp_args);
  auto sys_mpo = gluon.sys_mpo();
  auto left_env = gluon.left_env();
  auto right_env = gluon.right_env();
  auto init_state = gluon.random_init_state();

  CHECK(order(sys_mpo(1)) == 4);
  CHECK(order(sys_mpo(length(sites))) == 4);
  CHECK(order(init_state(1)) == 3);
  CHECK(order(init_state(length(sites))) == 3);

  std::list<IndexSet> indices = {
      commonInds(sys_mpo(1), left_env), commonInds(sys_mpo(length(sites)), right_env),
      commonInds(init_state(1), left_env),
      commonInds(init_state(length(sites)), right_env)};

  for (IndexSet idx : indices) {
    CHECK(length(idx) == 1);
    CHECK(hasTags(idx(1), "Link"));
  }
}
