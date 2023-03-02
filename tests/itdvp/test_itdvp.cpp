#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "itdvp/GenMPO.h"
#include "itdvp/iTDVP.h"
#include "itensor/all.h"
#include "kbasis/OneParticleBasis.h"

using namespace itensor;
using namespace Catch;

TEST_CASE("Check ", "[FixPointTensor]") {
  int N = 8;
  auto t = 0.5;
  auto mu = 0.0;
  auto time_steps = 30;
  auto dt = 1e-12;
  auto D = 6;
  auto ErrGoalInit = 1e-8;
  auto MaxIter = 40;
  auto MaxIterInit = 20;
  auto SeedInit = 0;
  auto ErrGoal = 1e-12;
  // Why the arg ConserveQNs is not needed? Otherwise it will raise the error:
  // `Can't allocate quantum ITensor with undefined divergence`
  auto sites = Fermion(N, {"ConserveQNs", false});
  auto ampo = AutoMPO(sites);

  for (int i = 1; i <= N; ++i) {
    ampo += -mu, "N", i;
  }
  for (int i = 1; i < N; ++i) {
    ampo += -t, "Cdag", i, "C", i + 1;
    ampo += -t, "Cdag", i + 1, "C", i;
  }
  auto H = toMPO(ampo);

  auto [W, is, iwl, iwr] = get_W(H);

  auto A = ITensor();  // ill-defined tensor
  auto [AL, AR, AC, C, La0, Ra0] =
      itdvp_initial(W, is, iwl, iwr, A, D, ErrGoalInit, MaxIterInit, SeedInit);

  Args args = {"ErrGoal=", 1e-4, "MaxIter", MaxIter};
  ITensor LW, RW;
  Real en, err;
  for (int i = 1; i <= time_steps; i++) {
    std::cout << "time step " << i << std::endl;
    std::tie(en, err, LW, RW) = itdvp(W, AL, AR, AC, C, La0, Ra0, dt, args);
    std::cout << "energy, error = " << en << " " << err << std::endl;
    if (args.getReal("ErrGoal") > ErrGoal) args.add("ErrGoal=", err * 0.1);
  }

  println(inds(LW));
  println(inds(RW));
  println(inds(W));

  auto joystick = LW * W * RW;
  println(inds(joystick));
  println(inds(A));
  println(inds(AL));
  println(inds(AR));
  println(inds(AC));
  println(inds(C));
  // TODO: set prime index
  A.replaceInds();
  auto res = A * joystick;
  println(inds(res));
}
