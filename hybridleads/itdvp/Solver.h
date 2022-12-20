#ifndef __SOLVER_H_CMC__
#define __SOLVER_H_CMC__

#include "GlobalIndices.h"
#include "itensor/all.h"

using namespace itensor;
using namespace std;

void apply_Heff(const ITensor& LW, const ITensor& RW, const ITensor& W,
                ITensor& AC, ITensor& C) {
  global::IS.check("LW", LW);
  global::IS.check("RW", RW);
  global::IS.check("W", W);
  global::IS.check("AC", AC);
  global::IS.check("C", C);

  auto iwl = global::IS.iwl();
  auto iwr = global::IS.iwr();

  // Update AC
  AC *= LW;
  AC *= W;
  AC *= RW;
  AC.noPrime();

  // Update C
  auto RW0 = replaceInds(RW, {iwr}, {iwl});
  C.noPrime();
  C *= LW;
  C *= RW0;
  C.prime();

  C /= norm(C);
  AC /= norm(AC);

  global::IS.check("AC", AC);
  global::IS.check("C", C);
}

void solve_gs(const ITensor& LW, const ITensor& RW, const ITensor& W,
              ITensor& AC, ITensor& C, const Args& args = Args::global())
// args: ErrGoal, MaxIter
{
  global::IS.check("LW", LW);
  global::IS.check("RW", RW);
  global::IS.check("W", W);
  global::IS.check("AC", AC);
  global::IS.check("C", C);

  auto iwl = global::IS.iwl();
  auto iwr = global::IS.iwr();

  // Update AC
  LocalOp Heff1(W, LW, RW, {"numCenter=", 1});
  davidson(Heff1, AC, args);

  // Update C
  auto RW0 = replaceInds(RW, {iwr}, {iwl});
  LocalOp Heff0(LW, RW0, {"NumCenter=", 0});
  auto C0 = noPrime(C);
  davidson(Heff0, C0, args);
  C = prime(C0, 2);

  C /= norm(C);
  AC /= norm(AC);

  global::IS.check("AC", AC);
  global::IS.check("C", C);
}

template <typename TimeType>
void time_evolve(const ITensor& LW, const ITensor& RW, const ITensor& W,
                 ITensor& AC, ITensor& C, TimeType dt,
                 const Args& args = Args::global())
// args: ErrGoal, MaxIter
{
  global::IS.check("LW", LW);
  global::IS.check("RW", RW);
  global::IS.check("W", W);
  global::IS.check("AC", AC);
  global::IS.check("C", C);

  auto iwl = global::IS.iwl();
  auto iwr = global::IS.iwr();

  // Update AC
  LocalOp Heff1(W, LW, RW, {"numCenter=", 1});
  applyExp(Heff1, AC, -dt, args);

  // Update C
  auto RW0 = replaceInds(RW, {iwr}, {iwl});
  LocalOp Heff0(LW, RW0, {"NumCenter=", 0});
  auto C0 = noPrime(C);
  applyExp(Heff0, C0, -dt, args);
  C = prime(C0, 2);

  global::IS.check("AC", AC);
  global::IS.check("C", C);
}

#endif
