#ifndef __ITDVP_H_CMC__
#define __ITDVP_H_CMC__

#include <iomanip>
#include <typeinfo>

#include "FixedPointTensor.h"
#include "RandomUtility.h"
#include "Solver.h"
#include "itensor/all.h"
#include "uGauge.h"

using namespace itensor;
using namespace std;

/**
 * @brief Randomly initialize iMPS in the mixed canonical form.
 *
 * @param W The uniform MPO tensor.
 * @param is Physical index of MPO tensor `W`.
 * @param iwl Left index of MPO tensor `W`.
 * @param iwr Right index of MPO tensor `W`.
 * @param A The initial MPS tensor. If an empty tensor (as an ill-defined
 * ITensor) is provided, random elements will be filled in in accordance with
 * `is` and `D`.
 * @param D The maximum bond dimension.
 * @param ErrGoal The tolerance for orthonormal procedure.
 * @param MaxIter Maximum number of iteration for orthonormal procedure.
 * @param seed Random seed used when `A` is empty.
 * @returns tuple <ITensor, ITensor, ITensor, ITensor, ITensor, ITensor> - {AL,
 * AR, AC, C, La0, Ra0}.
 * @retval AL - Left-normalized iMPS tensor in the mixed canonical form.
 * @retval AR - Right-normalized iMPS tensor in the mixed canonical form.
 * @retval AC - The center site tensor, `AL` @f$\times@f$ `C`.
 * @retval C - The center matrix of iMPS in the mixed canonical form.
 * @retval La0 - Left fixed-point tensor of `A`.
 * @retval Ra0 - Right fixed-point tensor of `A`.
 *
 * @see https://scipost.org/SciPostPhysLectNotes.7
 */
tuple<ITensor, ITensor, ITensor, ITensor, ITensor, ITensor> itdvp_initial(
    const ITensor& W, const Index& is, const Index& iwl, const Index& iwr,
    ITensor& A, int D, Real ErrGoal, int MaxIter, RandGen::SeedType seed = 0) {
  Index il;
  if (!A) {
    il = Index(D, "Link");
    A = ITensor(is, il, prime(il, 2));
    auto rand = RandGen(seed);
    auto gen = [&rand]() { return rand.real(); };
    A.generate(gen);
    // A.randomize();
    // A.set(1,1,1,1.);
  } else {
    il = findIndex(A, "Link,0");
    if (!hasIndex(A, prime(il, 2))) {
      cout << "Error: " << __FUNCTION__ << ": A has wrong index structure"
           << endl;
      throw;
    }
  }
  auto ir = sim(il);
  global::IS._iwl = iwl;
  global::IS._iwr = iwr;
  global::IS._is = is;
  global::IS._il = il;
  global::IS._ir = ir;

  auto [AL, AR, C] = MixedCanonical(A, ErrGoal, MaxIter);
  auto AC = get_AC(AL, C);
  auto La0 = randomITensor(il, prime(il));
  auto Ra0 = randomITensor(ir, prime(ir));
  return {AL, AR, AC, C, La0, Ra0};
}

/**
 * @brief
 *
 * @param AL
 * @param C
 * @param AC
 * @return Real
 */
inline Real diff_ALC_AC(const ITensor& AL, const ITensor& C,
                        const ITensor& AC) {
  auto ALC = get_AC(AL, C);
  auto d = norm(ALC - AC);
  return d;
}

/**
 * @brief Variational uniform MPS (VUMPS) algorithm.
 *
 * @tparam TimeType
 * @param W The uniform MPO tensor.
 * @param AL Left-normalized iMPS tensor in the mixed canonical form.
 * @param AR Right-normalized iMPS tensor in the mixed canonical form.
 * @param AC The center site tensor, `AL` @f$\times@f$ `C`.
 * @param C The center matrix of iMPS in the mixed canonical form.
 * @param La0 Left fixed-point tensor of initial iMPS `A`.
 * @param Ra0 Right fixed-point tensor of initial iMPS `A`.
 * @param dt Length of time step. If `dt` is real, imaginary time evolution will
 * be performed; otherwise if `dt` is imaginary, real time evolution will be
 * performed.
 * @param args ErrGoal, MaxIter, used in applyExp and arnoldi.
 * @returns tuple <Real, Real, ITensor, ITensor> {en, err, LW, RW}.
 * @retval en - The variational ground state energy.
 * @retval err - Error measure.
 * @retval LW - Left fixed-point tensor of `W`.
 * @retval RW - Right fixed-point tensor of `W`.
 *
 * @note SciPost Phys. Lect. Notes 7 (2019), Algorithm 4 and Sec 5.3.
 * @see https://scipost.org/SciPostPhysLectNotes.7
 */
template <typename TimeType>
tuple<Real, Real, ITensor, ITensor> itdvp(const ITensor& W, ITensor& AL,
                                          ITensor& AR, ITensor& AC, ITensor& C,
                                          ITensor& La0, ITensor& Ra0,
                                          TimeType dt,
                                          Args& args = Args::global()) {
  // C --> L, R
  auto L = get_LR<LEFT>(C);
  auto R = get_LR<RIGHT>(C);

  // AL, AR, L, R --> LW, RW, enL, enR
  auto [LW, enL] = get_LRW<LEFT>(AL, W, R, La0, args);
  auto [RW, enR] = get_LRW<RIGHT>(AR, W, L, Ra0, args);
  auto en = 0.5 * (enL + enR);

  // LW, RW, W, AC, C --> AC, C
  if constexpr (is_same_v<TimeType, Real>) {
    if (isinf(dt))
      solve_gs(LW, RW, W, AC, C, args);
    else
      time_evolve(LW, RW, W, AC, C, dt, args);
  } else {
    time_evolve(LW, RW, W, AC, C, dt, args);
  }
  // AC, C --> AL, AR
  AL = get_AL(AC, C);
  AR = get_AR(AC, C);
  // AL = get_ALR_2 <LEFT>  (AC, C);
  // AR = get_ALR_2 <RIGHT> (AC, C);

  C /= norm(C);
  AC /= norm(AC);
  auto errL = diff_ALC_AC(AL, C, AC);
  auto errR = diff_ALC_AC(AR, C, AC);
  auto err = (errL > errR ? errL : errR);
  return {en, err, LW, RW};
}

#endif
