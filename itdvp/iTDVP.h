#ifndef __ITDVP_H_CMC__
#define __ITDVP_H_CMC__
#include <typeinfo>
#include <iomanip>
#include "itensor/all.h"
#include "Solver.h"
#include "RandomUtility.h"
using namespace itensor;
using namespace std;

tuple <ITensor, ITensor, ITensor, ITensor, ITensor, ITensor>
itdvp_initial (const ITensor& W, const Index& is, const Index& iwl, const Index& iwr, ITensor& A, int D, Real ErrGoal, int MaxIter,
               RandGen::SeedType seed=0)
{
    Index il;
    if (!A)
    {
        il = Index(D,"Link");
        A = ITensor (is, il, prime(il,2));
        auto rand = RandGen (seed);
        auto gen = [&rand]() { return rand.real(); };
        A.generate (gen);
        //A.randomize();
        //A.set(1,1,1,1.);
    }
    else
    {
        il = findIndex (A, "Link,0");
        if (!hasIndex (A, prime(il,2)))
        {
            cout << "Error: " << __FUNCTION__ << ": A has wrong index structure" << endl;
            throw;
        }
    }
    auto ir = sim(il);

    global::IS._iwl = iwl;
    global::IS._iwr = iwr;
    global::IS._is = is;
    global::IS._il = il;
    global::IS._ir = ir;

    auto [AL, AR, C] = MixedCanonical (A, ErrGoal, MaxIter);

    auto AC = get_AC (AL, C);

    auto La0 = randomITensor (il, prime(il));
    auto Ra0 = randomITensor (ir, prime(ir));

    return {AL, AR, AC, C, La0, Ra0};
}

inline Real diff_ALC_AC (const ITensor& AL, const ITensor& C, const ITensor& AC)
{
    auto ALC = get_AC (AL, C);
    auto d = norm (ALC - AC);
    return d;
}

template <typename TimeType>
tuple <Real, Real, ITensor, ITensor>
itdvp
(const ITensor& W, ITensor& AL, ITensor& AR, ITensor& AC, ITensor& C, ITensor& La0, ITensor& Ra0, TimeType dt,
 Args& args=Args::global())
// args: ErrGoal, MaxIter, used in applyExp and arnoldi
{
    // C --> L, R
    auto L = get_LR <LEFT>  (C);
    auto R = get_LR <RIGHT> (C);

    // AL, AR, L, R --> LW, RW, enL, enR
    auto [LW, enL] = get_LRW <LEFT>  (AL, W, R, La0, args);
    auto [RW, enR] = get_LRW <RIGHT> (AR, W, L, Ra0, args);

    auto en = 0.5*(enL + enR);

    // LW, RW, W, AC, C --> AC, C
    if constexpr(is_same_v<TimeType,Real>)
    {
        if (isinf (dt))
            solve_gs (LW, RW, W, AC, C, args);
        else
            time_evolve (LW, RW, W, AC, C, dt, args);
    }
    else
    {
        time_evolve (LW, RW, W, AC, C, dt, args);
    }
    // AC, C --> AL, AR
    AL = get_AL (AC, C);
    AR = get_AR (AC, C);
    //AL = get_ALR_2 <LEFT>  (AC, C);
    //AR = get_ALR_2 <RIGHT> (AC, C);

    C /= norm(C);
    AC /= norm(AC);

    auto errL = diff_ALC_AC (AL, C, AC);
    auto errR = diff_ALC_AC (AR, C, AC);
    auto err = (errL > errR ? errL : errR);

    return {en, err, LW, RW};
}
#endif
