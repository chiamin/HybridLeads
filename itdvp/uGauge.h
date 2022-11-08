#ifndef __uGauge_H_CMC__
#define __uGauge_H_CMC__
#include "itensor/all.h"
#include "IUtility.h"
#include "GlobalIndices.h"
using namespace itensor;
using namespace std;

inline ITensor applyTransfer (const ITensor& x, const ITensor& A1, const ITensor& A2, int pLV=2)
{
    auto re = x * A1;
    re *= prime(dag(A2),"Link");
    re.prime(-pLV,"Link");
    return re;
}

inline ITensor applyTransfer (const ITensor& x, const ITensor& A1, const ITensor& A2, const ITensor& W, int pLV=2)
{
    auto re = x * A1;
    re *= W;
    re *= prime(dag(A2));
    re.prime(-pLV,"Link");
    return re;
}

class TransferMatrix
{
    public:
        TransferMatrix (const ITensor& A1, const ITensor& A2, int pLV=2)
        : _A1 (A1)
        , _A2 (A2)
        , _pLV (pLV)
        {
            auto i1s = findInds (A1, "Link");
            auto i2s = findInds (A2, "Link");
            _dim = i1s(1).dim() * i2s(1).dim();
        }

        int size () const { return _dim; }
        void product (const ITensor& x, ITensor& xE) const
        {
            xE = applyTransfer (x, _A1, _A2, _pLV);
            assert (order(x) == order(xE));
            assert (hasIndex (xE, x.inds()(1)));
            assert (hasIndex (xE, x.inds()(2)));
        }

    private:
        ITensor _A1, _A2;
        int _dim, _pLV;
};

inline Real applyT_diff (const ITensor& A1, const ITensor& A2, const ITensor& X)
{
    auto reX = prime(X,2) * A1;
    reX *= prime (dag(A2), "Link");
    auto dd = reX - X;
    auto d = norm(dd);
    return d;
}

inline bool check_leading_eigen (const ITensor& A1, const ITensor& A2, const ITensor& X, Real crit=1e-12)
{
    auto d = applyT_diff (A1, A2, X);
    if (d > crit)
    {
        cout << __FUNCTION__ << ": err = " << d << endl;
        return false;
    }
    return true;
}
// SciPost Phys. Lect. Notes 7 (2019), Algorithm 1
//
// Solve
//
// --L--A-- = --AL--L--
//      |       |
//
// where
//
// ---     L------                                           ---A----
// |       |                                                    |
// l    =  |            is the left leading eigenvector of      |
// |       |                                                    |
// ---     dagL---                                           --dagA--
//
// and
//
// ----AL---    ---
// |   |      = |    = iut::Identity
// --dagAL--    ---
//
// Input:
//
//   --A--
//     |
//
// Output:
//
//   --AL--
//     |
//
//   --L--
//
tuple<ITensor,ITensor> Orthogonalize (const ITensor& A, Real errGoal=1e-15, int maxIter=10000)
{
    global::IS.check("A",A);
    auto is = global::IS.is();
    auto il = global::IS.il();
    auto il1 = prime(il);
    auto il2 = prime(il,2);

    // A = AL * L
    auto args = Args ({"MinDim",il.dim(),"MaxDim",il.dim()});
    auto [AL, L] = denmatDecomp (A, {is, il}, Fromleft, args);
    L /= norm(L);
    // Set indices
    auto ic = commonIndex (AL, L);
    AL.replaceInds ({ic}, {il2});
    L.replaceInds ({ic,il2}, {il1,il});

    // First error
    auto L_old = iut::Identity (L.inds());
    auto err = norm (L - L_old);

    // Iterate L*A = AL*L
    auto drop = ITensor();
    for(int i = 0; i < maxIter; i++)
    {
        // Arnoldi
        auto E = TransferMatrix (A, AL);
        arnoldi (E, L, {"ErrGoal=",0.1*err});
        L.takeReal();
        //
        tie (drop, L) = polar (L, {il1}, args);
        ic = commonIndex (L, drop);
        L.replaceInds ({ic}, {il1});
        L /= norm(L);
        // QR
        L_old = L;
        tie (AL, L) = polar (L*A, {is,il1}, args);
        L /= norm(L);
        // Replace indices
        ic = commonIndex (AL, L);
        AL.replaceInds ({il1,ic}, {il,il2});
        L.replaceInds ({ic,il2}, {il1,il});
        err = norm (L - L_old);
        cout << "err = " << err << endl;
        if (err < errGoal) break;
    }
    return {AL, L};
}

inline void Rotate_basis (ITensor& C, ITensor& AL, ITensor& AR)
{
    global::IS.check("C",C);
    global::IS.check("AL",AL);
    global::IS.check("AR",AR);
    auto is = global::IS.is();
    auto il = global::IS.il();
    auto il2 = prime(il,2);
    auto ir = global::IS.ir();
    auto ir2 = prime(ir,2);

    assert (abs(1.-norm(C)) < 1e-12);

    ITensor U(il2), S, V;
    svd (C, U, S, V);
    auto iU = commonIndex (U,S);
    auto iV = commonIndex (V,S);

    auto apply_rot = [&is, &C] (ITensor& U, ITensor& A, const Index& iU)
    {
        auto iU2 = prime(iU,2);
        U.setPrime (2, iU);
        A *= U;
        A *= noPrime(dag(U));
    };
    apply_rot (U, AL, iU);
    apply_rot (V, AR, iV);
    AL.replaceInds ({iU, prime(iU,2)}, {il, il2});
    AR.replaceInds ({iV, prime(iV,2)}, {ir, ir2});

    C = S;
    C = iut::hard_copy (C);
    C.replaceInds ({iU,iV}, {il2,ir2});

    global::IS.check("C",C);
    global::IS.check("AL",AL);
    global::IS.check("AR",AR);
}

inline tuple<ITensor,ITensor,ITensor> MixedCanonical (const ITensor& A, Real errGoal=1e-12, int maxIter=10000)
{
    global::IS.check("A",A);
    auto il = global::IS.il();
    auto ir = global::IS.ir();
    auto il2 = prime(il,2);
    auto ir2 = prime(ir,2);

    Args args = {"Cutoff",1e-14};

    auto [AL, L]  = Orthogonalize (A, errGoal, maxIter);
    auto ALtmp = swapInds (AL, {il}, {il2});
    auto [AR, C] = Orthogonalize (ALtmp, errGoal, maxIter);
    AR.replaceInds ({il,il2}, {ir,ir2});
    C.replaceInds ({il,prime(il)}, {il2,ir2});

    global::IS.check("AL",AL);
    global::IS.check("AR",AR);
    global::IS.check("C",C);

    Rotate_basis (C, AL, AR);

    assert (check_ortho (AL, il2));
    assert (check_ortho (AR, ir2));
    return {AL, AR, C};
}

//----------------------------------------------------

inline ITensor get_AC (const ITensor& ALR, ITensor C)
{
    global::IS.check("C",C);

    auto AC = ALR * C;
    AC.noPrime();

    global::IS.check("AC",AC);
    return AC;
}

inline ITensor get_polarU (const ITensor& T, const IndexSet& iUs)
{
    ITensor U (iUs), S, V;
    svd (T, U, S, V);

    V *= delta (S.inds());
    auto pU = U * V;
    return pU;
}

template <enumLR dir>
inline ITensor get_ALR_2 (const ITensor& AC, const ITensor& C)
{
    assert (abs(norm(AC)-1) < 1e-12);
    assert (abs(norm(C)-1) < 1e-12);
    global::IS.check("AC",AC);
    global::IS.check("C",C);

    auto is = global::IS.is();
    auto i1 = global::IS.il();
    auto i2 = global::IS.ir();
    if constexpr (dir == RIGHT)
        swap (i1,i2);

    auto UA = get_polarU (AC, {is,i1});
    auto UC = get_polarU (C, {prime(i1,2)});
    auto UCdag = dag(UC);
    UCdag.noPrime (prime(i2,2));
    auto ALR = UA * UCdag;

    assert (check_ortho (ALR, prime(i1,2)));
    if constexpr (dir == LEFT)
        global::IS.check("AL",ALR);
    else
        global::IS.check("AR",ALR);
    return ALR;
}

inline ITensor get_AL (const ITensor& AC, const ITensor& C)
{
    global::IS.check("AC",AC);
    global::IS.check("C",C);

    auto is = global::IS.is();
    auto il = global::IS.il();
    auto ir = global::IS.ir();

    auto Cdag = dag(C);
    Cdag.noPrime (prime(ir,2));
    auto AC_Cdag = AC * Cdag;

    auto AL = get_polarU (AC_Cdag, {is,il});

    global::IS.check("AL",AL);
    assert (check_ortho (AL, prime(il,2)));
    return AL;
}

inline ITensor get_AR (const ITensor& AC, const ITensor& C)
{
    global::IS.check("AC",AC);
    global::IS.check("C",C);

    auto is = global::IS.is();
    auto il = global::IS.il();
    auto ir = global::IS.ir();

    auto Cdag = dag(C);
    Cdag.noPrime (prime(il,2));
    auto Cdag_AC = AC * Cdag;

    auto AR = get_polarU (Cdag_AC, {prime(ir,2)});

    global::IS.check("AR",AR);
    assert (check_ortho (AR, prime(ir,2)));
    return AR;
}

template <enumLR dir>
inline ITensor get_LR (const ITensor& C)
{
    global::IS.check("C",C);

    auto ic = global::IS.il();
    auto iu = global::IS.ir();
    if constexpr (dir == RIGHT)
        swap (ic, iu);
    auto ic2 = prime(ic,2);
    auto iu2 = prime(iu,2);

    auto Cdag = dag(C);
    Cdag.prime (iu2);
    auto LorR = C * Cdag;
    LorR.prime(-2);
    //LorR.replaceInds ({iu,prime(iu)}, {ic,prime(ic)});
    auto trace = LorR * dag(delta(LorR.inds()));
    LorR /= iut::toReal (trace);

    if constexpr (dir == RIGHT)
        global::IS.check("R",LorR);
    else
        global::IS.check("L",LorR);
    return LorR;
}

template <enumLR dir>
ITensor get_LR_2 (const ITensor& AL, ITensor R, Real crit=1e-15)
{
    if constexpr (dir == RIGHT)
    {
        global::IS.check("AL",AL);
        global::IS.check("R",R);
    }
    else
    {
        global::IS.check("AR",AL);
        global::IS.check("L",R);
    }

    R.prime(2);
    auto E = TransferMatrix (AL, AL, -2);
    arnoldi (E, R, {"ErrGoal=",crit});
    R.prime(-2);

    R.takeReal();

    auto traceT = R * dag(delta(R.inds()));
    R /= iut::toReal (traceT);

//    assert (check_leading_eigen (AL, AL, R));

    if constexpr (dir == RIGHT)
        global::IS.check("R",R);
    else
        global::IS.check("L",R);
    return R;
}
#endif
