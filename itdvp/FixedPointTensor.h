#ifndef __FIXEDPOINTTENSOR_H_CMC__
#define __FIXEDPOINTTENSOR_H_CMC__
#include <iomanip>
#include "itensor/all.h"
#include "IUtility.h"
#include "uGauge.h"

class SpecialTransferMatrix
// For computing (x|[1 - E + |R)(1|]
{
    public:
        SpecialTransferMatrix (const ITensor& AL, const ITensor& R)
        : _AL (AL)
        , _R (R)
        {
            assert (order(R) == 2);
        }

        int size () const { return _R.inds()(1).dim() * _R.inds()(2).dim(); }
        void product (const ITensor& x, ITensor& re) const
        {
            assert (order(x) == 2);
            // (x|1
            re = x;

            // (x|E
            ITensor xE = applyTransfer (x, _AL, _AL);
            re -= xE;

            // (x|R)(1|
            auto xR = x * _R;
            ITensor xRI;
            if (isReal (xR))
                xRI = iut::Identity (re.inds(), elt(xR));
            else
                xRI = iut::Identity (re.inds(), eltC(xR));
            re += xRI;

            assert (order(x) == order(re));
            assert (hasIndex (re, x.inds()(1)));
            assert (hasIndex (re, x.inds()(2)));
        }

    private:
        ITensor _AL, _R;
};

template <enumLR dir>
tuple<ITensor,Real> get_LRW (const ITensor& AL, const ITensor& W, const ITensor& R, ITensor& La0, const Args& args)
// PHYSICAL REVIEW B 97, 045145 (2018), Algorithm 6
// W should have indices: 1. up site index, 2. down site index, 3. left link index, 4. right link index
// AL should have indices: 1. site index, 2. left link index, 3. right link index
//
//  --------AL----
//  |       |
//  L--(b)--W--(a)
//  |       |
//  --------AL----
//
//  (2)--AL--(3)--C--
//       |
//      (1)
//
//      (1)
//       |
//  (3)--W--(4)
//       |
//      (2)
//
{
    // Check W is of Schur form:
    // Wba != 0 only for b >= a (lower triangled)
    global::IS.check("W",W);

    Index ia, iWa, iWb;
    if constexpr (dir == LEFT)
    {
        global::IS.check("AL",AL);
        global::IS.check("R",R);
        ia = global::IS.il();
        iWa = global::IS.iwr();
        iWb = global::IS.iwl();
    }
    else
    {
        global::IS.check("AR",AL);
        global::IS.check("L",R);
        ia = global::IS.ir();
        iWa = global::IS.iwl();
        iWb = global::IS.iwr();
    }

    // For LW: solve from a=dW to 1
    // For RW: solve from a=1 to dW
    vector<int> solve_order;
    solve_order.push_back(2);
    for(int i = iWa.dim(); i >= 1; i--)
    {
        if (i != 2)
            solve_order.push_back (i);
    }
    if constexpr (dir == RIGHT)
        reverse (solve_order.begin(), solve_order.end());

    auto M = SpecialTransferMatrix (AL, R);

    auto LRW = ITensor (ia, iWb, prime(ia));
    // 1) b == a == (dW or 1), Waa = I  -->  La = I
    for(int i = 1; i <= ia.dim(); i++)
        LRW.set (i, solve_order.front(), i, 1.);

    // 2)
    ITensor Ya;
    for(int i = 1; i < solve_order.size(); i++)
    {
        int a = solve_order.at(i);
        auto Wa = W * setElt (dag(iWa)=a);
        auto Waa = Wa * setElt (dag(iWb)=a);
        assert (order(Waa) == 2);

        Ya = applyTransfer (LRW, AL, AL, Wa);

        // Waa == 0
        if (i != solve_order.size()-1)
        {
            assert (norm(Waa) == 0.);
            auto Y = Ya * setElt (dag(iWb)=a);

            assert (order(Y) == 3);
            LRW += Y;
        }
        else
        //solve (L1|[1 - E + |R)(1|] = (Ya| - (Ya|R)(1|
        {
            auto YR = Ya * R;                               // (Ya|R)
            ITensor YRI;                                    // (Ya|R)(1|
            if (isReal (YR))
                YRI = iut::Identity (Ya.inds(), elt(YR));
            else
                YRI = iut::Identity (Ya.inds(), eltC(YR));
            auto Ya2 = Ya - YRI;                            // (Ya| - (Ya|R)(1|
            gmres (M, Ya2, La0, args);

            LRW += La0 * setElt (dag(iWb)=a);
        }

#ifdef DEBUG
        // Check fixed point
        auto crit = args.getReal("debug_err_crit",1e-8);
        auto LTW = applyTransfer (LRW, AL, AL, Wa);
        auto LTWa = LRW * setElt(dag(iWb)=a);
        auto d = LTW - LTWa;
        if (a != solve_order.back())
        {
            if (norm(d) > crit)
            {
                cout << "d = " << norm(d) << endl;
                throw;
            }
        }
        else
        {
            auto e = Ya * R;
            ITensor Ie;
            if (isReal(e))
                Ie = iut::Identity (d.inds(), elt(e));
            else
                Ie = iut::Identity (d.inds(), eltC(e));
            auto dI = d - Ie;
            auto d2 = norm(dI);
            cout << "d = " << d2 << endl;
            /*if (d2 > 1e-4)
            {
                cout << "d = " << d2 << endl;
                //throw;
            }*/
        }
#endif
    }
    auto en = iut::toReal (Ya * R); // energy density
    if constexpr (dir == LEFT)
        global::IS.check("LW",LRW);
    else
        global::IS.check("RW",LRW);

    return {LRW, en};
}

#endif
