#ifndef __MPSUTILITY_H_CMC__
#define __MPSUTILITY_H_CMC__
#include <iomanip>
#include "itensor/mps/mpo.h"
#include "GeneralUtility.h"
using namespace std;
using namespace itensor;

namespace iut {

IndexSet siteInds (const MPO& mpo)
{
    vector<Index> ss;
    for(int i = 1; i <= length(mpo); i++)
    {
        auto s = siteIndex (mpo, i);
        ss.push_back (s);
    }
    return IndexSet (ss);
}

bool operator== (const IndexSet& s1, const IndexSet& s2)
{
    if (length(s1) != length(s2))
        return false;
    for(int i = 1; i <= length(s1); i++)
        if (s1(i) != s2(i))
            return false;
    return true;
}

inline bool operator!= (const IndexSet& s1, const IndexSet& s2)
{
    return !operator==(s1,s2);
}

// Apply <gate> on sites (<i>,<i>+1) on <psi>.
// <dir> = Fromleft:  the orthogonality center of the resulting MPS is on <i>+1, or
//         Fromright: the orthogonality center of the resulting MPS is on <i>
// The truncation settings in <args>: MaxDim, Cutoff
Spectrum apply_gate (MPS& psi, int i, const ITensor& gate, Direction dir, const Args& args=Args::global())
{
    int oc = orthoCenter(psi);
    assert (oc == i || oc == i+1);
    auto AA = psi(i) * psi(i+1) * gate;
    AA.noPrime();
    auto spec = psi.svdBond (i, AA, dir, args);
    return spec;
}

MPO Make_Nup_MPO (const Electron& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"Nup",i;
    }
    return toMPO (ampo);
}

MPO Make_Ndn_MPO (const Electron& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"Ndn",i;
    }
    return toMPO (ampo);
}

MPO Make_NMPO (const Electron& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"Ntot",i;
    }
    return toMPO (ampo);
}

MPO Make_NMPO (const Fermion& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"N",i;
    }
    return toMPO (ampo);
}

template <typename MPSType>
void check_indices (const MPSType& psi)
{
    int N = length(psi);
    for(int i = 1; i <= N; i++)
    {
        if (i != N)
            mycheck (commonIndex (psi(i), psi(i+1)), "MPS/MPO link-index check failed");
        auto is = findIndex (psi(i), "Site,0");
        mycheck (is, "MPS/MPO site-index check failed");
        if constexpr (is_same_v <MPSType,MPO>)
            mycheck (prime(is), "MPS/MPO site-index check failed");
    }
}

void print_MPO_tensors (const MPO& mpo, int i)
{
    cout << "site " << i << endl << endl;
    auto il = leftLinkIndex (mpo, i);
    auto ir = rightLinkIndex (mpo, i);
    for(int i1 = 1; i1 <= il.dim(); i1++)
        for(int i2 = 1; i2 <= ir.dim(); i2++)
        {
            auto Wij = mpo(i);
            auto is = findIndex (Wij, "Site,0");
            if (il) Wij *= setElt(dag(il)=i1);
            if (ir) Wij *= setElt(dag(ir)=i2);
            if (norm(Wij) != 0.)
            {
                if (!il)
                    cout << "iright = " << i2 << endl;
                else if (!ir)
                    cout << "ileft = " << i1 << endl;
                else
                    cout << "ileft, iright = " << i1 << " " << i2 << endl;
                Wij.permute ({is, prime(is)});
                //PrintData(Wij);

                // Print in matrix form
                cout << setprecision(4);
                cout << endl;
                for(int j1 = 1; j1 <= dim(is); j1++)
                {
                    for(int j2 = 1; j2 <= dim(is); j2++)
                        cout << setw(10) << elt (Wij, j2, j1);
                    cout << endl;
                }
                cout << endl;
            }
        }
}

void print_MPO (const MPO& mpo)
{
    for(int i = 1; i <= length(mpo); i++)
    {
        cout << "------------------------" << endl;
        print_MPO_tensors (mpo, i);
        cout << "------------------------" << endl;
    }
}

template <typename MPSTYPE>
inline Index leftIndex (const MPSTYPE& mps, int i)
{
    if (i == 1)
        return uniqueIndex (mps(1), mps(2), "Link");
    else
        return leftLinkIndex (mps, i);
}

template <typename MPSTYPE>
inline Index rightIndex (const MPSTYPE& mps, int i)
{
    int N = length (mps);
    if (i == N)
        return uniqueIndex (mps(N), mps(N-1), "Link");
    else
        return rightLinkIndex (mps, i);
}

Real inner_imps (const MPS& mps1, const MPS& mps2)
{
    int N = length(mps1);
    auto L = delta (dag(leftIndex(mps1,1)), prime(leftIndex(mps2,1)));
    for(int i = 1; i < N; i++)
    {
        L *= mps1(i);
        L *= prime(dag(mps2(i)),"Link");
    }
    L *= mps1(N);
    auto il = leftIndex (mps1, N);
    L *= prime(dag(mps2(N)),{il});
    return elt(L);
}

void renewLinkInds (MPS& mps, int ibeg, int iend)
{
    auto ir0 = rightLinkIndex (mps, ibeg);
    auto ir = sim(ir0);
    mps.ref(ibeg).replaceInds ({ir0}, {ir});
    for(int i = ibeg+1; i < iend; i++)
    {
        // Left index
        auto il0 = leftLinkIndex (mps, i);
        mps.ref(i).replaceInds ({il0}, {dag(ir)});
        // Right index
        ir0 = rightLinkIndex (mps, i);
        ir = sim(ir0);
        mps.ref(i).replaceInds ({ir0}, {ir});
    }
    auto il0 = leftLinkIndex (mps, iend);
    mps.ref(iend).replaceInds ({il0}, {dag(ir)});
}

void MPSReplaceTensors (MPS& psi, const MPS& subpsi, int ibeg)
{
    for(int i = 1; i <= length(subpsi); i++)
    {
        int i0 = i+ibeg-1;
        auto iis = findIndex (psi(i0), "Site");
        auto iis2 = findIndex (subpsi(i), "Site");
        Index iil;
        if (i == 1)
            iil = leftLinkIndex (psi, i0);
        else if (i == length(subpsi))
            iil = rightLinkIndex (psi, i0);
        auto A = subpsi(i);
        A.replaceInds ({iis2}, {iis});
        psi.ref(i0) = A;
        if (iil)
        {
            mycheck (dim(iil) == 1, "not dummy index");
            psi.ref(i0) *= setElt(iil=1);
        }
    }
}
} // namespace
#endif
