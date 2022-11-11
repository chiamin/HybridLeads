#ifndef __ENTANGLEMENT_H_CMC__
#define __ENTANGLEMENT_H_CMC__
#include "ContainerUtility.h"
#include "GeneralUtility.h"
using namespace itensor;

Real EntangEntropy (const Spectrum& spec)
{
    Real S = 0.;
    for(int i = 1; i <= spec.size(); i++)
    {
        auto p = spec.eig(i);
        S += -p * log(p);
    }
    return S;
}


Real EntangEntropy_singular (const ITensor& S)
{
    assert (order(S) == 2);
    assert (dim(S.inds()(1)) == dim(S.inds()(2)));

    int m = S.inds()(1).dim();
    Real entropy = 0.;
    for(int i = 1; i <= m; i++)
    {
        Real s = elt (S,i,i);
        Real p = s*s;
        entropy += -p * log(p);
    }
    return entropy;
}

Real EntangEntropy_rho (const ITensor& rho, Real cutoff=1e-12)
{
    auto [U,D] = diagHermitian (rho, {"Cutoff",cutoff});

    int m = D.inds()(1).dim();
    Real S = 0.;
    for(int i = 1; i <= m; i++)
    {
        auto p = elt (D,i,i);
        if (p > cutoff)
        {
            S += -p * log(p);
        }
    }
    return S;
}

ITensor get_rho_brute_force (const MPS& psi, vector<int>& sites)
{
    std::sort (sites.begin(), sites.end());

    // Left
    ITensor L(1.);
    int imid = sites.at(sites.size()/2);
    for(int i = 1; i <= length(psi); i++)
    {
        ITensor A = psi(i);
        ITensor Adag;
        if (iut::in_vector (sites, i))
            Adag = prime(dag(A));
        else
            Adag = prime(dag(A),"Link");

        L *= A;
        L *= Adag;
        if (i == imid)
            break;
    }
    // Right
    ITensor R(1.);
    for(int i = length(psi); i >= 1; i--)
    {
        if (i == imid)
            break;
        ITensor A = psi(i);
        ITensor Adag;
        if (iut::in_vector (sites, i))
            Adag = prime(dag(A));
        else
            Adag = prime(dag(A),"Link");

        R *= A;
        R *= Adag;
    }
    ITensor rho = L * R;
    mycheck (order(rho) == 2*sites.size(), "size not match");
    return rho;
}

Real get_entang_entropy_brute_force (const MPS& psi, vector<int>& sites)
{
    auto rho = get_rho_brute_force (psi, sites);
    auto EE = EntangEntropy_rho (rho);
    return EE;
}

// Measure the entanglement entropy of a system between site i1 and i2
// <args> is used 
Real get_entang_entropy (MPS psi, int i1, int i2, const Args& args=Args::global())
{
    if (i1 > i2)
        swap (i1,i2);
    // Pull out a sub-chain of MPS between site i1 and i2
    //
    //        |  |  |  |  |       |  |
    //      --O--C--O--O--O--...--O--O--
    //
    //
    //        |  |  |  |  |       |  |
    //  =>  --O--(AA)--O--O--...--O--O--
    //
    //
    //           |
    //           U
    //            \
    //            S 
    //        |    \|  |  |       |  |
    //  =>  --O-----V--O--O--...--O--O--
    //
    //           |
    //           U-
    //        |    \|  |  |       |  |
    //  =>  --O-----(AA)--O--...--O--O--
    //
    //
    //           |  |
    //           O--U
    //               \
    //               S 
    //        |       \|  |       |  |
    //  =>  --O--------V--O--...--O--O--
    //
    //
    //  =>  ...
    //
    //
    //           |  |  |      |
    //           O--O--O--..--U
    //                         \
    //                          S 
    //        |                  \|  |
    //  =>  --O-------------------V--O--
    //
    // S^2 will be the eigenvalues of the reduced density matrix for the sub-chain
    // The entanglement entropy can be computed from the singular values S
    //
    psi.position(i1);
    ITensor C = psi(i1);
    Index ii;
    for(int i = i1; i < i2; i++)
    {
        auto AA = C * psi(i+1);
        auto inds = vector<Index> { siteIndex(psi,i) };
        if (i != i1)
            inds.push_back (ii);
        auto [U,S,V] = svd (AA, inds, args);
        C = S * V;
        ii = commonIndex (C, U);
    }
    auto AA = C * psi(i2+1);
    auto [U,S,V] = svd (AA, {siteIndex(psi,i2), ii}, args);
    return EntangEntropy_singular (S);
}
#endif
