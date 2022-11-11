#ifndef __FERMIMEASURE_H_CMC__
#define __FERMIMEASURE_H_CMC__
#include "itensor/all.h"
using namespace std;
using namespace itensor;

ITensor mult_onsite (const ITensor& op1, ITensor op2)
{
    op2 *= prime(op1);
    op2.mapPrime(2,1);
    return op2;
}

ITensor LeftTensor (const MPS& psi, const MPS& phi, int i, const ITensor& O=ITensor())
{
    ITensor re = psi.A(i);
    if (O)
    {
        re *= O;
        re.noPrime();
    }
    re *= dag (prime (psi.A(i), rightLinkIndex(psi,i)));
    return re;
}

// Multiply psi.A(i) * Oi * phi.A(i) to re.
// If <close> is "left" or "right", close the link on the left or right.
void psiA_Oi_phiA (ITensor& re, const MPS& psi, const MPS& phi, int i, const ITensor& O=ITensor(), string close="")
{
    re *= psi.A(i);
    if (O)
        re *= O;
    else
        re.prime ("Site");

    auto tmp = prime(phi.A(i));
    if (close == "right" && i != length(phi))
    {
        tmp.noPrime (prime(rightLinkIndex(phi,i)));
    }
    else if (close == "left" && i != 1)
    {
        tmp.noPrime (prime(rightLinkIndex(phi,i-1)));
    }
    re *= dag(tmp);
}

// i < j for j in js
vector<Real> Measure_Delta_ij (const SiteSet& st, int i, const vector<int>& js, const MPS& psi)
{
    vector<Real> re (js.size());

    ITensor Ciup = mult_onsite (st.op("F",i), st.op("Cdagup",i)),
             Cidn = mult_onsite (st.op("F",i), st.op("Cdagdn",i));

    ITensor L1 = LeftTensor (psi, psi, i, Ciup),
             L2 = LeftTensor (psi, psi, i, Cidn);

    int      ji = 0;
    int      j  = js.at(ji);

    Real sqrt2 = sqrt(2.);
    for(int k = i+1; k <= js.back(); k++)
    {
        // Check orthogonality center
        if (orthoCenter(psi) < i || orthoCenter(psi) > j)
        {
            cout << "Error: Measure_Fermi_Ops: orthoCenter(psi) has between i and j" << endl;
            cout << "       orth,i,j = " << orthoCenter(psi) << " " << i << " " << j << endl;
            throw;
        }

        // Measure
        if (k == j)
        {
            ITensor Cjdn = st.op("Cdagdn",j),
                     Cjup = st.op("Cdagup",j);
            ITensor tmp1 = L1,
                     tmp2 = L2;
            psiA_Oi_phiA (tmp1, psi, psi, j, Cjdn, "right");
            psiA_Oi_phiA (tmp2, psi, psi, j, Cjup, "right");

            re.at (ji) = (tmp1.real()-tmp2.real()) / sqrt2;

            if (k != js.back())
                j = js.at(++ji);
            else
                return re;
        }

        // Apply Electronic string operator
        psiA_Oi_phiA (L1, psi, psi, k, st.op("F",k));
        psiA_Oi_phiA (L2, psi, psi, k, st.op("F",k));
    }
    return re;
}

// Return re[i][j]. For the i-th site, measure the bond (i,j).
template <typename Lattice>
vector<int> Measured_bonds (int i, const Lattice& latt)
{
    vector<int> bi;
    for(int j : latt(i).nearest_nb)
    {
        if (j > i)
            bi.push_back (j);
    }
    std::sort(bi.begin(), bi.end());
    return bi;
}

template <typename Lattice>
void Measure_Delta (const SiteSet& sites, const MPS& psi, const Lattice& latt, int i, vector<int>& iv, vector<int>& jv, vector<Real>& Delta)
{
    if (i == length(psi)) return;
    vector<int>  js = Measured_bonds (i, latt);
    vector<Real> Ds = Measure_Delta_ij (sites, i, js, psi);

    for(int k = 0; k < js.size(); k++)
    {
        iv.push_back (i);
        jv.push_back (js.at(k));
        Delta.push_back (Ds.at(k));
    }
}
#endif
