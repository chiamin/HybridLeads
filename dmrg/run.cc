#include "itensor/all.h"
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "../itdvp/GlobalIndices.h"
#include "GenMPO.h"
using namespace itensor;
using namespace std;

template <typename SitesType>
MPO current_correlation (const SitesType& sites, int i)
{
    int N = length(sites);
    AutoMPO ampo (sites);
    if constexpr (is_same_v <SitesType, Fermion>)
    {
        ampo += -1_i,"Cdag",i,"C",i+1;
        ampo +=  1_i,"Cdag",i+1,"C",i;
        ampo +=  1_i,"Cdag",N-i,"C",N-i+1;
        ampo += -1_i,"Cdag",N-i+1,"C",N-i;
    }

    return toMPO (ampo);
}

Cplx myinner (const MPO& mpo, const MPS& psi)
{
    auto L = ITensor(1.);
    int N = length(psi);
    for(int i = 1; i <= N; i++)
    {
        L *= psi(i) * mpo(i);
        auto dagA = dag(psi(i));
        dagA.prime ("Site");
        if (i == 1)
        {
            auto ii = commonIndex (psi(i), psi(i+1), "Link");
            dagA.prime (ii);
        }
        else if (i == N)
        {
            auto ii = commonIndex (psi(i), psi(i-1), "Link");
            dagA.prime (ii);
        }
        else
        {
            dagA.prime ("Link");
        }
        L *= dagA;
    }
    return eltC(L);
}

MPS toMPS (const SiteSet& sites,
           const ITensor& AL, const Index& iALl, const Index& iALr, const Index& iALs,
           const ITensor& AR, const Index& iARl, const Index& iARr, const Index& iARs,
           const ITensor& AC, const Index& iACl, const Index& iACr, const Index& iACs,
           int oc=1)
{
    MPS psi (sites);
    int N = length (psi);

    // Set tensors
    vector<Index> ils(N+1), irs(N+1), iss(N+1);
    for(int i = 1; i <= N; i++)
    {
        Index il, ir;
        if (i < oc)
        {
            psi.ref(i) = AL;
            ils.at(i) = iALl;
            irs.at(i) = iALr;
            iss.at(i) = iALs;
        }
        else if (i == oc)
        {
            psi.ref(i) = AC;
            ils.at(i) = iACl;
            irs.at(i) = iACr;
            iss.at(i) = iACs;
        }
        else
        {
            psi.ref(i) = AR;
            ils.at(i) = iARl;
            irs.at(i) = iARr;
            iss.at(i) = iARs;
        }
        assert (hasIndex (psi(i), ils.at(i)));
        assert (hasIndex (psi(i), irs.at(i)));
    }

    // Replace indices
    auto inew = noPrime (sim (irs.at(1)));
    psi.ref(1).replaceInds ({iss.at(1), irs.at(1)}, {sites(1), inew});
    irs.at(1) = inew;
    for(int i = 2; i < N; i++)
    {
        inew = noPrime(sim(irs.at(i)));
        psi.ref(i).replaceInds ({iss.at(i), ils.at(i), irs.at(i)}, {sites(i), irs.at(i-1), inew});
        ils.at(i) = irs.at(i-1);
        irs.at(i) = inew;
    }
    psi.ref(N).replaceInds ({iss.at(N), ils.at(N)}, {sites(N), irs.at(N-1)});
    ils.at(N) = irs.at(N-1);

    // Check
    for(int i = 1; i <= N; i++)
    {
        assert (hasIndex (psi(i), sites(i)));
        if (i != N)
            assert (commonIndex (psi(i), psi(i+1)));
    }

    return psi;
}

int main(int argc, char* argv[])
{
    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto L          = input.getInt("L");
    auto imp_site   = input.getInt("imp_site");
    auto t          = input.getReal("t");
    auto t_impL     = input.getReal("t_impL");
    auto t_impR     = input.getReal("t_impR");
    auto muL        = input.getReal("muL");
    auto muR        = input.getReal("muR");
    auto mu_imp     = input.getReal("mu_imp");
    auto V          = input.getReal("V");
    auto V_impL     = input.getReal("V_impL");
    auto V_impR     = input.getReal("V_impR");
    auto read_dir_L = input.getString("read_dir_L");
    auto read_dir_R = input.getString("read_dir_R");

    auto WriteDim   = input.getInt("WriteDim",-1);
    auto do_write   = input.getYesNo("write_to_file");
    auto out_dir    = input.getString("outdir",".");
    auto out_minm   = input.getInt("out_minm",0);
    auto H_file     = input.getString("H_outfile","H.mpo");
    auto ConserveQNs = input.getYesNo("ConserveQNs",false);
    auto sweeps     = Read_sweeps (infile);

    // Read the tensors
    ITensor AL_L, AR_L, AC_L, LW_L, RW_L;
    ITensor AL_R, AR_R, AC_R, LW_R, RW_R;
    readFromFile (read_dir_L+"/AL.itensor", AL_L);
    readFromFile (read_dir_L+"/AR.itensor", AR_L);
    readFromFile (read_dir_L+"/AC.itensor", AC_L);
    readFromFile (read_dir_L+"/LW.itensor", LW_L);
    readFromFile (read_dir_L+"/RW.itensor", RW_L);
    readFromFile (read_dir_R+"/AL.itensor", AL_R);
    readFromFile (read_dir_R+"/AR.itensor", AR_R);
    readFromFile (read_dir_R+"/AC.itensor", AC_R);
    readFromFile (read_dir_R+"/LW.itensor", LW_R);
    readFromFile (read_dir_R+"/RW.itensor", RW_R);
    GlobalIndices ISL, ISR;
    ISL.read (read_dir_L+"/global.inds");
    ISR.read (read_dir_R+"/global.inds");

    // Site set
    using SitesType = Fermion;
    auto sites = SitesType (L, {"ConserveQNs",ConserveQNs});

    // Initialze MPS
    MPS psi = toMPS (sites, AL_L, ISL.il(),          prime(ISL.il(),2), ISL.is(),
                            AR_R, prime(ISR.ir(),2), ISR.ir(),          ISR.is(),
                            AC_L, ISL.il(),          ISL.ir(),          ISL.is());
    psi.position(1);
    assert (commonIndex (LW_L, psi(1)));
    assert (commonIndex (RW_R, psi(L)));

    // Make MPO
    MPO H = single_empurity_mpo (sites, imp_site, t, t_impL, t_impR, muL, muR, mu_imp, V, V_impL, V_impR);
    to_inf_mpo (H, ISL.iwl(), ISR.iwr());
    assert (commonIndex (LW_L, H(1)));
    assert (commonIndex (RW_R, H(L)));

    // Write to files
    writeToFile (out_dir+"/"+H_file, H);

    // DMRG
    MyObserver<SitesType> myobs (sites, psi, {"Write",do_write,"out_dir",out_dir,"out_minm",out_minm});
    dmrg (psi, H, LW_L, RW_R, sweeps, myobs, {"WriteDim",WriteDim});
/*
    int N = length(sites);
    for(int i = 1; i <= N/2; i++)
    {
        auto J = current_correlation (sites, i);
        auto j = myinner (J, psi);
        cout << "<J,J> " << i << " = " << j << endl;
    }
*/
    return 0;
}
