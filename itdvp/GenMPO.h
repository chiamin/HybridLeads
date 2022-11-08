#ifndef __GENMPO_H_CMC__
#define __GENMPO_H_CMC__
#include "itensor/mps/mpo.h"
using namespace itensor;
using namespace std;

bool check_Schur (const ITensor& W, const Index& is, const Index& irow, const Index& icol)
{
    assert (irow.dim() == icol.dim());
    int dW = irow.dim();


    for(int a = dW; a >= 1; a--)
    {
        auto Wa = W * setElt (dag(icol)=a);
        auto Waa = Wa * setElt (dag(irow)=a);
        if (a == 2 || a == 1)
        {
            auto I = iut::Identity (Waa.inds());
            if (norm(I-Waa) > 1e-12)
                throw;
        }
        else
        {
            if (norm(Waa) > 1e-12)
                throw;
        }

        for(int b = 1; b <= dW; b++)
        {
            int a2=a, b2=b;
            if (a == 2) a2 = dW;
            else if (a == dW) a2 = 2;
            if (b == 2) b2 = dW;
            else if (b == dW) b2 = 2;

            if (a2 > b2)
            {
                auto Wab = Wa * setElt (dag(irow)=b);
                if (norm(Wab) > 1e-12)
                    throw;
            }
        }
    }
    return true;
}

inline ITensor to_Schur (const ITensor& W1, const Index& is, const Index& irow, const Index& icol)
{
    ITensor W;
    int D = irow.dim();
    assert (irow.dim() == icol.dim());

    auto new_ind = [&D] (int i1)
    {
        if (i1 == 1)
            return i1;
        else if (i1 == 2)
            return D;
        else
            return i1-1;
    };

    for(int i1 = 1; i1 <= D; i1++)
        for(int j1 = 1; j1 <= D; j1++)
        {
            int i = new_ind (i1);
            int j = new_ind (j1);
            auto Wab1 = W1 * setElt (dag(irow)=i1) * setElt (dag(icol)=j1);
            auto Wab  = Wab1 * setElt (irow=i) * setElt (icol=j);
            W += Wab;
        }
    check_Schur (W, is, irow, icol);
    return W;
}

inline tuple <ITensor,Index,Index,Index> to_Schur (const MPO& H)
{
    auto W = H(2);
    auto irow = commonIndex (H(2),H(1));
    auto icol = commonIndex (H(2),H(3));
    auto is = findIndex (W, "Site,0");
    auto W_re = to_Schur (W, is, irow, icol);
    return {W_re, is, irow, icol};
}

inline tuple <ITensor,Index,Index,Index> get_MPO_Tensor_TFI (Real gx)
{
    auto sites = SpinHalf (3,{"ConserveSz",false});
    AutoMPO ampo (sites);
    for(int j = 1; j < 3; ++j)
    {
        ampo += -4,"Sz",j,"Sz",j+1;
        ampo += 2*gx,"Sx",j;
    }
    ampo += 2*gx,"Sx",3;
    auto H = toMPO (ampo);

    auto W = H(2);
    auto iwl = commonIndex (H(2),H(1));
    auto iwr = commonIndex (H(2),H(3));
    auto is = findIndex (W, "Site,0");
    return {W, is, iwl, iwr};
}

inline tuple <ITensor,Index,Index,Index> get_W (const MPO& H)
{
    auto W = H(2);
    auto iwl = commonIndex (H(2),H(1));
    auto iwr = commonIndex (H(2),H(3));
    auto is = findIndex (W, "Site,0");
    check_Schur (W, is, iwl, iwr);
    return {W, is, iwl, iwr};
}

tuple <ITensor,Index,Index,Index> get_MPO_Tensor_Hubbard_spinless (Real t, Real mu, const Args& args=Args::global())
{
    auto sites = Spinless (3, args);
    AutoMPO ampo (sites);
    for(int i = 1; i <= 2; i++)
    {
        ampo += -t,"Cdag",i,"C",i+1;
        ampo += -t,"Cdag",i+1,"C",i;
        ampo += -mu,"N",2;
    }
    ampo += -mu,"N",3;
    auto H = toMPO (ampo);
    return get_W (H);
}

tuple <ITensor,Index,Index,Index> get_MPO_Tensor_XXZ_spin1 (Real Delta, const Args& args=Args::global())
{
    auto sites = SpinOne (3, args)  ;
    AutoMPO ampo (sites);
    for(int i = 1; i <= 2; i++)
    {
        ampo += 0.5,"S+",i,"S-",i+1;
        ampo += 0.5,"S+",i+1,"S-",i;
        ampo += Delta,"Sz",i,"Sz",i+1;
    }
    auto H = toMPO (ampo);
    return get_W (H);
}

void to_inf_mpo (MPO& mpo, Index iL=Index(), Index iR=Index())
{
    int N = length (mpo);
    if (N < 4)
    {
        cout << "Error: " << __FUNCTION__ << ": MPO length must >= 4" << endl;
        throw;
    }

    auto is1 = findIndex (mpo(1), "Site,0");
    auto is2 = findIndex (mpo(2), "Site,0");
    auto i1r = commonIndex (mpo(1), mpo(2));
    auto i2r = commonIndex (mpo(2), mpo(3));
    auto i2l = dag(i1r);
    if (!iL) iL = sim(i2l);

    mpo.ref(1) = replaceInds (mpo(2), {is2, prime(is2), i2l, i2r}, {is1, prime(is1), iL, i1r});
    assert (commonIndex (mpo(1), mpo(2)));

    auto isN  = findIndex (mpo(N), "Site,0");
    auto isN2 = findIndex (mpo(N-1), "Site,0");
    auto iNl  = commonIndex (mpo(N), mpo(N-1));
    auto iN2l = commonIndex (mpo(N-1), mpo(N-2));
    auto iN2r = dag(iNl);
    if (!iR) iR = sim(iN2r);

    mpo.ref(N) = replaceInds (mpo(N-1), {isN2, prime(isN2), iN2l, iN2r}, {isN, prime(isN), iNl, iR});
    assert (commonIndex (mpo(N-1), mpo(N)));
}

MPO single_empurity_mpo (const SiteSet& sites, int imp_site, Real t, Real t_impL, Real t_impR, Real muL, Real muR, Real mu_imp,
                         Real V, Real V_impL, Real V_impR)
{
    auto N = length (sites);
    AutoMPO ampo (sites);
    for(int i = 1; i <= N; ++i)
    {
        Real mui, ti, Vi;
        // Set mu
        if      (i < imp_site)  mui = muL;
        else if (i > imp_site)  mui = muR;
        else                    mui = mu_imp;
        // Set t
        if      (i == imp_site-1)   ti = t_impL;
        else if (i == imp_site)     ti = t_impR;
        else                        ti = t;
        // Set V
        if      (i == imp_site-1)   Vi = V_impL;
        else if (i == imp_site)     Vi = V_impR;
        else                        Vi = V;

        if (i != N)
        {
            ampo += -ti,"Cdag",i,"C",i+1;
            ampo += -ti,"Cdag",i+1,"C",i;
        }
        if (mui != 0.)
            ampo += -mui,"N",i;
        if (Vi != 0. && i != N)
            ampo += Vi,"N",i,"N",i+1;
    }
    auto H = toMPO (ampo);

    return H;
}

tuple <AutoMPO,int,int>
t_mu_V_ampo
(const SiteSet& sites, int L_device,
 Real t_leadL, Real t_leadR, Real t_device, Real t_contactL, Real t_contactR,
 Real mu_leadL, Real mu_leadR, Real mu_device,
 Real V_leadL, Real V_leadR, Real V_device, Real V_contactL, Real V_contactR)
{
    auto N = length (sites);
    int first_site_device = (N - L_device) / 2 + 1;
    int last_site_device = first_site_device + L_device - 1;
    AutoMPO ampo (sites);
    for(int i = 1; i <= N; ++i)
    {
        Real mui, ti, Vi;
        // Set mu
        if      (i < first_site_device)  mui = mu_leadL;
        else if (i > last_site_device)   mui = mu_leadR;
        else                             mui = mu_device;
        // Set t and V
        if (i == first_site_device-1)
        {
            ti = t_contactL;
            Vi = V_contactL;
        }
        else if (i == last_site_device)
        {
            ti = t_contactR;
            Vi = V_contactR;
        }
        else if (i < first_site_device-1)
        {
            ti = t_leadL;
            Vi = V_leadL;
        }
        else if (i > last_site_device)
        {
            ti = t_leadR;
            Vi = V_leadR;
        }
        else
        {
            ti = t_device;
            Vi = V_device;
        }

        if (i != N)
        {
            ampo += -ti,"Cdag",i,"C",i+1;
            ampo += -ti,"Cdag",i+1,"C",i;
cout << "t: " << i << " " << i+1 << " " << ti << endl;
        }
        if (mui != 0.)
        {
            ampo += -mui,"N",i;
cout << "mu: " << i << " " << mui << endl;
        }
        if (Vi != 0. && i != N)
        {
            ampo += Vi,"N",i,"N",i+1;
cout << "V: " << i << " " << i+1 << " " << Vi << endl;
        }
    }
    return {ampo, first_site_device, last_site_device};
}

void ampo_add_mu (AutoMPO& ampo, const vector<tuple<int,Real>>& site_mu, bool verbose=false)
{
    if (verbose)
        cout << "i mu" << endl;
    for(auto [site, mu] : site_mu)
    {
        if (mu != 0.)
        {
            ampo += -mu,"N",site;
            if (verbose)
                cout << site << " " << mu << endl;
        }
    }
}

void ampo_add_tV (AutoMPO& ampo, const vector<tuple<int,int,Real,Real>>& link_t_V, bool verbose=false)
{
    if (verbose)
        cout << "i j  t V" << endl;
    for(auto [i, j, t, V] : link_t_V)
    {
        if (t != 0.)
        {
            ampo += -t,"Cdag",i,"C",j;
            ampo += -t,"Cdag",j,"C",i;
        }
        if (V != 0.)
        {
            ampo += V,"N",i,"N",j;
            ampo += V,"N",j,"N",i;
        }
        if (verbose)
            cout << i << " " << j << " " << t << " " << V << endl;
    }
}
#endif
