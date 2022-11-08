#ifndef __MYOBSERVER_H_CMC__
#define __MYOBSERVER_H_CMC__
#include <iomanip>
#include "itensor/all.h"
#include "Entanglement.h"

class MyObserver : public DMRGObserver
{
    public:
        MyObserver (const Fermion& sites, const MPS& psi, const Args& args = Args::global())
        : DMRGObserver (psi, args)
        , _sites (sites)
        , _ns (length(psi)+1,0.)
        , _Npar (0.)
        , _sites_obs (2, args)
        , _specs (length(psi))
        {
            _write = args.getBool ("Write",false);
            _out_dir = args.getString("out_dir",".");

            // Current operator
            AutoMPO ampo (_sites_obs);
            ampo += -2_i,"Cdag",1,"C",2;
            auto mpo = toMPO (ampo);
            _current_op = mpo(1) * mpo(2);
        }

        void measure (const Args& args);

        Real Npar () const { return _Npar; }
        const Spectrum& spec (int i) const { return _specs.at(i); }

    private:
        bool                _write;
        string              _out_dir;	// empty string "" if not write
        Fermion             _sites;

        vector<Real>        _ns;
        Real                _Npar;

        // Observables
        Fermion             _sites_obs;
        ITensor             _current_op;
        vector<Spectrum>    _specs;
};

inline Real Onsite_mea (const ITensor& A, const ITensor& op)
{
    ITensor re = A * op;
    re.noPrime ("Site");
    re *= dag(A);
    return iut::toReal (re);
}

void MyObserver :: measure (const Args& args)
{
    DMRGObserver::measure (args);

    cout << scientific << setprecision(14);
    // Define your measurements below
    // Call psi() to access the MPS
    //
    auto N = length(psi());
    auto b = args.getInt("AtBond");
    auto sw = args.getInt("Sweep");
    auto ha = args.getInt("HalfSweep");
    auto energy = args.getReal("Energy",0);

    if (b != N)
        _specs.at(b) = spectrum();

    int oc = orthoCenter(psi());
    int nc = args.getInt("NumCenter");
    // measure during the second half of sweep
    if ((nc == 2 && oc == N) || ha == 2)
    {
        // Density
        ITensor n_op = _sites.op("N",oc);
        Real ni = Onsite_mea (psi().A(oc), n_op);
        cout << "\tn " << oc << " = " << ni << endl;

        // Current
        if (oc != N)
        {
            auto const& i1 = _sites_obs (1);
            auto const& i2 = _sites_obs (2);
            auto const& j1 = _sites (oc);
            auto const& j2 = _sites (oc+1);
            auto Jop = replaceInds (_current_op, {i1, prime(i1), i2, prime(i2)}, {j1, prime(j1), j2, prime(j2)});

            auto phi = psi()(oc) * psi()(oc+1);
            auto phiJ = noPrime (phi * Jop, "Site");
            Cplx j = eltC (phiJ * dag(phi));
            cout << "\tcurrent " << oc << " " << oc+1 << " = " << j << endl;
        }

        // Entanglement entropy
        Real S = iut::EntangEntropy (spectrum());
        cout << "EE " << oc << " = " << S << endl;
    }

    if (oc == 1 && ha == 2)
    {
        if (_write)
        {
            cout << "write MPS" << endl;
            writeToFile (_out_dir+"/psi.mps", psi());
        }
    }
}

#endif
