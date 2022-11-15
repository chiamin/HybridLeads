#ifndef __MYOBSERVER_H_CMC__
#define __MYOBSERVER_H_CMC__
#include "itensor/all.h"
using namespace itensor;

template <typename SitesType>
class MyObserver : public DMRGObserver
{
    public:
        MyObserver (const SitesType& sites, const MPS& psi, const Args& args = Args::global())
        : DMRGObserver (psi, args)
        , _sites (sites)
        , _ns (length(psi)+1,0.)
        , _Npar (0.)
        {
            _write = args.getBool ("Write",false);
            _write_minm = args.getInt ("out_minm",0);
            _out_dir = args.getString("out_dir",".");
            _ConserveNf = args.getBool ("ConserveNf",true);
        }

        void measure (const Args& args = Args::global());

        Real Npar () const { return _Npar; }

    private:
        bool                 _write;
        string               _out_dir;	// empty string "" if not write
        int                  _write_minm;
        vector<int>          _iDel, _jDel;
        vector<Real>         _Delta;
        bool                 _ConserveNf;
        SitesType            _sites;

        vector<Real>         _ns;
        Real                 _Npar;
};

Real Onsite_mea (const ITensor& A, const ITensor& op)
{
    ITensor re = A * op;
    re.noPrime ("Site");
    re *= dag(A);
    return re.real();
}

template <typename SitesType>
void MyObserver<SitesType> :: measure (const Args& args)
{
    DMRGObserver::measure (args);

    // Define your measurements below
    // Call psi() to access the MPS
    //
    auto N = length(psi());
    auto b = args.getInt("AtBond",1);
    auto sw = args.getInt("Sweep",0);
    auto ha = args.getInt("HalfSweep",0);
    auto energy = args.getReal("Energy",0);

    // On-site measure
    int oc = orthoCenter(psi());
/*    if (oc == N || ha == 2) // measure during the second half of sweep
    {
        ITensor nup_op = _sites.op("Nup",oc),
                ndn_op = _sites.op("Ndn",oc);
        Real n_up = Onsite_mea (psi().A(oc), nup_op),
             n_dn = Onsite_mea (psi().A(oc), ndn_op);

        Real ni = n_up + n_dn,
             sz = 0.5*(n_up-n_dn);

        _Npar += ni - _ns.at(oc);
        _ns.at(oc) = ni;
    
        cout << "    Measure on site " << oc << endl;
        cout << "      n_up = " << n_up << endl;
        cout << "      n_dn = " << n_dn << endl;
        cout << "      N_tot = " << ni << endl;
        cout << "      sz = " << sz << endl;
    }*/

    if (oc == 1 && ha == 2)
    {
        int m = args.getInt("MaxDim");
        // Write MPS
        if (_write && m >= _write_minm)
        {
            writeToFile (_out_dir+"/psi"+"_m"+to_string(m)+".mps", psi());
        }
    }
}

#endif
