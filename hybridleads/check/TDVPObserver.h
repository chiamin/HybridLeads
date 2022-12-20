#ifndef __TDVPOBSERVER_H_CMC__
#define __TDVPOBSERVER_H_CMC__
#include <iomanip>
#include <map>

#include "ContainerUtility.h"
#include "Entanglement.h"
#include "itensor/all.h"
using namespace vectool;
using namespace iutility;

template <typename SitesType>
class TDVPObserver : public DMRGObserver {
 public:
  TDVPObserver(const SitesType& sites, const MPS& psi,
               const Args& args = Args::global())
      : DMRGObserver(psi, args),
        _sites(sites),
        _ns(length(psi), 0.),
        _Npar(0.),
        _specs(length(psi)) {
    _write = args.getBool("Write", false);
    _out_dir = args.getString("out_dir", ".");
  }

  void measure(const Args& args);

  Real Npar() const { return _Npar; }
  auto const& ns() const { return _ns; }
  const Spectrum& spec(int i) const { return _specs.at(i); }

 private:
  bool _write;
  string _out_dir;  // empty string "" if not write
  SitesType _sites;

  // Observables
  vector<Real> _ns;
  Real _Npar;
  vector<Spectrum> _specs;
};

template <typename SitesType>
void TDVPObserver<SitesType>::measure(const Args& args) {
  DMRGObserver::measure(args);

  cout << scientific << setprecision(14);
  // Define your measurements below
  // Call psi() to access the MPS
  //
  auto N = length(psi());
  auto b = args.getInt("AtBond");
  auto sw = args.getInt("Sweep");
  auto ha = args.getInt("HalfSweep");
  auto energy = args.getReal("Energy", 0);

  if (b != N) _specs.at(b) = spectrum();

  int oc = orthoCenter(psi());
  int nc = args.getInt("NumCenter");
  // measure during the second half of sweep
  if ((nc == 2 && oc == N) || ha == 2) {
    // Density
    ITensor n_tmp = psi().A(oc) * _sites.op("N", oc);
    n_tmp.noPrime("Site");
    n_tmp *= dag(psi().A(oc));
    Real ni = toReal(n_tmp);
    cout << "\tn " << oc << " = " << ni << endl;
    _ns.at(oc - 1) = ni;

    // Entanglement entropy
    Real S = EntangEntropy(spectrum());
    cout << "\tentang entropy " << oc << " = " << S << endl;
  }

  // At the end of a sweep
  if (oc == 1 && ha == 2 && b == 1) {
    for (int i = 1; i < N; i++) {
      cout << "\tbond dim " << i << " = " << dim(rightLinkIndex(psi(), i))
           << endl;
    }

    if (_write) {
      cout << "write MPS" << endl;
      writeToFile(_out_dir + "/psi.mps", psi());
    }
  }
}
#endif
