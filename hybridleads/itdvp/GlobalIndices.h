#ifndef __GLOBALINDCIES_H_CMC__
#define __GLOBALINDCIES_H_CMC__
#include "itensor/all.h"
using namespace itensor;

enum enumLR { LEFT, RIGHT };

struct GlobalIndices {
  const Index& il() const { return _il; }
  const Index& ir() const { return _ir; }
  const Index& is() const { return _is; }
  const Index& iwl() const { return _iwl; }
  const Index& iwr() const { return _iwr; }

  void write(const string& filename) const {
    ofstream ofs(filename);
    write(ofs);
  }
  void write(ostream& ofs) const {
    itensor::write(ofs, _is);
    itensor::write(ofs, _il);
    itensor::write(ofs, _ir);
    itensor::write(ofs, _iwl);
    itensor::write(ofs, _iwr);
  }
  void read(const string& filename) {
    ifstream ifs(filename);
    if (!ifs) {
      cout << __FUNCTION__ << ": Cannot open file: " << filename << endl;
      throw;
    }
    read(ifs);
  }
  void read(istream& ifs) {
    itensor::read(ifs, _is);
    itensor::read(ifs, _il);
    itensor::read(ifs, _ir);
    itensor::read(ifs, _iwl);
    itensor::read(ifs, _iwr);
  }

  vector<Index> LW_inds() const { return {_il, prime(_il), _iwl}; }
  vector<Index> RW_inds() const { return {_ir, prime(_ir), _iwr}; }

  bool check(string name, const ITensor& T) const {
    // assert (isReal (T));
    if (name == "A") {
      assert(order(T) == 3);
      assert(hasIndex(T, _is));
      assert(hasIndex(T, _il));
      assert(hasIndex(T, prime(_il, 2)));
    } else if (name == "W") {
      assert(order(T) == 4);
      assert(hasIndex(T, _is));
      assert(hasIndex(T, prime(_is)));
      assert(hasIndex(T, _iwl));
      assert(hasIndex(T, _iwr));
    } else if (name == "AL") {
      assert(order(T) == 3);
      assert(hasIndex(T, _is));
      assert(hasIndex(T, _il));
      assert(hasIndex(T, prime(_il, 2)));
    } else if (name == "AR") {
      assert(order(T) == 3);
      assert(hasIndex(T, _is));
      assert(hasIndex(T, _ir));
      assert(hasIndex(T, prime(_ir, 2)));
    } else if (name == "AC") {
      assert(order(T) == 3);
      assert(hasIndex(T, _is));
      assert(hasIndex(T, _il));
      assert(hasIndex(T, _ir));
    } else if (name == "C") {
      assert(order(T) == 2);
      assert(hasIndex(T, prime(_il, 2)));
      assert(hasIndex(T, prime(_ir, 2)));
    } else if (name == "R") {
      assert(order(T) == 2);
      assert(hasIndex(T, _il));
      assert(hasIndex(T, prime(_il)));
    } else if (name == "L") {
      assert(order(T) == 2);
      assert(hasIndex(T, _ir));
      assert(hasIndex(T, prime(_ir)));
    } else if (name == "LW") {
      assert(order(T) == 3);
      assert(hasIndex(T, _il));
      assert(hasIndex(T, prime(_il)));
      assert(hasIndex(T, _iwl));
    } else if (name == "RW") {
      assert(order(T) == 3);
      assert(hasIndex(T, _ir));
      assert(hasIndex(T, prime(_ir)));
      assert(hasIndex(T, _iwr));
    }
    return true;
  }

  Index _il, _ir, _is, _iwl, _iwr;
};

namespace global {
auto IS = GlobalIndices();
}
#endif
