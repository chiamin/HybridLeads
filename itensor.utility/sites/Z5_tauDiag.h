//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_Z5_tauDiag_H
#define __ITENSOR_Z5_tauDiag_H
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class Z5_tauDiagSite;

using Z5_tauDiag = BasicSiteSet<Z5_tauDiagSite>;

class Z5_tauDiagSite
    {
    Index s;
    public:

    Z5_tauDiagSite(Index I) : s(I) { }

    Z5_tauDiagSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,Z5");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        if(args.getBool("ConserveQNs",true))
          {
          s = Index(QN({"T",0,5}),1,
                    QN({"T",1,5}),1,
                    QN({"T",2,5}),1,
                    QN({"T",3,5}),1,
                    QN({"T",4,5}),1,
                    Out,ts);
          }
        else
          {
          s = Index(5,ts);
          }
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "0") { return s(1); }
        else
        if(state == "1") { return s(2); }
        else
        if(state == "2") { return s(3); }
        else
        if(state == "3") { return s(4); }
        else
        if(state == "4") { return s(5); }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Zer = s(1);
        auto ZerP = sP(1);
        auto One = s(2);
        auto OneP = sP(2);
        auto Two = s(3);
        auto TwoP = sP(3);
        auto Thr = s(4);
        auto ThrP = sP(4);
        auto Fou = s(5);
        auto FouP = sP(5);

        auto Op = ITensor(dag(s),sP);

        if(opname == "N")
            {
            Op.set(One,OneP,1);
            Op.set(Two,TwoP,2);
            Op.set(Thr,ThrP,3);
            Op.set(Fou,FouP,4);
            }
        else
        if(opname == "Sig")
            {
            Op.set(Zer,FouP,1);
            Op.set(One,ZerP,1);
            Op.set(Two,OneP,1);
            Op.set(Thr,TwoP,1);
            Op.set(Fou,ThrP,1);
            }
        else
        if(opname == "SigSqr")
            {
            Op.set(Zer,ThrP,1);
            Op.set(One,FouP,1);
            Op.set(Two,ZerP,1);
            Op.set(Thr,OneP,1);
            Op.set(Fou,TwoP,1);
            }
        else
        if(opname == "SigDag")
            {
            Op.set(Fou,ZerP,1);
            Op.set(Zer,OneP,1);
            Op.set(One,TwoP,1);
            Op.set(Two,ThrP,1);
            Op.set(Thr,FouP,1);
            }
        else
        if(opname == "Tau")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/5.)+sin(2.*Pi/5.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/5.)+sin(4.*Pi/5.)*1_i);
            Op.set(Thr,ThrP,cos(6.*Pi/5.)+sin(6.*Pi/5.)*1_i);
            Op.set(Fou,FouP,cos(8.*Pi/5.)+sin(8.*Pi/5.)*1_i);
            }
        else
        if (opname == "LogTau")
            {
            Real a = 1.25663706;
            Real b = 2.51327412;
            Op.set(Zer,ZerP, 0);
            Op.set(One,OneP,-a);
            Op.set(Two,TwoP,-b);
            Op.set(Thr,ThrP, b);
            Op.set(Fou,FouP, a);
            }
        else
        if(opname == "TauSqr")
            {
            Op.set(Zer,ZerP,1);
            Cplx c1 = cos(2.*Pi/5.)+sin(2.*Pi/5.)*1_i;
            Cplx c2 = cos(4.*Pi/5.)+sin(4.*Pi/5.)*1_i;
            Cplx c3 = cos(6.*Pi/5.)+sin(6.*Pi/5.)*1_i;
            Cplx c4 = cos(8.*Pi/5.)+sin(8.*Pi/5.)*1_i;
            Op.set(One,OneP,c1*c1);
            Op.set(Two,TwoP,c2*c2);
            Op.set(Thr,ThrP,c3*c3);
            Op.set(Fou,FouP,c4*c4);
            }
        else
        if(opname == "TauDag")
            {
            Op.set(Zer,ZerP,1);
            Op.set(One,OneP,cos(2.*Pi/5.)-sin(2.*Pi/5.)*1_i);
            Op.set(Two,TwoP,cos(4.*Pi/5.)-sin(4.*Pi/5.)*1_i);
            Op.set(Thr,ThrP,cos(6.*Pi/5.)-sin(6.*Pi/5.)*1_i);
            Op.set(Fou,FouP,cos(8.*Pi/5.)-sin(8.*Pi/5.)*1_i);
            }
        else
        if(opname == "Proj0")
            {
            Op.set(Zer,ZerP,1);
            }
        else
        if(opname == "Proj1")
            {
            Op.set(One,OneP,1);
            }
        else
        if(opname == "Proj2")
            {
            Op.set(Two,TwoP,1);
            }
        else
        if(opname == "Proj3")
            {
            Op.set(Thr,ThrP,1);
            }
        else
        if(opname == "Proj4")
            {
            Op.set(Fou,FouP,1);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    Z5_tauDiagSite(int n, Args const& args = Args::global())
        {
        *this = Z5_tauDiagSite({args,"SiteNumber=",n});
        }

    };

} //namespace itensor

#endif
