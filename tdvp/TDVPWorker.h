#ifndef __ITENSOR_TDVP_H
#define __ITENSOR_TDVP_H

#include "itensor/iterativesolvers.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"
#include "MyLocalmpo.h"

template <typename TimeType>
void back_propagate
(int ha, int b, int numCenter, TimeType t, MPS& psi, const ITensor& phi1, MyLocalMPO& H, Spectrum& spec, const Args& args)
{
    auto b1 = (ha == 1 ? b+1 : b);

    ITensor phi0;
    if(numCenter == 2)
    {
        phi0 = psi(b1);
    }
    else if(numCenter == 1)
    {
        Index l;
        if(ha == 1) l = commonIndex(psi(b),psi(b+1));
        else        l = commonIndex(psi(b-1),psi(b));
        ITensor U,S,V(l);
        spec = svd(phi1,U,S,V,args);
        psi.ref(b) = U;
        phi0 = S*V;
    }

    H.numCenter(numCenter-1);
    H.position(b1,psi);

    applyExp(H,phi0,t/2,args);

    if(args.getBool("DoNormalize",true))
        phi0 /= norm(phi0);
    
    if(numCenter == 2)
    {
        psi.ref(b1) = phi0;
        psi.leftLim (b1-1);
        psi.rightLim (b1+1);
    }
    if(numCenter == 1)
    {
        if(ha == 1)
        {
            psi.ref(b+1) *= phi0;
            psi.leftLim(b);
            psi.rightLim(b+2);
        }
        else
        {
            psi.ref(b-1) *= phi0;
            psi.leftLim(b-2);
            psi.rightLim(b);
        }
    }
}


template <Direction dir>
void
TDVPWorker(MPS & psi,
           MyLocalMPO& H,
           Cplx t,
           Sweeps const& sweeps,
           DMRGObserver& obs,
           Args args)
{ 
    // Truncate blocks of degenerate singular values (or not)
    args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

    const bool silent = args.getBool("Silent",false);
    if(silent)
    {
        args.add("Quiet",true);
        args.add("PrintEigs",false);
        args.add("NoMeasure",true);
        args.add("DebugLevel",-1);
    }
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? -1 : 0));
    const int numCenter = args.getInt("NumCenter",2);
    if(numCenter != 1)
        args.add("Truncate",args.getBool("Truncate",true));
    else
        args.add("Truncate",args.getBool("Truncate",false));

    const int N = length(psi);

    int oc = orthoCenter (psi);

    args.add("DebugLevel",debug_level);

    // The bonds to sweep
    vector<int> bs;
    if (dir == Fromleft)
    {
        for(int i = oc; i < N; i++)
            bs.push_back (i);
        if (numCenter == 1)
            bs.push_back (N);
    }
    else
    {
        for(int i = oc-numCenter+1; i >= 1; i--)
            bs.push_back (i);
    }

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
    {
        cpu_time sw_time;
        args.add("Sweep",sw);
        args.add("NSweep",sweeps.nsweep());
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("MaxIter",sweeps.niter(sw));

        if(!H.doWrite()
           && args.defined("WriteDim")
           && sweeps.maxdim(sw) >= args.getInt("WriteDim"))
        {
            if(!quiet)
            {
                println("\nTurning on write to disk, write_dir = ",
                        args.getString("WriteDir","./"));
            }
  
            //psi.doWrite(true);
            H.doWrite(true,args);
        }

        // 0, 1 and 2-site wavefunctions
        ITensor phi0,phi1;
        Spectrum spec;
        const int ha = (dir == Fromleft ? 1 : 2);
        for(int b : bs)
        {
            if(!quiet)
                printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));

            H.numCenter(numCenter);
            H.position(b,psi);

            if(numCenter == 2)
                phi1 = psi(b)*psi(b+1);
            else if(numCenter == 1)
                phi1 = psi(b);

            applyExp(H,phi1,-t/2,args);

            if(args.getBool("DoNormalize",true))
                phi1 /= norm(phi1);

            if(numCenter == 2)
                spec = psi.svdBond(b,phi1,(ha==1 ? Fromleft : Fromright),H,args);
            else if(numCenter == 1)
                psi.ref(b) = phi1;

            if((ha == 1 && b+numCenter-1 != N) || (ha == 2 && b != 1))
            {
                back_propagate (ha, b, numCenter, t, psi, phi1, H, spec, args);
            }
 
             if(!quiet)
             { 
                 printfln("    Truncated to Cutoff=%.1E, Min_dim=%d, Max_dim=%d",
                          sweeps.cutoff(sw),
                          sweeps.mindim(sw), 
                          sweeps.maxdim(sw) );
                 printfln("    Trunc. err=%.1E, States kept: %s",
                          spec.truncerr(),
                          showDim(linkIndex(psi,b)) );
             }

            obs.lastSpectrum(spec);
            args.add("AtBond",b);
            args.add("HalfSweep",ha);
            args.add("Truncerr",spec.truncerr());
            obs.measure (args);

        } //for loop over b

        if(!silent)
        {
            auto sm = sw_time.sincemark();
            printfln("    Sweep %d/%d CPU time = %s (Wall time = %s)",
                     sw,sweeps.nsweep(),showtime(sm.time),showtime(sm.wall));
        }
    } //for loop over sw
  
    if(args.getBool("DoNormalize",true))
    {
        psi.normalize();
    }
}

#endif
