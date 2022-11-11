#ifndef __SQUARELATTICE_H_CMC__
#define __SQUARELATTICE_H_CMC__
#include <vector>
using namespace std;

struct LatticeSite
{
    int         index;
    int         x,y;
    vector<int> nearest_nb;
};

class SquareLattice
{
    public:
        SquareLattice (int lx, int ly, bool periodic_x, bool periodic_y);

              int          N       ()         const { return _N; }
        const LatticeSite& operator() (int i) const { return _sites.at(i); }

        // x and y start from 1
        // Returned index starts from 1
        // (1,1)=1, (1,2)=2,...
        static int mix_ind (int x, int y, int Lx, int Ly) { return (x-1)*Ly + y; }
        // If x in the lattice, return the adjusted site
        // else, return -1
        static int apply_bc (int x, int L, bool pbc)
        {
            if (x > L)
            {
                if (pbc) x -= L;
                else x = -1;
            }
            if (x <= 0)
            {
                if (pbc) x += L;
                else x = -1;
            }
            return x;
        }

    private:
        int                 _lx, _ly, _N, _periodic_x, _periodic_y;
        vector<LatticeSite> _sites;
};

SquareLattice :: SquareLattice (int lx, int ly, bool periodic_x, bool periodic_y)
: _lx (lx)
, _ly (ly)
, _N (lx*ly)
, _periodic_x (periodic_x)
, _periodic_y (periodic_y)
, _sites (lx*ly+1)
{
    for(int x = 1; x <= lx; x++)
    {
        int xp = x+1,
            xn = x-1;
        for(int y = 1; y <= ly; y++)
        {
            int yp = y+1,
                yn = y-1;

            int i = SquareLattice::mix_ind (x,y,lx,ly);
            _sites.at(i).index = i;
            _sites.at(i).x = x;
            _sites.at(i).y = y;

            // Set boundary neighbors
            xp = SquareLattice::apply_bc (xp, lx, periodic_x);
            xn = SquareLattice::apply_bc (xn, lx, periodic_x);
            yp = SquareLattice::apply_bc (yp, ly, periodic_y);
            yn = SquareLattice::apply_bc (yn, ly, periodic_y);

            // Set nearest-neighbor indices
            if (xp != -1)
            {
                int j = SquareLattice::mix_ind (xp,y,lx,ly);
                _sites.at(i).nearest_nb.push_back (j);
            }
            if (xn != -1)
            {
                int j = SquareLattice::mix_ind (xn,y,lx,ly);
                _sites.at(i).nearest_nb.push_back (j);
            }
            if (yp != -1)
            {
                int j = SquareLattice::mix_ind (x,yp,lx,ly);
                _sites.at(i).nearest_nb.push_back (j);
            }
            if (yn != -1)
            {
                int j = SquareLattice::mix_ind (x,yn,lx,ly);
                _sites.at(i).nearest_nb.push_back (j);
            }
        }
    }
}
#endif
