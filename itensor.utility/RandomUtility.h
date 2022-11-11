#ifndef __GANDGAN_H_CMC__
#define __GANDGAN_H_CMC__
#include <random>

class RandGen
{
    public:
        using RealType = double;
        using SeedType = unsigned long;
        using UnifReal = std::uniform_real_distribution<>;

        RandGen (SeedType seed_=0, double min=0., double max=1.)
        : gen (rd())
        , dis (std::make_unique<UnifReal>(min, max))
        {
            if (seed_) seed_real (seed_);
        }

        RealType real () { return dis->operator()(gen); }

        void range_real  (RealType min, RealType max) { dis = std::make_unique<UnifReal>(min,max); }
        void seed_real   (SeedType value)             { gen.seed (value); }

    private:
         std::random_device             rd;
         std::mt19937_64                gen;
         std::unique_ptr<UnifReal>      dis;
};
#endif
