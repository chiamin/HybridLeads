#ifndef __PERMUT_TABLE_H_CMC__
#define __PERMUT_TABLE_H_CMC__
#include <vector>
#include "ContainerUtility.h"
#include "CountSwap.h"
using namespace std;

tuple <vector<vector<int>>, vector<int>> get_permutations (const vector<int>& pos)
// <pos> is assumed to be sorted
// Return 1) all the permutations of pos. 2) the number of neighboring swaps required
// The size of pos can be 2, 3, or 4
// Did not care about the efficiency. Can be slow.
{
    assert (pos.size() >=2 && pos.size() <= 4);

    vector<vector<int>> re;
    vector<int> swaps;

    auto add_permute = [&re, &swaps] (const vector<int>& ppos)
    {
        if (!in_vector (re, ppos))
        {
            re.push_back (ppos);
            int nswap = countSwaps (ppos);
            swaps.push_back (nswap);
        }
    };

    for(int p1 : pos)
    {
        auto tmp2 = pos;
        remove_element (tmp2, p1);
        for(int p2 : tmp2)
        // p2 != p1
        {
            if (pos.size() == 2)
            // Two positions
            {
                add_permute ({p1,p2});
            }
            else
            {
                auto tmp3 = tmp2;
                remove_element (tmp3, p2);
                for(int p3 : tmp3)
                // p3 != p1,p2
                {
                    if (pos.size() == 3)
                    // Three positions
                    {
                        add_permute ({p1,p2,p3});
                    }
                    else
                    // Four positions
                    {
                        auto tmp4 = tmp3;
                        remove_element (tmp4, p3);
                        for(int p4 : tmp4)
                        // p4 != p1,p2,p3
                        {
                            add_permute ({p1,p2,p3,p4});
                        }
                    }
                }
            }
        }
    }
    return {re, swaps};
}

class PermutTable
{
    public:
        PermutTable () {}

        // Important: The input <pos> is assumed to be sorted.
        tuple <vector<vector<int>>, vector<int>> permut (const vector<int>& pos);

    private:
        // Key: site positions
        // Value: 1) all the permutations, 2) corresponding number of swaps
        map <vector<int>, tuple <vector<vector<int>>, vector<int>>> _ptable;
};

tuple <vector<vector<int>>, vector<int>> PermutTable :: permut (const vector<int>& pos)
{
    // Get the representation of <pos>, which basically cares only about the numbers are different or the same
    // Example: if pos=3,6,6,8, rep=0,1,1,2
    //          if pos=1,2,2,2, rep=0,1,1,1
    //          if pos=1,3,7,9, rep=0,1,2,3
    vector<int> rep (pos.size(), 0);
    for(int i = 1; i < pos.size(); i++)
    {
        if (pos.at(i) != pos.at(i-1))
            rep.at(i) = rep.at(i-1) + 1;
    }

    // Would check keys in _ptable.
    // Optimize it if necessary.
    if (_ptable.count (rep) == 0)
        _ptable.insert ({rep, get_permutations (rep)});

    // Compute the permutations of the actual positions.
    // This would be less efficient than return the representations, and then convert each of them to positions in the loop outside.
    // However it is conceptually simpler, and may not be important in the whole algorithm.
    // Optimize if found necessary.
    auto [permut_poss, swaps] = _ptable.at (rep);
    for(auto& ppos : permut_poss)
    {
        for(int& pi : ppos)
        {
            pi = pos.at(pi);
        }
    }

    return {permut_poss, swaps};
}
#endif
