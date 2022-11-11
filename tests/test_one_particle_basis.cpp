#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "itensor/all.h"

#include "../kbasis/OneParticleBasis.h"

using namespace itensor;
using namespace Catch;


TEST_CASE ( "Check matrix elements of tight-binding Hamiltonian", "[tight_binding_Hamilt]" ){
    int mat_dim = 4;
    auto mat = tight_binding_Hamilt(mat_dim, 0.5, 1.);
    float expected_mat[4][4] = {
        {-1.0, -0.5, 0.0, 0.0},
        {-0.5, -1.0, -0.5, 0.0},
        {0.0, -0.5, -1.0, -0.5},
        {0.0, 0.0, -0.5, -1.0}
    };

    for (int i = 0; i < mat_dim; i++){
        for (int j = 0; j < mat_dim; j++){
            CHECK( mat(i, j) == Approx(expected_mat[i][j]).epsilon(1e-12) );
        }
    }
}


// TEST_CASE ( "Check one particle basis", ){

// }
