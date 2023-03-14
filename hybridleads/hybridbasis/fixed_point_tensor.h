#ifndef HYDRIDBASIS_FIXEDPOINTTENSOR_H_
#define HYDRIDBASIS_FIXEDPOINTTENSOR_H_

#include <glog/logging.h>

#include "hybridbasis/utils.h"
#include "itdvp/iTDVP.h"
#include "itensor/all.h"

class FixedPointTensor {
 public:
  /**
   * @brief Construct a new Fixed-Point Tensor object.
   *
   * @param mpo The instance of class: `itensor::MPO`.
   * @param uniform_site The site to which there are least one identical MPO
   * tensor in adjacent.
   * @param args Arguments containing the iTDVP parameters with these keywords:
   * (1) time evolution: `"time_steps"`, `"dt"`, `"tdvp_tol"`, `"tdvp_max_iter"`
   * (2) TDVP preworks: `"max_bond_dim"`, `"ortho_tol"`, `"ortho_max_iter"`,
   * `"seed"` (3) logging: `"log_to_std"`, default true
   * @throws `std::invalid_argument` when `uniform_site` does not give any
   * identical MPO tensors as its neighbour.
   */
  FixedPointTensor(itensor::MPO const& mpo, int uniform_site,
                   itensor::Args const& args = itensor::Args::global()) {
    _mpo = mpo;
    _uniform_site = uniform_site;
    mpo_checker();
    int time_steps = args.getInt("time_steps", 30);
    itensor::Real dt = args.getReal("dt", 1e-12);
    int max_bond_dim = args.getInt("max_bond_dim", 1);
    itensor::Real tdvp_tol = args.getReal("tdvp_tol", 1e-12);
    int tdvp_max_iter = args.getInt("tdvp_max_iter", 40);
    itensor::Real ortho_tol = args.getReal("ortho_tol", 1e-12);
    int ortho_max_iter = args.getInt("ortho_max_iter", 20);
    RandGen::SeedType seed = args.getInt("seed", 0);
    FLAGS_logtostderr = args.getBool("log_to_std", true);
    google::InitGoogleLogging("fixed_point_tensor");
    itdvp_routine(time_steps, dt, max_bond_dim, tdvp_tol, tdvp_max_iter,
                  ortho_tol, ortho_max_iter, seed);
  }

  itensor::ITensor get(std::string side) {
    std::map<std::string, int> mapper = {{"Left", -1}, {"Right", 1}};
    return (mapper[side] < 0) ? left_fixpt_tensor : right_fixpt_tensor;
  }

 protected:
  itensor::MPO _mpo;
  int _uniform_site;
  itensor::Index left_link, right_link, phys_bond;
  itensor::ITensor left_fixpt_tensor, right_fixpt_tensor;
  itensor::Real en, err;

  void mpo_checker() {
    int neighbour_site;
    if (order(_mpo(_uniform_site - 1)) == order(_mpo(_uniform_site))) {
      neighbour_site = _uniform_site - 1;
    } else if (order(_mpo(_uniform_site)) == order(_mpo(_uniform_site + 1))) {
      neighbour_site = _uniform_site + 1;
    } else {
      throw std::invalid_argument(
          "Cannot find any neighbouring sites with same MPO tensor order.");
    }
    if (!ALLCLOSE(_mpo(_uniform_site), _mpo(neighbour_site))) {
      throw std::invalid_argument(
          "The `uniform site` should be picked from the bulk, with at least "
          "one neighbouring MPO tensor being identical with MPO tensor on this "
          "site.");
    }
  }

  void get_indices() {
    left_link = commonIndex(_mpo(_uniform_site - 1), _mpo(_uniform_site));
    right_link = commonIndex(_mpo(_uniform_site), _mpo(_uniform_site + 1));
    phys_bond = findIndex(_mpo(_uniform_site), "Site,0");
  }

  void itdvp_routine(int time_steps, itensor::Real dt, int max_bond_dim,
                     itensor::Real tdvp_tol, int tdvp_max_iter,
                     itensor::Real ortho_tol, int ortho_max_iter,
                     RandGen::SeedType seed) {
    get_indices();
    auto impo = _mpo(_uniform_site);
    auto imps = ITensor();  // ill-defined tensor
    auto [imps_left, imps_right, imps_left_center, imps_center, imps_left_idty,
          imps_right_idty] =
        itdvp_initial(impo, phys_bond, left_link, right_link, imps,
                      max_bond_dim, ortho_tol, ortho_max_iter, seed);
    itensor::Args args = {"ErrGoal=", tdvp_tol, "MaxIter", tdvp_max_iter};
    for (int i = 1; i <= time_steps; i++) {
      std::tie(en, err, left_fixpt_tensor, right_fixpt_tensor) =
          itdvp(impo, imps_left, imps_right, imps_left_center, imps_center,
                imps_left_idty, imps_right_idty, dt, args);
      LOG(INFO) << std::printf("In time step %o, energy, error = %.3e, %.3e\n",
                               i, en, err);
      if (args.getReal("ErrGoal") > tdvp_tol) {
        args.add("ErrGoal=", err * 0.1);
      }
    }
  }
};

#endif
