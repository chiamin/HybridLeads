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
  FixedPointTensor(
      itensor::MPO const& mpo, int uniform_site,
      itensor::Args const& args = itensor::Args::global()
  ) {
    mpo_ = mpo;
    uniform_site_ = uniform_site;
    mpo_checker();
    int time_steps = args.getInt("time_steps", 30);
    itensor::Real dt = args.getReal("dt", 1e-12);
    int max_bond_dim = args.getInt("max_bond_dim", 1);
    itensor::Real tdvp_tol = args.getReal("tdvp_tol", 1e-12);
    int tdvp_max_iter = args.getInt("tdvp_max_iter", 40);
    itensor::Real ortho_tol = args.getReal("ortho_tol", 1e-12);
    int ortho_max_iter = args.getInt("ortho_max_iter", 20);
    RandGen::SeedType seed = args.getInt("seed", 0);
    itdvp_routine(
        time_steps, dt, max_bond_dim, tdvp_tol, tdvp_max_iter, ortho_tol,
        ortho_max_iter, seed
    );
  }

  itensor::ITensor get(std::string side) {
    std::map<std::string, int> mapper = {{"Left", -1}, {"Right", 1}};
    return (mapper[side] < 0) ? left_fixpt_tensor_ : right_fixpt_tensor_;
  }

  itensor::Index get_mpo_virtual_idx(std::string side) {
    return itensor::commonIndex(get(side), mpo_(uniform_site_));
  }

  itensor::Index get_mps_virtual_idx(std::string side) {
    return itensor::uniqueInds(get(side), mpo_(uniform_site_))(1);
  }

  int uniform_site() { return uniform_site_; }

 protected:
  itensor::MPO mpo_;
  int uniform_site_;
  itensor::Index mpo_left_idx_, mpo_right_idx_, phys_idx_;
  itensor::ITensor left_fixpt_tensor_, right_fixpt_tensor_;
  itensor::Real en_, err_;

  void mpo_checker() {
    int neighbour_site;
    if (itensor::order(mpo_(uniform_site_ - 1)) ==
        itensor::order(mpo_(uniform_site_))) {
      neighbour_site = uniform_site_ - 1;
    } else if (itensor::order(mpo_(uniform_site_)) == itensor::order(mpo_(uniform_site_ + 1))) {
      neighbour_site = uniform_site_ + 1;
    } else {
      throw std::invalid_argument(
          "Cannot find any neighbouring sites with same MPO tensor order."
      );
    }
    if (!ALLCLOSE(mpo_(uniform_site_), mpo_(neighbour_site))) {
      throw std::invalid_argument(
          "The `uniform site` should be picked from the bulk, with at least "
          "one neighbouring MPO tensor being identical with MPO tensor on this "
          "site."
      );
    }
  }

  void get_indices() {
    mpo_left_idx_ = commonIndex(mpo_(uniform_site_ - 1), mpo_(uniform_site_));
    mpo_right_idx_ = commonIndex(mpo_(uniform_site_), mpo_(uniform_site_ + 1));
    phys_idx_ = findIndex(mpo_(uniform_site_), "Site,0");
  }

  void itdvp_routine(
      int time_steps, itensor::Real dt, int max_bond_dim, itensor::Real tdvp_tol,
      int tdvp_max_iter, itensor::Real ortho_tol, int ortho_max_iter,
      RandGen::SeedType seed
  ) {
    get_indices();
    auto impo = mpo_(uniform_site_);
    auto imps = ITensor();  // ill-defined tensor
    auto
        [imps_left, imps_right, imps_left_center, imps_center, imps_left_idty,
         imps_right_idty] =
            itdvp_initial(
                impo, phys_idx_, mpo_left_idx_, mpo_right_idx_, imps, max_bond_dim,
                ortho_tol, ortho_max_iter, seed
            );
    itensor::Args args = {"ErrGoal=", tdvp_tol, "MaxIter", tdvp_max_iter};
    for (int i = 1; i <= time_steps; i++) {
      std::tie(en_, err_, left_fixpt_tensor_, right_fixpt_tensor_) = itdvp(
          impo, imps_left, imps_right, imps_left_center, imps_center, imps_left_idty,
          imps_right_idty, dt, args
      );
      LOG(INFO
      ) << std::printf("In time step %o, energy, error = %.3e, %.3e\n", i, en_, err_);
      if (args.getReal("ErrGoal") > tdvp_tol) {
        args.add("ErrGoal=", err_ * 0.1);
      }
    }
  }
};

#endif
