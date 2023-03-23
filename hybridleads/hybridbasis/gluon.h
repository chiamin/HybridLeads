#ifndef HYDRIDBASIS_GLUON_H_
#define HYDRIDBASIS_GLUON_H_

#include <armadillo>

#include "hybridbasis/fixed_point_tensor.h"
#include "itensor/all.h"

class Gluon {
 public:
  Gluon(
      itensor::MPO const& mpo, itensor::SiteSet const& sites, int left_lead_size,
      int right_lead_size, itensor::Args const& itdvp_args = itensor::Args::global()
  )
      : left_fxpts_(mpo, 2, itdvp_args),
        right_fxpts_(mpo, itensor::length(mpo) - 1, itdvp_args) {
    mpo_ = mpo;
    sites_ = sites;
    left_lead_size_ = left_lead_size;
    right_lead_size_ = right_lead_size;
    sys_size_ = itensor::length(sites_);
    if (left_lead_size_ + sys_size_ + right_lead_size_ != itensor::length(mpo_)) {
      throw std::invalid_argument("Total size does not match with the given MPO.");
    }
  }

  int sys_size() { return sys_size_; }

  itensor::MPO sys_mpo() {
    itensor::MPO sys_mpo(sites_);
    for (int site = 1; site <= sys_size_; ++site) {
      itensor::IndexSet new_site_inds = itensor::siteInds(sys_mpo, site);
      sys_mpo.ref(site) = mpo_(left_lead_size_ + site);
      sys_mpo.ref(site).replaceInds(
          itensor::siteInds(mpo_, left_lead_size_ + site), new_site_inds
      );
    }
    sys_mpo.ref(1).replaceInds(
        {itensor::leftLinkIndex(mpo_, left_lead_size_ + 1)},
        {left_fxpts_.get_mpo_virtual_idx("Left")}
    );
    sys_mpo.ref(sys_size_).replaceInds(
        {itensor::rightLinkIndex(mpo_, left_lead_size_ + sys_size_)},
        {right_fxpts_.get_mpo_virtual_idx("Right")}
    );
    return sys_mpo;
  }

  itensor::MPS random_init_state() {
    itensor::InitState state(sites_);
    for (int i : range1(sys_size_)) {
      if (i % 2 == 1)
        state.set(i, "1");
      else
        state.set(i, "0");
    }
    itensor::MPS mps = itensor::randomMPS(state);
    itensor::ITensor first_site_rand_ts = itensor::randomITensor(itensor::IndexSet(
        {left_fxpts_.get_mps_virtual_idx("Left")}, itensor::inds(mps(1))
    ));
    itensor::ITensor last_site_rand_ts = itensor::randomITensor(itensor::IndexSet(
        itensor::inds(mps(sys_size_)), {right_fxpts_.get_mps_virtual_idx("Right")}
    ));
    mps.set(1, first_site_rand_ts);
    mps.set(sys_size_, last_site_rand_ts);
    return mps;
  }

  itensor::ITensor left_env() { return left_fxpts_.get("Left"); }

  itensor::ITensor right_env() { return right_fxpts_.get("Right"); }

 protected:
  itensor::MPO mpo_;
  itensor::SiteSet sites_;
  int left_lead_size_, right_lead_size_, sys_size_;
  FixedPointTensor left_fxpts_, right_fxpts_;
};

#endif
