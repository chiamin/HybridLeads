#ifndef HYDRIDBASIS_UTILS_H_
#define HYDRIDBASIS_UTILS_H_

#include "itensor/all.h"

/**
 * @brief Element-wise comparison of 2 given itensors.
 *
 * @param t1 The first tensor.
 * @param t2 The second tensor.
 * @param atol The absolute tolerance. Default to 1e-12.
 * @return bool
 * @throws `std::logic_error` when the order of two tensors mismatches.
 */
bool ALLCLOSE(itensor::ITensor t1, itensor::ITensor t2, double atol = 1e-12) {
  if (order(t1) != order(t2)) {
    throw std::logic_error("The order of two tensors mismatch.");
  }
  auto check_close_zero = [&atol](itensor::Real r) {
    if (abs(r) > atol) {
      throw std::logic_error("Two tensors are not close.");
    }
  };
  try {
    t2.replaceInds(inds(t2), inds(t1));
    auto diff_t = t1 - t2;
    diff_t.visit(check_close_zero);  // visit() only accepts lambda func
  } catch (std::logic_error& e) {
    return false;
  }
  return true;
}

#endif
