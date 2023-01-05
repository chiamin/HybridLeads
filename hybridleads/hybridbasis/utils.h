#include "itensor/all.h"

using namespace itensor;

/**
 * @brief Element-wise comparison of 2 given itensors.
 *
 * @param t1 The first tensor.
 * @param t2 The second tensor.
 * @param atol The absolute tolerance. Default to 1e-12.
 * @return bool
 */
bool ALLCLOSE(ITensor t1, ITensor t2, double atol = 1e-12) {
  REQUIRE(order(t1) == order(t2));
  auto check_close_zero = [&atol](Real r) {
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
