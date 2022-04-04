#include "array_bridge.h"
#include <aero/array/array.hpp>

namespace aero {

  Array::Array(std::size_t number_of_elements) :
    values_(number_of_elements, 0.0) {}

  Array:: Array(std::size_t number_of_elements,
      Real initial_value) : values_(number_of_elements, initial_value) {}

  Array::Array(const std::vector<Real> &values) : values_(values) {}

  Array& Array::operator=(const std::vector<Real> &values) {
    this->values_ = values;
    return *this;
  }

  Array* Array::clone() const {
    return new Array(*this);
  }

  void Array::copy_in(const Real *input) {
    for (int i=0; i<this->values_.size(); ++i) this->values_[i] = input[i];
  }

  void Array::copy_in(const std::vector<Real> &input) {
    this->values_ = input;
  }

  void Array::copy_out(Real *output) const {
    for (int i=0; i<this->values_.size(); ++i) output[i] = this->values_[i];
  }

  void Array::copy_out(std::vector<Real> &output) const {
    output = this->values_;
  }

  std::size_t Array::size() const {
    return this->values_.size();
  }

} // namespace aero
