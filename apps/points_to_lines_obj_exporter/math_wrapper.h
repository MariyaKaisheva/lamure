#ifndef MATH_WRAPPER_H
#define MATH_WRAPPER_H

#include <vector>

namespace npr {
namespace math {

    std::vector<double> fit_polynomial(std::vector<std::pair<double, double>> const& coordinates, int const polynomial_degree);
} //namespace math
} //namespace npr


#endif //MATH_WRAPPER_H