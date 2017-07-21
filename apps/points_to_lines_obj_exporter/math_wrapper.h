#ifndef MATH_WRAPPER_H
#define MATH_WRAPPER_H

#include <vector>

	namespace math {

	    std::vector<double> fit_polynomial(std::vector<std::pair<double, double>> const& coordinates, int const polynomial_degree);
    }


#endif //MATH_WRAPPER_H