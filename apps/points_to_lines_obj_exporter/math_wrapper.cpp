#include "math_wrapper.h"
	
#include <alglib/cpp/include/interpolation.h>

namespace npr {
namespace math {

	std::vector<double> 
	fit_polynomial (std::vector<std::pair<double, double>> const& coordinates, int const polynomial_degree) {

		auto num_points = coordinates.size();
		std::vector<double> x_coordinates;
		std::vector<double> y_coordinates;

		for (auto const& current_point_coords : coordinates){
			x_coordinates.push_back(current_point_coords.first);
			y_coordinates.push_back(current_point_coords.second);
		}
		//assign point coord-component to alglib array
		alglib::real_1d_array    point_x_coords;
		point_x_coords.setcontent(num_points, (double*) &x_coordinates[0]);
		alglib::real_1d_array    point_y_coords;
		point_y_coords.setcontent(num_points, (double*) &y_coordinates[0]);

    	std::vector<double> weights(num_points, 1.0); //usually all=1
    	alglib::real_1d_array    point_weights;
    	point_weights.setcontent(num_points, (double*) &weights[0]); 

    	std::vector<double> x_coords_contrains_vec;
    	x_coords_contrains_vec.push_back(point_x_coords[0]);
    	x_coords_contrains_vec.push_back(point_x_coords[num_points-1]); 
    	alglib::real_1d_array    point_x_coord_constraints;
    	point_x_coord_constraints.setcontent(2, (double*) &x_coords_contrains_vec[0]);

    	std::vector<double> y_coords_contrains_vec;
    	y_coords_contrains_vec.push_back(point_y_coords[0]);
    	y_coords_contrains_vec.push_back(point_y_coords[num_points-1]); 
    	alglib::real_1d_array    point_y_coord_constraints;
    	point_y_coord_constraints.setcontent(2, (double*) &y_coords_contrains_vec[0]);


    	alglib::integer_1d_array point_constraint_types; //0 for C^0 Constraints; 1 for C^1 Constraints.
    	std::vector<alglib::ae_int_t> contraint_types(2, 0);
    	point_constraint_types.setcontent(2, (alglib::ae_int_t*) &contraint_types[0]);

    	alglib::ae_int_t         number_of_poly_degree_plus_1 = polynomial_degree + 1; //--> 1: linear func, 2: quadratic, etc.
    	alglib::ae_int_t         error_info;
    	alglib::barycentricinterpolant bary_interp; //needs to be converted with the alglib::PolynomialBar2Pow() function to be in standard form
    	alglib::polynomialfitreport fit_report; //contains: RMSError, AvgError, AvgRelError and MaxError

    	alglib::polynomialfitwc(point_x_coords, point_y_coords, 
                            point_weights, point_x_coord_constraints, point_y_coord_constraints, point_constraint_types,
                            number_of_poly_degree_plus_1, error_info, bary_interp, fit_report);

    	//cast barycentric interpolant to coefficients for x^0 ... x^(n-1))
    	alglib::real_1d_array fit_coefficients;
    	alglib::polynomialbar2pow(bary_interp, fit_coefficients);
    	std::vector<double> coefficients (number_of_poly_degree_plus_1); 
    	for (int coef_idx = 0; coef_idx < number_of_poly_degree_plus_1; ++coef_idx){
    		coefficients.at(coef_idx) = fit_coefficients[coef_idx];
    	}
    	return coefficients;
	};

    /*
    alglib::polynomialfitwc(point_x_coords, point_y_coords, 
                            point_weights, point_x_coord_constraints, point_y_coord_constraints, point_constraint_types,
                            number_of_poly_degree_plus_1, error_info, bary_interp, fit_report);

    //cast barycentric interpolant to coefficients for x^0 ... x^(n-1))
    alglib::real_1d_array fit_coefficients;
    alglib::polynomialbar2pow(bary_interp, fit_coefficients);*/
} //namespace math
} //namespace npr