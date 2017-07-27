#ifndef BINNING_H
#define BINNING_H


#include <vector>

#include "utils.h"

	struct plane{
		float height_;
		std::vector<std::vector<bool>> binary_image_;
	};


	using bins_t = std::vector<xyzall_surfel_t>;
	namespace binning {

		std::vector<bins_t> generate_bins(std::vector<xyzall_surfel_t> const& all_surfels, float inital_bin_thickness);
		bool evaluate_similarty(plane const& plane_A, plane const& plane_B);
		plane convert_bin_to_plane(bins_t const& bin_of_surfels, uint const resolution, bounding_rect const& bounding_corners);

	}		

#endif //BINNING_H