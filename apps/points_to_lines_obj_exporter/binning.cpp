#include "binning.h"

#include "utils.h"

bool binning::
 evaluate_similarty(plane const& plane_A, plane const& plane_B){
	return false;
}

std::vector<bins_t> binning::
 generate_bins(std::vector<xyzall_surfel_t> const& all_surfels, float inital_bin_thickness){
 	std::vector<bins_t> layers;
 	auto bounding_corners = utils::compute_bounding_corners(all_surfels);

 	return layers;
}

plane binning::
 convert_bin_to_plane(bins_t const& bin_of_surfels, uint resolution, bounding_rect const& bounding_corners){
 	plane plane_from_bin;


 	return plane_from_bin;
}