#ifndef BINNING_H
#define BINNING_H


#include <vector>
#include <algorithm>

#include "utils.h"

namespace npr {
namespace binning {

	class bin{
		public:
			bin(std::vector<xyzall_surfel_t> const& input_surfels, float upper_bound, float lower_bound, float axis_location) : 
																																pos_along_slicing_axis_(axis_location), 
																																lower_bound_size_ (lower_bound),
																																upper_bound_size_(upper_bound){
				generate_single_bin(input_surfels);
			}

			void evaluate_content_to_binary(bounding_rect const& bounding_corners, uint const resolution){
				uint const bool_vec_size = resolution * resolution; 
				std::vector<uint32_t> temp_binary_image (bool_vec_size, 0);
				auto cell_length = (bounding_corners.max_z - bounding_corners.min_z) / resolution;
	        	auto cell_width = (bounding_corners.max_x - bounding_corners.min_x) / resolution;
	        	uint const resolution_minus_one = resolution - 1;
		        for (auto const& current_surfel : content_){ 
		          auto x_index = std::min(resolution_minus_one, uint(std::max(0, int( (current_surfel.pos_coordinates[0] - bounding_corners.min_x) / cell_width)) ) );
		          auto z_index = std::min(resolution_minus_one, uint(std::max(0, int( (current_surfel.pos_coordinates[2]- bounding_corners.min_z) / cell_length)) ) );
		          int64_t cell_index = z_index * resolution + x_index;
		          ++temp_binary_image[cell_index];// = true;
		        }

		        binary_image_ = temp_binary_image;
			}

			void shrink_to_half_upper_bound(){
				auto new_upper_bound = upper_bound_size_ / 2.0;
				upper_bound_size_ = new_upper_bound;
				auto remove_lambda = [&](xyzall_surfel_t const& surfel){
					return(surfel.pos_coordinates[1] > (pos_along_slicing_axis_ + upper_bound_size_) );
				};
				content_.erase(std::remove_if(content_.begin(), content_.end(), remove_lambda), content_.end());
			}

			void shrink_to_half_lower_bound(){
				auto new_lower_bound = lower_bound_size_ / 2.0;
				lower_bound_size_ = new_lower_bound;
				auto remove_lambda = [&](xyzall_surfel_t const& surfel){
					return(surfel.pos_coordinates[1] < (pos_along_slicing_axis_ - lower_bound_size_) );
				};
				content_.erase(std::remove_if(content_.begin(), content_.end(), remove_lambda), content_.end());
			}

			void clear_binary_image(){
				binary_image_.clear();
				binary_image_.shrink_to_fit();
			}

		public:
			float pos_along_slicing_axis_;
			float lower_bound_size_;
			float upper_bound_size_;
			std::vector<xyzall_surfel_t> content_;
			std::vector<uint32_t> binary_image_;
			uint32_t bin_depth = 0;


		private:
			void generate_single_bin(std::vector<xyzall_surfel_t> const& input_surfels){
				content_.reserve(input_surfels.size());
				auto copy_lambda = [&]( xyzall_surfel_t const& surfel){return (surfel.pos_coordinates[1] >= (pos_along_slicing_axis_ - lower_bound_size_) ) && (surfel.pos_coordinates[1] <= (pos_along_slicing_axis_ + upper_bound_size_) );};
		    	auto it = std::copy_if(input_surfels.begin(), input_surfels.end(), content_.begin(), copy_lambda);
		    	content_.resize(std::distance(content_.begin(), it));
			}
	};

	struct evaluation_job{
		uint top_bin_id_ ;
		uint bottom_bin_id_ ;
		evaluation_job (uint top_id, uint bottom_id) : top_bin_id_(top_id), bottom_bin_id_(bottom_id) {}
	};



	bool evaluate_similarity(bin const& bin_A, 
							 bin const& bin_B,
							 bool using_merge_binning);

	bool evaluate_proximity(bin const& bin_A, 
                        	bin const& bin_B,
                        	float max_distance);

	std::vector<bin> generate_all_bins(std::vector<xyzall_surfel_t> const& all_surfels, 
									   float const inital_bin_half_height, 
									   uint& max_num_loops,
									   float max_distance_between_two_neighbouring_bins = -1.0); //TODO remove defaut value after split binning is removed
} //namespace binning
} //namespace npr

#endif //BINNING_H