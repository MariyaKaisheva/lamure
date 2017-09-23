#ifndef BINNING_H
#define BINNING_H


#include <vector>
#include <algorithm>

#include <scm/gl_core/math.h>

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

			bin(std::vector<xyzall_surfel_t> const& input_surfels,
				float range_threshols, 
				scm::math::vec3f plane_origin,
				float bounding_sphere_radius,
				float rotation_offset_angle) :  pos_along_slicing_axis_(plane_origin.y),
												lower_bound_size_(range_threshols),
												upper_bound_size_(range_threshols) {

				generate_single_radial_bin(input_surfels, 
										   plane_origin,
										   bounding_sphere_radius,
										   rotation_offset_angle);
			}

			void evaluate_content_to_binary(bounding_rect const& bounding_corners,
											uint const resolution,
											scm::math::vec3f const& bounding_sphere_center,
										    bool radial_slicing){

				uint const bool_vec_size = resolution * resolution; 
				std::vector<uint32_t> temp_binary_image (bool_vec_size, 0);
				uint const resolution_minus_one = resolution - 1; 

				if(!radial_slicing){
					auto cell_length = (bounding_corners.max_z - bounding_corners.min_z) / resolution;
		        	auto cell_width = (bounding_corners.max_x - bounding_corners.min_x) / resolution;
		        		
	        		for (auto const& current_surfel : content_){
			          auto x_index = std::min(resolution_minus_one, uint(std::max(0, int( (current_surfel.pos_coordinates[0] - bounding_corners.min_x) / cell_width)) ) );
			          auto z_index = std::min(resolution_minus_one, uint(std::max(0, int( (current_surfel.pos_coordinates[2]- bounding_corners.min_z) / cell_length)) ) );
			          int64_t cell_index = z_index * resolution + x_index;
			          ++temp_binary_image[cell_index];// = true;
			        }
	        	} else{
	        		float length_of_bounded_slicing_plane = utils::compute_distance(bounding_corners_.first, bounding_sphere_center);
	        		auto cell_length = (length_of_bounded_slicing_plane) / resolution;
	        		float width_of_bounded_slicing_plane =  utils::compute_distance(bounding_corners_.second, bounding_sphere_center);
	        		auto cell_width = (width_of_bounded_slicing_plane) / resolution;
	        		for (auto const& current_surfel : content_){
	        		  //std::cout << "current_surfel_pos.x " << current_surfel.pos_coordinates[0] << " bounding_corners_MIN.x " << bounding_corners_.first.x  << " bounding_corners_MAX.x " << bounding_corners_.second.x <<  "\n";
	        		  //std::cout << "cell_width " << cell_width<<  "\n";
	        		  auto x_index = std::min(resolution_minus_one, uint(std::max(0, int( (current_surfel.pos_coordinates[0] - bounding_corners_.first.x) / cell_width)) ) );
			          auto z_index = std::min(resolution_minus_one, uint(std::max(0, int( (current_surfel.pos_coordinates[2]- bounding_corners_.first.z) / cell_length)) ) );
			          //std::cout << "x_index " << x_index << " z_index " << z_index <<  "\n";
			          int64_t cell_index = z_index * resolution + x_index;
			          ++temp_binary_image[cell_index];// = true;
	        		}
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

			/*bounding_rect get_bounding_corners(){
				auto min_corner = bounding_corners_.first;
				auto max_corner = bounding_corners_.second;
				return utils::bounding_rect min_max_corner_coordinates(min_corner.x, min_corner.y, min_corner.z, max_corner.x, max_corner.y, max_corner.z);
			}*/

		public:
			float pos_along_slicing_axis_;
			float lower_bound_size_;
			float upper_bound_size_;
			std::vector<xyzall_surfel_t> content_;
			std::vector<uint32_t> binary_image_;
			uint32_t bin_depth = 0;
			std::pair<scm::math::vec3f, scm::math::vec3f> bounding_corners_;


		private:
			void generate_single_bin(std::vector<xyzall_surfel_t> const& input_surfels){
				content_.reserve(input_surfels.size());
				auto copy_lambda = [&]( xyzall_surfel_t const& surfel){return (surfel.pos_coordinates[1] >= (pos_along_slicing_axis_ - lower_bound_size_) ) && (surfel.pos_coordinates[1] <= (pos_along_slicing_axis_ + upper_bound_size_) );};
		    	auto it = std::copy_if(input_surfels.begin(), input_surfels.end(), content_.begin(), copy_lambda);
		    	content_.resize(std::distance(content_.begin(), it));
			}

			void generate_single_radial_bin(std::vector<xyzall_surfel_t> const& input_surfels, 
										    scm::math::vec3f const& plane_origin,
										    float bounding_sphere_radius,
										    float rotation_offset_angle) {


				scm::math::vec3f rotation_orientation(0.0f, 0.0f, 1.0f);
				scm::math::mat4f rot_mat = scm::math::make_rotation(rotation_offset_angle, rotation_orientation);
				scm::math::vec4f hom_rot_plane_normal = rot_mat * scm::math::vec4(0.0f, 1.0f, 0.0f, 0.0f);
				scm::math::vec3f plane_normal     = scm::math::vec3f(hom_rot_plane_normal.x, hom_rot_plane_normal.y, hom_rot_plane_normal.z);


				//min and max corners of the 'flat' bounding search region
				scm::math::vec3f min_corner(plane_origin.x, 
											plane_origin.y, 
											plane_origin.z - bounding_sphere_radius);

				scm::math::vec3f max_corner(plane_origin.x + bounding_sphere_radius,
											plane_origin.y,
											plane_origin.z + bounding_sphere_radius);
				bounding_corners_ = std::make_pair(min_corner, max_corner); //store 'flat' values for binary image computation

				//rotate max corner to match the actual orientation of the slicing place associated wtih this bin
				//but first translate corner to the origin of the local coordinate system for correct rotaion results
				scm::math::mat4f origin_transl_mat = scm::math::make_translation(-plane_origin.x, 
																				 -plane_origin.y,
																				 -plane_origin.z);
				scm::math::vec4f corner_at_origin = origin_transl_mat * scm::math::vec4f(max_corner, 1.0f);
				scm::math::vec4f rotated_corner = rot_mat * corner_at_origin; //actual rotation
				scm::math::vec4f rotated_corner_at_initial_location = scm::math::inverse(origin_transl_mat) * rotated_corner;
				max_corner = scm::math::vec3f(rotated_corner_at_initial_location.x, rotated_corner_at_initial_location.y, rotated_corner_at_initial_location.z);


				//expand the 'flat' search region to a search box

				std::cout << "PLANE NORMAL: " << plane_normal << "\n";
				std::cout << "bounding sphere_radus: " << bounding_sphere_radius << "\n"; 
				scm::math::vec3f trans_vec (bounding_sphere_radius * plane_normal * 1.0);
				scm::math::mat4f translation_mat = scm::math::make_translation(trans_vec);
				scm::math::vec4f translated_corner = translation_mat * scm::math::vec4f(max_corner, 1.0f);
				max_corner = scm::math::vec3f(translated_corner.x, translated_corner.y, translated_corner.z); 

				
				scm::math::vec3f trans_vec_for_min_corner (bounding_sphere_radius * ( plane_normal * -1 *  1.0));
				scm::math::mat4f translation_mat_for_min_corner = scm::math::make_translation(trans_vec_for_min_corner);
				scm::math::vec4f translated_min_corner = translation_mat_for_min_corner * scm::math::vec4f(min_corner, 1.0f);
				min_corner = scm::math::vec3f(translated_min_corner.x, translated_min_corner.y, translated_min_corner.z); 

				std::cout << "MIN corner (" << min_corner.x << ", " << min_corner.y << ", "  << min_corner.z << ")\n ";
				std::cout << "MAX corner (" << max_corner.x << ", " << max_corner.y << ", "  << max_corner.z << ")\n ";
				std::cout << "BS RADIUS: " << bounding_sphere_radius << "\n";
				content_.reserve(input_surfels.size());


				auto copy_lambda = [&]( xyzall_surfel_t const& surfel) {
					
										bool check_x_coord = (surfel.pos_coordinates[0] >= min_corner.x && surfel.pos_coordinates[0] <= max_corner.x);
										bool check_y_coord = (surfel.pos_coordinates[1] >= min_corner.y && surfel.pos_coordinates[1] <= max_corner.y);
										bool check_z_coord = (surfel.pos_coordinates[2] >= min_corner.z && surfel.pos_coordinates[2] <= max_corner.z);
										return (check_x_coord && check_y_coord && check_z_coord);

									};

				std::cout << "BEFORE COPYING: " << input_surfels.size() << "\n";
		    	auto it = std::copy_if(input_surfels.begin(), input_surfels.end(), content_.begin(), copy_lambda);
		    	content_.resize(std::distance(content_.begin(), it));
				std::cout << "AFTER COPYING: " << content_.size() << "\n";

				/*
		    	//project content
		    	for(auto & surfel : content_){
		    		scm::math::vec3f surfel_original_pos = scm::math::vec3f(surfel.pos_coordinates[0], surfel.pos_coordinates[1], surfel.pos_coordinates[2]);
		    		//vector between boundig sphere center and current surfel 
		    		scm::math::vec3f origin_surfel_vec = scm::math::vec3f(surfel_original_pos.x - plane_origin.x,
		    															  surfel_original_pos.y - plane_origin.y,
		    															  surfel_original_pos.z - plane_origin.z);
		    		//compute projection angel 
		    		scm::math::vec3f normalized_origin_surfel_vec = utils::normalize(origin_surfel_vec);
		    		scm::math::vec3f slicing_plane_vec = scm::math::vec3f(min_corner.x - plane_origin.x,
		    															  min_corner.y - plane_origin.y,
		    															  min_corner.z - plane_origin.z);
		    		scm::math::vec3f normalized_plane_vec = utils::normalize(slicing_plane_vec);
		    		float projection_angle =  acos( utils::dot(normalized_origin_surfel_vec, normalized_plane_vec));

		    		//surfel projection
		    		scm::math::mat4f projection_mat = //scm::math::inverse
		    			(scm::math::make_rotation(projection_angle, rotation_orientation));
		    		scm::math::vec4f rotated_origin_surfel_vec = projection_mat * origin_surfel_vec;

		    		//TODO check if this is correct!!!
		    		surfel.pos_coordinates[0] = rotated_origin_surfel_vec.x;
		    		surfel.pos_coordinates[1] = rotated_origin_surfel_vec.y;
		    		surfel.pos_coordinates[2] = rotated_origin_surfel_vec.z;
		    	}
		    	*/	
			}
	};

	struct evaluation_job{
		uint top_bin_id_ ;
		uint bottom_bin_id_ ;
		evaluation_job (uint top_id, uint bottom_id) : top_bin_id_(top_id), bottom_bin_id_(bottom_id) {}
	};



	bool evaluate_similarity(bin const& bin_A, 
							 bin const& bin_B);

	bool evaluate_proximity(bin const& bin_A, 
                        	bin const& bin_B,
                        	float max_distance);

	std::vector<bin> generate_all_bins(std::vector<xyzall_surfel_t> const& all_surfels, 
									   float const inital_bin_half_height, 
									   uint& max_num_loops,
									   float max_distance_between_two_neighbouring_bins = -1.0,
									   bool verbose = false); //TODO remove defaut value after split binning is removed
} //namespace binning
} //namespace npr

#endif //BINNING_H