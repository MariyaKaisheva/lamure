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
				float absolute_rotation_offset_angle,
				float relative_rotation_offset_angle) :  pos_along_slicing_axis_(plane_origin.y),
														 lower_bound_size_(range_threshols),
														 upper_bound_size_(range_threshols) {

				generate_single_radial_bin(input_surfels, 
										   plane_origin,
										   bounding_sphere_radius,
										   absolute_rotation_offset_angle,
										   relative_rotation_offset_angle);
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

		public:
			float pos_along_slicing_axis_;
			float lower_bound_size_;
			float upper_bound_size_;
			std::vector<xyzall_surfel_t> content_;
			std::vector<uint32_t> binary_image_;
			uint32_t bin_depth = 0;
			std::pair<scm::math::vec3f, scm::math::vec3f> bounding_corners_;
			scm::math::mat4f radial_rotation_mat_;


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
										    float absolute_rotation_offset_angle,
										    float relative_rotation_offset_angle) {


				scm::math::vec3f rotation_orientation(0.0f, 0.0f, 1.0f);
				scm::math::mat4f rot_mat = scm::math::make_rotation(relative_rotation_offset_angle, rotation_orientation);

				scm::math::mat4f origin_transl_mat = scm::math::make_translation(-plane_origin.x, 
																				 -plane_origin.y,
																				 -plane_origin.z);

				scm::math::mat4f final_plane_rotation_mat = scm::math::inverse(origin_transl_mat) * rot_mat * origin_transl_mat;
				radial_rotation_mat_ =  final_plane_rotation_mat; //save transformations to recover the original model after feature extraction is done

				scm::math::mat4f point_to_non_AABB_space_mat = scm::math::inverse(origin_transl_mat) * scm::math::inverse(rot_mat) * origin_transl_mat;


				content_.reserve(input_surfels.size());


				std::vector<xyzall_surfel_t> points_in_rotated_plane_space(input_surfels);
				
				// TRANSFORM points_in_rotated_plane_space BY point_to_non_AABB_space_mat
				utils::transform_surfels_by_matrix(points_in_rotated_plane_space, point_to_non_AABB_space_mat);


				#if 1//work with PRISM-shaped slicing volume
				//min and max corners of the 'boxy' bounding search region
				scm::math::vec3f min_box_corner(plane_origin.x, 
											plane_origin.y - bounding_sphere_radius * 0.01 , 
											plane_origin.z - bounding_sphere_radius);

				scm::math::vec3f max_box_corner(plane_origin.x + bounding_sphere_radius,
											plane_origin.y + bounding_sphere_radius * 0.01 ,
											plane_origin.z + bounding_sphere_radius);

				//obsolite?
				//bounding_corners_ = std::make_pair(min_box_corner, max_box_corner); //store 'flat' values for binary image computation

				auto copy_lambda = [&]( xyzall_surfel_t const& surfel) {
					bool check_x_coord = (surfel.pos_coordinates[0] >= min_box_corner.x && surfel.pos_coordinates[0] <= max_box_corner.x);
					bool check_y_coord = (surfel.pos_coordinates[1] >= min_box_corner.y && surfel.pos_coordinates[1] <= max_box_corner.y);
					bool check_z_coord = (surfel.pos_coordinates[2] >= min_box_corner.z && surfel.pos_coordinates[2] <= max_box_corner.z);
					return (check_x_coord && check_y_coord && check_z_coord);
				};

				#else //work with PYRAMID-shaped slicing volume
				//point-in-pyramid test
				/**pyramid volume is implicitly defined by 4 tringle sides;
			     **these are further specifed via 3 times 90-deg rotation of a single explicitly defined tringle

			      						*B

										
			     	*A (plane origin)


			     						*C
			     */
				scm::math::vec4f vertex_A_hom_coord(plane_origin.x, 
													plane_origin.y,
													plane_origin.z,
													1.0);

				float degree_as_radians = absolute_rotation_offset_angle * 3.14159265359 / 180.0;

				float pyramid_half_base_side = bounding_sphere_radius * tan(degree_as_radians / 2.0);



				scm::math::vec4f  vertex_B_hom_coord(plane_origin.x + bounding_sphere_radius, 
										   			 plane_origin.y + pyramid_half_base_side,
										   			 plane_origin.z + pyramid_half_base_side,
										   			 1.0);

				scm::math::vec4f  vertex_C_hom_coord(plane_origin.x + bounding_sphere_radius, 
						   				   			 plane_origin.y + pyramid_half_base_side,
						   				   			 plane_origin.z - bounding_sphere_radius,
						   				   			 1.0);

				scm::math::vec4f AB_vec = vertex_B_hom_coord - vertex_A_hom_coord;
				scm::math::vec4f AC_vec = vertex_C_hom_coord - vertex_A_hom_coord;


				std::vector<scm::math::vec4f> pyramid_side_normals_vec; 
				//TODO change variable names 
				scm::math::mat4f side_flip_rot_mat = scm::math::make_rotation(180.0f, 1.0f, 0.0f, 0.0f);
				scm::math::vec4f pyramid_side_normal = scm::math::cross(scm::math::vec3f(AB_vec), scm::math::vec3f(AC_vec) );
				pyramid_side_normals_vec.push_back(pyramid_side_normal);
				pyramid_side_normals_vec.push_back(side_flip_rot_mat * pyramid_side_normal);
				pyramid_side_normals_vec.push_back(scm::math::vec4f(0.0, 0.0, 1.0, 0.0));
				pyramid_side_normals_vec.push_back(scm::math::vec4f(0.0, 0.0, -1.0, 0.0));


				auto copy_lambda = [&]( xyzall_surfel_t const& surfel) {
					scm::math::vec4f origin_to_surfel_vec = scm::math::vec4f(surfel.pos_coordinates[0] - plane_origin.x, 
		    																 surfel.pos_coordinates[1] - plane_origin.y,
		    																 surfel.pos_coordinates[2] - plane_origin.z,
		    																 0.0);

					for (auto const& side_normal : pyramid_side_normals_vec){
						float dot_product_side_normal_point_vec = scm::math::dot(scm::math::normalize(scm::math::vec3f(side_normal) ), scm::math::normalize(scm::math::vec3f(origin_to_surfel_vec) ));
						if (dot_product_side_normal_point_vec > 0.0){
							return false;
						}
					}
					
					//bool is_distance_less_than_pyramid_height = surfel.pos_coordinates[0] >= plane_origin.x;

					//return is_distance_less_than_pyramid_height;

					return true;
				};
				#endif //switch between 2 types/shpes of  slicing volume 

		    	auto it = std::copy_if(points_in_rotated_plane_space.begin(), points_in_rotated_plane_space.end(), content_.begin(), copy_lambda);
		    	content_.resize(std::distance(content_.begin(), it));
		    	std::cout << "BIN CONTENT: " << content_.size() << "\n";

		    	//project bin content
		    	for(auto & surfel : content_){

		    		//drop y-coordinate to rough projection onto the slicing plane
		    		scm::math::vec3f straight_projected_surfel_pos = scm::math::vec3f(surfel.pos_coordinates[0], plane_origin.y, surfel.pos_coordinates[2]);
		    		scm::math::vec3f sphere_origin_to_straight_projected_surfel_vec =  scm::math::normalize(straight_projected_surfel_pos - plane_origin); //plane_origin = center of rotation

		    		//scale-corection to get correct radial projection
		    		scm::math::vec3f sphere_origin_to_unprojected_surfel_vec = scm::math::vec3f(surfel.pos_coordinates[0] - plane_origin.x, 
		    																					surfel.pos_coordinates[1] - plane_origin.y,
		    																					surfel.pos_coordinates[2] - plane_origin.z);

		    		if(0.0f == surfel.pos_coordinates[0] && 0.0f == surfel.pos_coordinates[1] && 0.0f == surfel.pos_coordinates[2] ) {
		    			continue;
		    		}

		    		float scaling_factor = scm::math::length(sphere_origin_to_unprojected_surfel_vec);
		    		scm::math::vec3f radial_projected_surfel = scaling_factor * sphere_origin_to_straight_projected_surfel_vec + plane_origin;

		    		surfel.pos_coordinates[0] = radial_projected_surfel.x;
		    		surfel.pos_coordinates[1] = radial_projected_surfel.y;
		    		surfel.pos_coordinates[2] = radial_projected_surfel.z;
		    	}
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
									   bool radial_slicing,
									   scm::math::vec3f & bounding_sphere_center,
									   float max_distance_between_two_neighbouring_bins = -1.0,
									   bool verbose = false); //TODO remove defaut value after split binning is removed
} //namespace binning
} //namespace npr

#endif //BINNING_H