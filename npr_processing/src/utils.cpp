#include <lamure/npr/utils.h>

#include <lamure/ren/bvh.h>

#include <memory>

namespace npr {
namespace utils {
	float compute_avg_min_distance(std::vector<xyzall_surfel_t> const& input_data, uint32_t const num_cells_pro_dim_x, uint32_t const num_cells_pro_dim_y, uint32_t const num_cells_pro_dim_z) {
		float average_min_distance = 0.0f;
		float search_radius = 0.0f;
		float inital_search_radius = 0.0f;

		std::vector<grid_cell> cells_vec;
		generate_cells(input_data, inital_search_radius, cells_vec, num_cells_pro_dim_x, num_cells_pro_dim_y, num_cells_pro_dim_z);
		
		#pragma omp parallel for reduction (+ : average_min_distance)
		  for(uint point_counter = 0; point_counter < input_data.size(); ++point_counter){
		  	auto const& current_point = input_data[point_counter];
			search_radius = inital_search_radius; 
			auto start_coord = lamure::vec3f(current_point.pos_coordinates[0], current_point.pos_coordinates[1], current_point.pos_coordinates[2]);
			float min_distance = std::numeric_limits<float>::max();

			std::vector<grid_cell*> candidate_cells;
			find_candidate_cells(cells_vec, current_point.pos_coordinates, search_radius, candidate_cells);
			bool not_enough_candidates = test_for_sufficency(candidate_cells);

			while(not_enough_candidates){
				search_radius +=  1.5*inital_search_radius;
				candidate_cells.clear();
				find_candidate_cells(cells_vec, current_point.pos_coordinates, search_radius, candidate_cells);
				not_enough_candidates = test_for_sufficency(candidate_cells);

				//if degenerated: break
				if(search_radius > 1000000.0f) {
					break;
				}
			}

			for (auto const& curent_cell : candidate_cells) {
				for(auto const& cell_point : curent_cell->content_){

					auto end_coord = lamure::vec3f(cell_point.pos_coordinates[0], cell_point.pos_coordinates[1], cell_point.pos_coordinates[2]);
					if(!(end_coord.x == start_coord.x && end_coord.y == start_coord.y && end_coord.z == start_coord.z)){
						auto current_distance = compute_distance(start_coord, end_coord);
						min_distance = std::min(min_distance, current_distance);
					}
				}

			}

			if( min_distance == std::numeric_limits<float>::max() ) {
				min_distance = 0.0f;
			}

			average_min_distance += min_distance;
		}

		return average_min_distance / input_data.size();
	}


	void transform(std::vector<line>& line_data_vec, scm::math::mat4f const& transformation_mat) {

	    size_t num_lines_in_vector = line_data_vec.size();
	    #pragma omp parallel for
	    for(size_t line_idx = 0; line_idx < num_lines_in_vector; ++line_idx) {
			auto& line = line_data_vec[line_idx];
			scm::math::vec3f inital_start_coord = scm::math::vec4f(line.start.pos_coordinates_[0], line.start.pos_coordinates_[1], line.start.pos_coordinates_[2], 1.0f);
			float transformed_start_x = (transformation_mat * inital_start_coord).x;
			float transformed_start_y = (transformation_mat * inital_start_coord).y;
			float transformed_start_z = (transformation_mat * inital_start_coord).z;
			scm::math::vec3f inital_end_coord = scm::math::vec4f(line.end.pos_coordinates_[0], line.end.pos_coordinates_[1], line.end.pos_coordinates_[2], 1.0f);
			float transformed_end_x = (transformation_mat * inital_end_coord).x;
			float transformed_end_y = (transformation_mat * inital_end_coord).y;
			float transformed_end_z = (transformation_mat * inital_end_coord).z;
			line.start.pos_coordinates_[0] = transformed_start_x;
			line.start.pos_coordinates_[1] = transformed_start_y;
			line.start.pos_coordinates_[2] = transformed_start_z;
			line.end.pos_coordinates_[0] = transformed_end_x;
			line.end.pos_coordinates_[1] = transformed_end_y;
			line.end.pos_coordinates_[2] = transformed_end_z;
		}
	}

	std::pair<float, float> estimate_binning_densities(std::vector<xyzall_surfel_t>& input_data, bool is_verbose){

		/*std::shared_ptr<lamure::ren::bvh> bvh = std::make_shared<lamure::ren::bvh>( lamure::ren::bvh(bvh_filename) );
		scm::math::vec3f upper_bound_coord = bvh->get_bounding_box(0).max_vertex();
		scm::math::vec3f lower_bound_coord = bvh->get_bounding_box(0).min_vertex();

		scm::math::vec4f homogen_upper_bound_coord{ upper_bound_coord[0], upper_bound_coord[1], upper_bound_coord[2], 1.0f };
		scm::math::vec4f homogen_lower_bound_coord{ lower_bound_coord[0], lower_bound_coord[1], lower_bound_coord[2], 1.0f };

		homogen_upper_bound_coord = inverse_rotation_mat * homogen_upper_bound_coord;
		homogen_lower_bound_coord = inverse_rotation_mat * homogen_lower_bound_coord;

		upper_bound_coord = scm::math::vec3f{ homogen_upper_bound_coord[0], homogen_upper_bound_coord[1], homogen_upper_bound_coord[2] };
		lower_bound_coord = scm::math::vec3f{ homogen_lower_bound_coord[0], homogen_lower_bound_coord[1], homogen_lower_bound_coord[2] };*/

		/*if(is_verbose){
			std::cout << "Model bounding box corners \n";
			std::cout << "\t - min: (" << lower_bound_coord.x << ", " << lower_bound_coord.y << ", " << lower_bound_coord.z << ") \n";
			std::cout << "\t - max: (" << upper_bound_coord.x << ", " << upper_bound_coord.y << ", " << upper_bound_coord.z << ") \n";
		}*/

		/*float bb_dims[3] {-1.0f, -1.0f, -1.0f};
		float max_length(0.0f);
		float min_length(std::numeric_limits<float>::max());

		for (int dim_idx = 0; dim_idx < 3; ++dim_idx) {
			bb_dims[dim_idx] = std::fabs(upper_bound_coord[dim_idx] - lower_bound_coord[dim_idx]);
			max_length = std::max(max_length, bb_dims[dim_idx]);
			min_length = std::min(min_length, bb_dims[dim_idx]);
		}*/
		uint32_t last_el = input_data.size() - 1;
		float slicing_axis_length =  std::fabs(input_data[0].pos_coordinates[1] - input_data[last_el].pos_coordinates[1]);
		float min_distance_between_two_bins = 0.05 * slicing_axis_length;
		float max_distance_between_two_bins = 0.1 * slicing_axis_length; 

		return std::make_pair(min_distance_between_two_bins, max_distance_between_two_bins);
	}

	std::pair<float, float> estimate_binning_densities(std::string bvh_filename, bool is_verbose){
		std::shared_ptr<lamure::ren::bvh> bvh = std::make_shared<lamure::ren::bvh>( lamure::ren::bvh(bvh_filename) );
	}

} //namespace utils
} //namespace npr
