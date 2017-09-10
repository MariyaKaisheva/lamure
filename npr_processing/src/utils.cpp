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

	//compute 3D bb of set of input surfels
	bounding_rect compute_bounding_corners (std::vector<xyzall_surfel_t> const& input_data){
			float min_x = std::numeric_limits<float>::max();
			float max_x = std::numeric_limits<float>::lowest();
	       	float min_y = std::numeric_limits<float>::max();
	        float max_y = std::numeric_limits<float>::lowest();
	        float min_z = std::numeric_limits<float>::max();
	        float max_z = std::numeric_limits<float>::lowest();

			//compute corner points for the gird structure 
	        for (auto const& current_point : input_data){
	          min_x = std::min(min_x, current_point.pos_coordinates[0]);
	          max_x = std::max(max_x, current_point.pos_coordinates[0]);
	          min_y = std::min(min_y, current_point.pos_coordinates[1]);
	          max_y = std::max(max_y, current_point.pos_coordinates[1]);
	          min_z = std::min(min_z, current_point.pos_coordinates[2]);
	          max_z = std::max(max_z, current_point.pos_coordinates[2]);
	        } 

	        bounding_rect corner_coordinates;
	        corner_coordinates.min_x = min_x;
	        corner_coordinates.max_x = max_x;
	        corner_coordinates.min_y = min_y;
	        corner_coordinates.max_y = max_y;
	        corner_coordinates.min_z = min_z;
	        corner_coordinates.max_z = max_z;

	        return corner_coordinates;
	}

	lamure::vec3f compute_cluster_centroid_position (std::vector<point> const& point_cluster) {
	  float average_x = 0.0;
	  float average_y = 0.0;
	  float average_z = 0.0;
	  for(auto const& point : point_cluster){
	    average_x += point.pos_coordinates_[0];
	    average_y += point.pos_coordinates_[1];
	    average_z += point.pos_coordinates_[2];
	  }
	  auto number_of_surfels_per_layer = point_cluster.size();
	  average_x /= number_of_surfels_per_layer;
	  average_y /= number_of_surfels_per_layer;
	  average_z /= number_of_surfels_per_layer;
	  lamure::vec3f average_position(average_x, average_y, average_z);
	  
	  return average_position;
	}

	float compute_global_average_line_length(std::vector<line> const& all_lines) {
		float avg_line_length = 0.0; 
		if(all_lines.size() > 1){
  			for (auto const& line : all_lines){
    			avg_line_length += line.length; 
  			}
  			return avg_line_length/all_lines.size();     
		}
		else{
  			//std::cout << "Total num lines < 1 \n"; 
  			return avg_line_length; 
		}
	}

	//computes dot product of 2 vectors
	float dot(lamure::vec3f const& vec_A, lamure::vec3f const& vec_B) {
	  return vec_A.x*vec_B.x + vec_A.y*vec_B.y + vec_A.z*vec_B.z;
	}

	//find cells intersected by a 3D sphere
	void find_candidate_cells(std::vector<grid_cell>& all_cells, const float* sphere_origin, float radius, std::vector<grid_cell*>& out_candidates) {
		for(auto& curent_cell : all_cells) {
			if(curent_cell.content_.size() > 0) {
				if(curent_cell.compute_intersection(sphere_origin, radius)){
					out_candidates.push_back(&curent_cell);
				}
			}
		}
	}

	//generate a 3D grid with fixed number of cells along input surfels
	void generate_cells(std::vector<xyzall_surfel_t> const& input_data, 
							   float& suggested_search_radius, 
							   std::vector<grid_cell>& out_cells, 
							   int const num_cells_pro_dim_x, 
							   int const num_cells_pro_dim_y,
							   int const num_cells_pro_dim_z) {

        auto bounding_rect = compute_bounding_corners(input_data);
        float min_x = bounding_rect.min_x;
		float max_x = bounding_rect.max_x;
		float min_y = bounding_rect.min_y;
		float max_y = bounding_rect.max_y;
		float min_z = bounding_rect.min_z;
		float max_z = bounding_rect.max_z;

        float cell_width = fabs(max_x - min_x) / num_cells_pro_dim_x;
        float cell_height = fabs(max_y - min_y) / num_cells_pro_dim_y;
        float cell_depth = fabs(max_z - min_z) / num_cells_pro_dim_z;
       //one forth of shortest cell edge

        float comparison_values[3] = {cell_width, cell_height, cell_depth};

        suggested_search_radius = std::numeric_limits<float>::max();

        for(uint32_t dim_idx = 0; dim_idx < 3; ++dim_idx) {
       	  if(0.0 == comparison_values[dim_idx]) {
            comparison_values[dim_idx] = std::numeric_limits<float>::max();
       	  }

       	  suggested_search_radius = std::min(suggested_search_radius, comparison_values[dim_idx]);
       }

       if( std::numeric_limits<float>::max() == suggested_search_radius ) {
       	suggested_search_radius = 1.0;
       }

       suggested_search_radius /= 4.0;


        size_t num_cells = num_cells_pro_dim_x * num_cells_pro_dim_y * num_cells_pro_dim_z;
        out_cells.reserve(num_cells);

        for( size_t z_idx = 0; z_idx < num_cells_pro_dim_z; ++z_idx) {
        	for( size_t y_idx = 0; y_idx < num_cells_pro_dim_y; ++y_idx) {
        		for( size_t x_idx = 0; x_idx < num_cells_pro_dim_x; ++x_idx) {

        			float cell_min_x = min_x + x_idx * cell_width;
        			float cell_max_x = cell_min_x + cell_width;

        			float cell_min_y = min_y + y_idx * cell_height;
        			float cell_max_y = cell_min_y + cell_height;

        			float cell_min_z = min_z + z_idx * cell_depth;
        			float cell_max_z = cell_min_z + cell_depth;

        		    grid_cell current_cell_to_position;
        		    current_cell_to_position.min_corner_coord_[0] = cell_min_x;
        		    current_cell_to_position.min_corner_coord_[1] = cell_min_y;
        		 	current_cell_to_position.min_corner_coord_[2] = cell_min_z;

        		    current_cell_to_position.max_corner_coord_[0] = cell_max_x;
        		    current_cell_to_position.max_corner_coord_[1] = cell_max_y;
        		 	current_cell_to_position.max_corner_coord_[2] = cell_max_z;

        		    out_cells.push_back(current_cell_to_position);
        		}
        	}      	
        }

       	for(auto& current_surfel : input_data){
       	  int32_t x_index = 0; 
       	  int32_t y_index = 0;
          int32_t z_index = 0;

          if( 0.0 != cell_width ) {
          	x_index = std::min(num_cells_pro_dim_x - 1, std::max(0, int32_t( (current_surfel.pos_coordinates[0] - min_x) / cell_width)));
      	  }

          if( 0.0 != cell_height ) {
          	y_index = std::min(num_cells_pro_dim_y - 1, std::max(0, int32_t( (current_surfel.pos_coordinates[1] - min_y) / cell_height)));
      	  }

          if( 0.0 != cell_depth ) {
          	z_index = std::min(num_cells_pro_dim_z - 1, std::max(0, int32_t( (current_surfel.pos_coordinates[2] - min_z) / cell_depth)));
      	  }

          //std::cout << "x_index " << x_index << " y_index " << y_index << " z_index " << z_index << "\n";
          
          int64_t cell_index = z_index * num_cells_pro_dim_z * num_cells_pro_dim_y + y_index * num_cells_pro_dim_y  + x_index;
          auto& current_cell = out_cells[cell_index];
          current_cell.content_.push_back(current_surfel);
       	}
	}

	//normalize lamure lamure::vec3f
	lamure::vec3f normalize (lamure::vec3f const& in_vector) {
	  auto vector_length = sqrt(in_vector.x*in_vector.x + 
	                            in_vector.y*in_vector.y + 
	                            in_vector.z*in_vector.z);
	  return in_vector / vector_length;
	}


	//sort points in ascending order (=> DEPRECATED)
	std::vector<point> order_points(std::vector<point> const& cluster_of_points, bool euclidean_distance) {
		auto input_size = cluster_of_points.size();
		std::vector<point> ordered_result(input_size);
		std::vector<int32_t> unused_nearest_neighbor_of(input_size, -1); //point_used_status_vec (input_size, false);
		uint nearest_point_index = 0;

		int32_t current_point_index = -1;

		for (uint point_counter = 0; point_counter < input_size; ++point_counter) {

			float current_shortest_distance = std::numeric_limits<float>::max();

			ordered_result[point_counter] = (cluster_of_points[nearest_point_index]);
			current_point_index = nearest_point_index;
			auto& current_point = cluster_of_points[nearest_point_index];
			auto start_coordinates = lamure::vec3f(current_point.pos_coordinates_[0], current_point.pos_coordinates_[1], current_point.pos_coordinates_[2]);

			for (uint cluster_iterator = 0; cluster_iterator < input_size; ++cluster_iterator){
				if( -1 == unused_nearest_neighbor_of[cluster_iterator]){
					if( (point_counter != cluster_iterator) && (cluster_iterator != unused_nearest_neighbor_of[current_point_index]) ) {
						
						auto& next_point = cluster_of_points[cluster_iterator];
						auto end_coordinates = lamure::vec3f(next_point.pos_coordinates_[0], next_point.pos_coordinates_[1], next_point.pos_coordinates_[2]);
						auto distance_to_next_point = compute_distance(start_coordinates, end_coordinates);

						if(current_shortest_distance > distance_to_next_point) {
								current_shortest_distance = distance_to_next_point;
								nearest_point_index = cluster_iterator;
							}
					}
				}
			}
			unused_nearest_neighbor_of[nearest_point_index] = current_point_index; // make the current point idx the predecessor of the found nearest point
		}


		return ordered_result;
	}

	//helper function to check if sufficient candidate cells were found
	bool test_for_sufficency(std::vector<grid_cell*> const& input_cells) {
		uint content_counter = 0;
		for (auto const& current_cell : input_cells){
			content_counter += current_cell->content_.size();
		}

		if (content_counter > 1){
			return false;
		}else{
			return true;
		}
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
