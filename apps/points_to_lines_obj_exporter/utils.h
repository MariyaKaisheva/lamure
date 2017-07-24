#ifndef UTILS_H
#define UTILS_H

#include <lamure/types.h>

#include <iostream>
#include <vector>
#include <limits>

			//surfel layout as it is used for rendering 
		struct xyzall_surfel_t {
		  xyzall_surfel_t() : radius_(0.0) {}
		  float pos_coordinates[3];
		  uint8_t r_, g_, b_, fake_;
		  float radius_;
		  float nx_, ny_, nz_;
		};

		struct point{
			point() : pos_coordinates_{0.0, 0.0, 0.0}, id_(0), r_(0), g_(150), b_(50), is_used_(false), is_member_of_cluster_(false) {}
			point(float* pos) : pos_coordinates_{pos[0], pos[1], pos[2]},  id_(0), r_(200), g_(150), b_(50), is_used_(false), is_member_of_cluster_(false) {}                        
			point(float* pos, int32_t id) : pos_coordinates_{pos[0], pos[1], pos[2]},  id_(id), r_(200), g_(150), b_(50), is_used_(false), is_member_of_cluster_(false) {}   
			point(float* pos, int32_t id, uint8_t red, uint8_t green, uint8_t blue, bool is_used = false, bool is_member_of_cluster = false) : 
			                          pos_coordinates_{pos[0], pos[1], pos[2]},
			                          id_(id),
			                          r_(red),
			                          g_(green),
			                          b_(blue),
			                          is_used_(is_used),
			                          is_member_of_cluster_(is_member_of_cluster) {}

			bool operator==(point const& rhs){
			  return (*this).id_ == rhs.id_;
			}

			void set_color(lamure::vec3b const& color_vector) {
			  r_ = color_vector[0];
			  g_ = color_vector[1];
			  b_ = color_vector[2];
			}

			float pos_coordinates_[3];
			int32_t id_;
			uint8_t r_, g_, b_;
			bool is_used_;
			bool is_member_of_cluster_;
		};

		struct line{
		    point start;
		    point end;
		    float length;
		};


		class grid_cell
		{
			public:
				float min_corner_coord_[3];
				float max_corner_coord_[3];
				float wall_size_;
				std::vector<point>	content_;
				bool compute_intersection(float* sphere_origin, float radius){
					float squared_distance = 0.0;
					for(uint i = 0; i < 3; ++i){
						float temp_value = 0.0;
						if(sphere_origin[i] < min_corner_coord_[i] ){
							temp_value = (min_corner_coord_[i] - sphere_origin[i]) * 
										 (min_corner_coord_[i] - sphere_origin[i]);
						}
				         
				        if ( sphere_origin[i] > max_corner_coord_[i] )
				        {
							temp_value = (sphere_origin[i] - max_corner_coord_[i]) * 
										 (sphere_origin[i] - max_corner_coord_[i]);
				        }

						squared_distance += temp_value;
					}

					return squared_distance <= (radius * radius);
				}		
		};


	namespace utils{

		inline lamure::vec3f normalize (lamure::vec3f const& centroid_surfel_vec) {
		  auto vector_length = sqrt(centroid_surfel_vec.x*centroid_surfel_vec.x + 
		                            centroid_surfel_vec.y*centroid_surfel_vec.y + 
		                            centroid_surfel_vec.z*centroid_surfel_vec.z);
		  lamure::vec3f result = centroid_surfel_vec / vector_length;
		  return result;
		}

		//computes dot product of 2 vectors
		inline float dot(lamure::vec3f const& vec_A, lamure::vec3f const& vec_B) {
		  float product = vec_A.x*vec_B.x + vec_A.y*vec_B.y + vec_A.z*vec_B.z;
		  return product;
		}


		inline float compute_distance(lamure::vec3f const& pos1, lamure::vec3f const& pos2) {
			lamure::vec3f distance_vector((pos1.x - pos2.x), (pos1.y - pos2.y), (pos1.z - pos2.z));
			float result = sqrt(distance_vector.x*distance_vector.x +
			                    distance_vector.y*distance_vector.y +
			                    distance_vector.z*distance_vector.z);
			return result;
		}

		inline lamure::vec3f compute_cluster_centroid_position (std::vector<point> const& point_cluster) {
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

		inline float compute_global_average_line_length(std::vector<line> const& all_lines) {
			float avg_line_length = 0.0; 
			if(all_lines.size() > 1){
	  			for (auto const& line : all_lines){
	    			avg_line_length += line.length; 
	  			}
	  			return avg_line_length/all_lines.size();     
			}
			else{
	  			std::cout << "Total num lines < 1 \n"; 
	  			return avg_line_length; 
			}
		}

		inline float compute_average_min_point_distance(std::vector<xyzall_surfel_t> const& input_data){
			float average_min_distance = 0.0f;

			#pragma omp parallel for reduction( + : average_min_distance)
			for(uint32_t outer_point_idx = 0; outer_point_idx < input_data.size(); ++outer_point_idx) {
				auto const& start_point = input_data[outer_point_idx];
				float current_min_distance = std::numeric_limits<float>::max();
				for(uint32_t inner_point_idx = 0; inner_point_idx < input_data.size(); ++inner_point_idx) {
					auto const& next_point = input_data[inner_point_idx];
					if( outer_point_idx != inner_point_idx ) {
						auto current_distance = compute_distance(lamure::vec3f(start_point.pos_coordinates[0], start_point.pos_coordinates[1], start_point.pos_coordinates[2]),
																 lamure::vec3f(next_point.pos_coordinates[0], next_point.pos_coordinates[1], next_point.pos_coordinates[2]));
						current_min_distance = std::min(current_distance, current_min_distance);
					}

				}

				average_min_distance += current_min_distance;
			}



			return average_min_distance / input_data.size();

		}

		inline std::vector<grid_cell> find_candidate_cells(std::vector<grid_cell> const& all_cells, float* sphere_origin, float radius) {
			std::vector<grid_cell> candidates;

			for(auto const& curren_cell : all_cells) {
				if(curren_cell.content_.size() > 0) {
					if(curent_cell.compute_intersection(sphere_origin, radius)){
						candidates.push_back(curren_cell);
					}
				}
			}
			std::cout << candidates.size() << " - num candidate cells \n";
			return candidates;
		}

		/*inline std::vector<grid_cell> generate_cells(std::vector<xyzall_surfel_t> const& input_data) {
			float cell_size = cbrt(ceil(input_data.size() / 1000));
			float min_x = std::numeric_limits<float>::max();
	        float max_x = std::numeric_limits<float>::lowest();
	       	float min_y = std::numeric_limits<float>::max();
	        float max_y = std::numeric_limits<float>::lowest();
	        float min_z = std::numeric_limits<float>::max();
	        float max_z = std::numeric_limits<float>::lowest();
	        for (auto & current_point : input_data){
	          min_x = std::min(min_x, current_point.pos_coordinates[0]);
	          max_x = std::max(max_x, current_point.pos_coordinates[0]);
	          min_y = std::min(min_z, current_point.pos_coordinates[1]);
	          max_y = std::max(max_z, current_point.pos_coordinates[1]);
	          min_z = std::min(min_z, current_point.pos_coordinates[2]);
	          max_z = std::max(max_z, current_point.pos_coordinates[2]);
	        } 

	       uint8_t num_cells_x_direction = fabs(max_x - min_x)/cell_size;
	       uint8_t num_cells_y_direction = fabs(max_y - min_y)/cell_size;
	       uint8_t num_cells_z_direction = fabs(max_z - min_z)/cell_size;

	       	std::vector<grid_cell> cells_vec(num_cells_x_direction * num_cells_y_direction * num_cells_z_direction);

	       	for(auto& current_surfel : input_data){

	       	}

	       	return cells_vec;
		} */

		inline float compute_average_min_point_distance_gridbased(std::vector<xyzall_surfel_t> const& input_data){
			float average_min_distance = 0.0f;
			auto cells_vec = generate_cells(input_data);
			auto  cell_wall_size = cells_vec[0]->wall_size_;
			float search_radius = sqrt(2)*cell_wall_size;
			for(auto const& current_point : input_data){
				auto start_coord = lamure::vec3f(current_point.pos_coordinates[0], current_point.pos_coordinates[1], current_point.pos_coordinates[2]);
				float current_min_distance = std::numeric_limits<float>::max();
				auto candidate_cells = find_candidate_cells(cells_vec, current_point.pos_coordinates, search_radius);
				for (auto const& curent_cell : candidate_cells) {
					for(auto const& cell_point : curren_cell->content_){
						auto end_coord = lamure::vec3f(cell_point.pos_coordinates[0], cell_point.pos_coordinates[1], cell_point.pos_coordinates[2]);
						auto current_distance = compute_distance(start_coord, end_coord);
						min_distance = std::min(min_distance, current_distance);
					}

				} 
				average_min_distance += min_distance;
			}

			return average_min_distance / input_data.size();
		}

		inline std::vector<point> order_points(std::vector<point> const& cluster_of_points, bool euclidian_distance){
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
		};
	}

#endif //UTILS_H