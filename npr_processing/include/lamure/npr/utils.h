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
			point(float pos_x, float pos_y, float pos_z) : pos_coordinates_{pos_x, pos_y, pos_z},  id_(0), r_(200), g_(0), b_(250), is_used_(false), is_member_of_cluster_(false) {}
			point(float* pos) : pos_coordinates_{pos[0], pos[1], pos[2]},  id_(0), r_(200), g_(0), b_(250), is_used_(false), is_member_of_cluster_(false) {}
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

namespace npr {
namespace utils {
			inline float compute_distance(lamure::vec3f const& pos1, lamure::vec3f const& pos2) {
				lamure::vec3f distance_vector((pos1.x - pos2.x), (pos1.y - pos2.y), (pos1.z - pos2.z));
				float result = std::sqrt(distance_vector.x*distance_vector.x +
				                    	 distance_vector.y*distance_vector.y +
				                    	 distance_vector.z*distance_vector.z);
				return result;
			}
		}


		struct line{
			line() : start(point()), end(point()), length(0.0f){}
			line(point const& start_point, point const& end_point, float length) : start(start_point), end(end_point), length(length) {}
			line(point const& start_point, point const& end_point) : start(start_point), end(end_point) {
				auto l = utils::compute_distance(lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]),
                                           		lamure::vec3f(end_point.pos_coordinates_[0], end_point.pos_coordinates_[1], end_point.pos_coordinates_[2]));
				length = l;
			}

		    point start;
		    point end;
		    float length;
		};

		struct bounding_rect{ //TODO substitute by scm::gl::boxf
			float min_x;
			float max_x;
			float min_y;
			float max_y;
			float min_z;
			float max_z;
		};


		class grid_cell
		{
			public:
				float min_corner_coord_[3];
				float max_corner_coord_[3];
				std::vector<xyzall_surfel_t>	content_;
				bool compute_intersection(const float* sphere_origin, float radius) const{
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
	  			//std::cout << "Total num lines < 1 \n"; 
	  			return avg_line_length; 
			}
		}

		inline void find_candidate_cells(std::vector<grid_cell>& all_cells, const float* sphere_origin, float radius, std::vector<grid_cell*>& out_candidates) {
			//std::vector<grid_cell> candidates;

			for(auto& curent_cell : all_cells) {
				if(curent_cell.content_.size() > 0) {
					if(curent_cell.compute_intersection(sphere_origin, radius)){
						out_candidates.push_back(&curent_cell);
					}
				}
			}
			//std::cout << candidates.size() << " - num candidate cells \n";
			//return candidates;
		}

		inline bounding_rect compute_bounding_corners (std::vector<xyzall_surfel_t> const& input_data){
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


		inline void generate_cells(std::vector<xyzall_surfel_t> const& input_data, 
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



	       	//std::vector<grid_cell> cells_vec(num_cells_pro_dim * num_cells_pro_dim * num_cells_pro_dim);

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

	       	//return cells_vec;
		}

		inline bool test_for_sufficency(std::vector<grid_cell*> const& input_cells) {
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

		float compute_avg_min_distance(std::vector<xyzall_surfel_t> const& input_data, uint32_t const num_cells_pro_dim_x, uint32_t const num_cells_pro_dim_y, uint32_t const num_cells_pro_dim_z);


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
		}

		void transform(std::vector<line>& line_data_vec, scm::math::mat4f const& transformation_mat);

		std::pair<float, float> estimate_binning_densities(std::string const& bvh_filename);

} //utils
} //npr



#endif //UTILS_H