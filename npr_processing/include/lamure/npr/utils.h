#ifndef UTILS_H
#define UTILS_H

#include <lamure/types.h>

#include <iostream>
#include <vector>
#include <limits>

	//general structs
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

//utils specific structs
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


		struct grid_cell
		{

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


	//free functions
	namespace utils{

	//compute average min distance between for surfels in provided data set (parallelized with OMP)
	float compute_avg_min_distance(std::vector<xyzall_surfel_t> const& input_data, 
								   uint32_t const num_cells_pro_dim_x, uint32_t const num_cells_pro_dim_y, uint32_t const num_cells_pro_dim_z);

	//compute 3D bb of set of input surfels
	bounding_rect compute_bounding_corners (std::vector<xyzall_surfel_t> const& input_data);

	//compute centroid position for set of 3D points
	lamure::vec3f compute_cluster_centroid_position (std::vector<point> const& point_cluster);

	//compute average line length for set of 3D input lines
	float compute_global_average_line_length(std::vector<line> const& all_lines);

	//computes dot product of 2 vectors
	float dot(lamure::vec3f const& vec_A, lamure::vec3f const& vec_B);

	std::pair<float, float> estimate_binning_densities(std::vector<xyzall_surfel_t>& input_data, bool is_verbose);
	std::pair<float, float> estimate_binning_densities(std::string bvh_filename, bool is_verbose);

	//find cells intersected by a 3D sphere
	void find_candidate_cells(std::vector<grid_cell>& all_cells, const float* sphere_origin, float radius, std::vector<grid_cell*>& out_candidates);

	//generate a 3D grid with fixed number of cells along input surfels
	void generate_cells(std::vector<xyzall_surfel_t> const& input_data, 
							   float& suggested_search_radius, 
							   std::vector<grid_cell>& out_cells, 
							   int const num_cells_pro_dim_x, 
							   int const num_cells_pro_dim_y,
							   int const num_cells_pro_dim_z);

	//normalize vec3f
	lamure::vec3f normalize (lamure::vec3f const& in_vector);

	//sort points in ascending order
	std::vector<point> order_points(std::vector<point> const& cluster_of_points, bool euclidean_distance);

	//helper function to check if sufficient candidate cells were found
	bool test_for_sufficency(std::vector<grid_cell*> const& input_cells);

	//transform line vector vertices to by provided transformation matrix
	void transform(std::vector<line>& line_data_vec, scm::math::mat4f const& transformation_mat);


} //utils
} //npr



#endif //UTILS_H