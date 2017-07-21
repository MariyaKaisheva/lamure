#ifndef UTILS_H
#define UTILS_H

#include <lamure/types.h>

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