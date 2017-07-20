#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <lamure/types.h>

#include <cstdint>
#include <vector>

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

	using clusters_t = std::vector<point>;
	using bins_t = std::vector<xyzall_surfel_t>;

	namespace clustering {

		float compute_distance(lamure::vec3f const& pos1, lamure::vec3f const& pos2);
		std::pair<float, point> find_nearest_neighbour (point const& start_point, std::vector<point> const& all_points);
		std::vector<point*> find_near_neighbours(float max_distance, point const& start_point, std::vector<point>& all_points);
		void expand_cluster(point const& current_point, 
							std::vector<point*>& neighbourhood, 
							std::vector<point>& all_points_per_layer, 
							clusters_t& current_cluster, 
							 float eps, 
							 uint8_t min_points);
		std::vector<clusters_t>  create_DBSCAN_clusters (bins_t& all_surfels_per_layer, float eps, uint8_t min_points);
		std::vector<clusters_t>  create_clusters (bins_t& all_surfels_per_layer);
		   
    }


#endif //CLUSTERING_H