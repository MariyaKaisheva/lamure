#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <lamure/types.h>

#include <cstdint>
#include <vector>

#include "utils.h"
#include "nurbscurve.h"



	using clusters_t = std::vector<point>;
	using bins_t = std::vector<xyzall_surfel_t>;

	namespace npr {
	namespace clustering {

		std::pair<float, point> find_nearest_neighbour (point const& start_point, std::vector<point> const& all_points);
		std::vector<point*> find_near_neighbours(float max_distance, point const& start_point, std::vector<point>& all_points);
		void expand_cluster(point const& current_point, 
							std::vector<point*>& neighbourhood, 
							std::vector<point>& all_points_per_layer, 
							clusters_t& current_cluster, 
							float eps, 
							uint8_t min_points);
		std::vector<clusters_t>  create_DBSCAN_clusters (bins_t const& all_surfels_per_layer, float eps, uint8_t min_points);
		gpucast::math::nurbscurve3d find_corresponding_cluster_curve(gpucast::math::nurbscurve3d const& current_cluster_curve,
                                                         			 std::vector<gpucast::math::nurbscurve3d > const& guiding_nurbs_in_adjacent_bin, 
                                                         			 std::vector<bool> & remaining_available_clusters);
		std::vector<clusters_t>  create_clusters (bins_t const& all_surfels_per_layer);
		   
    } //namespace clustering
	} //namespace npr


#endif //CLUSTERING_H