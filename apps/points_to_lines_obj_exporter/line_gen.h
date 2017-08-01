#ifndef LINE_GEN_H
#define LINE_GEN_H

#include "clustering.h"


//sort in descending order based on y coordinate value 
bool comparator (const xyzall_surfel_t& A, const xyzall_surfel_t& B) {
    return A.pos_coordinates[1] < B.pos_coordinates[1];
}

inline std::vector<line> generate_lines_from_curve (std::vector<point> const& ordered_points) {

  //std::cout << "incoming cluster size: " << cluster_of_points.size() << std::endl;
 
  //coppy cluster content to vector of control points

  std::vector<gpucast::math::point3d> control_points_vec(ordered_points.size());
  
  #pragma omp parallel for
  for (uint32_t cluster_index = 0; cluster_index < ordered_points.size(); ++cluster_index){
    auto x = ordered_points.at(cluster_index).pos_coordinates_[0];
    auto y = ordered_points.at(cluster_index).pos_coordinates_[1];
    auto z = ordered_points.at(cluster_index).pos_coordinates_[2];
    auto weight = 1.0;
    auto current_control_point = gpucast::math::point3d(x, y, z, weight);
    control_points_vec[cluster_index] = current_control_point;
  }

  //std::cout << "control_points_vec size: " << control_points_vec.size() << std::endl;

   uint32_t degree = 4;

   //num control points must be >= order (degree + 1)
   if (control_points_vec.size() < degree + 1) {
      //throw std::runtime_error("Insufficient number of control points");
      degree = control_points_vec.size() - 1;
   }

  //generate knot vector
  std::vector<double> knot_vec;
 // auto knot_vec_size = control_points_vec.size() + degree + 1; //knot vector should be of size: num. control points + order
  
  float last_knot_value = 0;
  for(uint32_t knot_counter = 0; knot_counter < control_points_vec.size(); ++knot_counter){ //
      knot_vec.push_back(double(knot_counter));
      last_knot_value = knot_counter;
  } 
  for(uint32_t i = 0; i <= degree; ++i){ 
      knot_vec.push_back(double(last_knot_value));
  } 

  //curve fitting
  gpucast::math::nurbscurve3d nurbs_curve;
  nurbs_curve.set_points(control_points_vec.begin(), control_points_vec.end());
  nurbs_curve.degree(degree);
  nurbs_curve.set_knotvector(knot_vec.begin(), knot_vec.end());

  //sample the curve inside the knot span
  std::vector<line> line_segments_vec;
  float parameter_t = degree;
  float sampling_step = 0.6;
  //int ppoint_id_counter = 0;
  while(parameter_t < last_knot_value - sampling_step){

    auto st_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
    parameter_t += sampling_step;
    auto end_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
    

    float pos_st_point[3] = {st_sampled_curve_point[0], st_sampled_curve_point[1], st_sampled_curve_point[2]};
    float pos_end_point[3] = {end_sampled_curve_point[0], end_sampled_curve_point[1], end_sampled_curve_point[2]};

    point current_start_point(pos_st_point);
    point current_end_point(pos_end_point);
    line current_line(current_start_point, current_end_point);
    /*current_line.start = current_start_point;
    current_line.end = current_end_point;
    current_line.length = utils::compute_distance(lamure::vec3f(current_line.start.pos_coordinates_[0], current_line.start.pos_coordinates_[1], current_line.start.pos_coordinates_[2]),
                                           lamure::vec3f(current_line.end.pos_coordinates_[0], current_line.end.pos_coordinates_[1], current_line.end.pos_coordinates_[2])); */
    line_segments_vec.push_back(current_line);
    //parameter_t += sampling_step;
  }

  //return sampled points as vector of line segments (lines)
  return line_segments_vec;
}

inline std::vector<line> generate_lines(std::vector<xyzall_surfel_t>& input_data, uint32_t& max_num_line_loops, bool use_nurbs, bool apply_alpha_shapes){
	uint32_t current_cluster_id = 0;
	uint8_t order = 5; //TODO cosider changeing this variable to user-defined one

	//parameters used for distance-based sampling; TODO make them use input dependent if still used in the long run
    uint32_t max_num_points = 40;
    bool naive_sampling = true;

    //sort input points according to their y-coordinate 
    std::sort(input_data.begin(), input_data.end(), comparator);
    lamure::vec3f direction_ref_vector (1.0, 0.0, 0.0);


    //inital global computation for whole model 
    //with size to holding appyimately 1000 data points (assuming uniform point distribution)
    uint8_t num_cells_pro_dim =  ceil(cbrt(input_data.size() / 1000)); 
    float avg_min_distance = utils::compute_avg_min_distance(input_data, num_cells_pro_dim, num_cells_pro_dim, num_cells_pro_dim);
    float distance_threshold =  avg_min_distance * 4.0; 

    std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size());

    std::vector<line> line_data;
 

    //adaptive binning (distributes input surfels into descrite num. bins and projects them onto 2d plane)
    auto bins_vec = binning::generate_all_bins(input_data, distance_threshold, max_num_line_loops);

    for (auto const& current_bin_of_surfels : bins_vec){    

        std::vector<clusters_t> all_clusters_per_bin_vector;

        //compute average minimal distance between surfels in a bin
        float avg_min_distance_per_bin = utils::compute_avg_min_distance(current_bin_of_surfels.content_, num_cells_pro_dim, 1, num_cells_pro_dim);

        //set parameters for DBSCAN
		float eps = avg_min_distance_per_bin * 20.0; // radis of search area
		uint8_t minPoints = 3; //minimal number of data points that should be located inside search ared

		//generate clusters 
		all_clusters_per_bin_vector = clustering::create_DBSCAN_clusters(current_bin_of_surfels.content_, eps, minPoints);
		std::cout << "all_clusters_per_bin: " <<  all_clusters_per_bin_vector.size() << " \n";


        //color clusters
        for( auto& current_cluster : all_clusters_per_bin_vector ) {
          uint32_t current_cluster_color_id = id_to_color_hash(current_cluster_id);
          lamure::vec3b current_cluster_color = color_array[current_cluster_color_id];
          for( auto& point_in_current_cluster : current_cluster) {
            point_in_current_cluster.set_color(current_cluster_color);
          }
          ++current_cluster_id;
        }

        line current_line; 

        //create oulines for underlying shape of a cluster 
        if(all_clusters_per_bin_vector.size() > 0){

          for(auto& current_cluster : all_clusters_per_bin_vector){
            uint32_t cluster_size = current_cluster.size();


            if(cluster_size > order) { //at least 2 vertices per cluster are need for one complete line 
              std::cout << "Cluster size " << cluster_size << std::endl;
        		
              //no cluster reduction yet => all cluster members are kept;
              auto& sampled_cluster = current_cluster;

        	  //remove unnecessary points in a cluster
              if(naive_sampling){ //use  max-distance approach to select predifinned number of cluster members
              	std::cout << "applying naive sampling with " << max_num_points << "\n";
              	sampled_cluster = sampling::apply_distance_optimization_sampling (current_cluster, max_num_points); 
              }else if(apply_alpha_shapes){ //use alpha_shapes to select only points that make up the concave hull of a cluster
				float cluster_y_coord = current_cluster[0].pos_coordinates_[1];
				std::vector<alpha::cgal_point_2> cgal_points(cluster_size);
				alpha::do_input_conversion(cgal_points, current_cluster);
				auto cgal_line_segments = alpha::generate_alpha_shape(cgal_points);
				sampled_cluster =  alpha::do_output_conversion(cgal_line_segments, cluster_y_coord);
              }

              //sort cluster content
              auto ordered_cluster = utils::order_points(sampled_cluster, true);


              if(!use_nurbs){ //create straightforward line segments  
                uint32_t num_lines_to_push = (ordered_cluster.size()) - 1;

                for (uint line_idx = 0; line_idx < num_lines_to_push; ++line_idx) { 
                  current_line.start = ordered_cluster.at(line_idx);
                  current_line.end = ordered_cluster.at(line_idx+1);
                  current_line.length = utils::compute_distance(lamure::vec3f(current_line.start.pos_coordinates_[0], current_line.start.pos_coordinates_[1], current_line.start.pos_coordinates_[2]),
                                                          lamure::vec3f(current_line.end.pos_coordinates_[0], current_line.end.pos_coordinates_[1], current_line.end.pos_coordinates_[2]));
                  line_data.push_back(current_line);
                }
              }else {//fit NURBS curve and evaluate it 
                std::vector<line> const line_data_from_sampled_curve = generate_lines_from_curve(ordered_cluster);
                for(auto& current_line : line_data_from_sampled_curve){
                  line_data.push_back(current_line);
                }
              }
            }
          }
        }
        else{
          std::cout << "no clusters in the current layer \n"; 
        }
    }
  
    return line_data;
}

#endif //LINE_GEN_H