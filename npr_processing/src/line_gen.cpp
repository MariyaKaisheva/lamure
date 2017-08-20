#include <lamure/npr/line_gen.h>

#include <lamure/npr/binning.h>
#include <lamure/npr/color_hash_map.h>
#include <lamure/npr/sampling.h>
#include <lamure/npr/alpha_shapes_wrapper.h>

namespace npr {
namespace line_gen{

//sort in descending order based on y coordinate value 
bool 
comparator (const xyzall_surfel_t& A, const xyzall_surfel_t& B) {
    return A.pos_coordinates[1] < B.pos_coordinates[1];
}

gpucast::math::point3d 
get_midpoint(gpucast::math::point3d const& start_point, gpucast::math::point3d const& end_point) {
	auto x = (start_point[0] + end_point[0]) / 2.0;
  	auto y = (start_point[1] + end_point[1]) / 2.0;
  	auto z = (start_point[2] + end_point[2]) / 2.0;
  	auto middle_point = gpucast::math::point3d(x, y, z, 1);
  	return middle_point;
}

void 
interpolate_cluster_input(std::vector<point> & combined_cluster_points, std::vector<uint32_t> const& original_cluster_sizes){
  auto num_clusters = original_cluster_sizes.size();
  auto first_point_index = 0;
  auto last_point_index = original_cluster_sizes[0] + original_cluster_sizes[1] - 1;


  for(int cluster_index = 0; cluster_index < num_clusters - 1; ++cluster_index){
    uint32_t num_points_lower_cluster = original_cluster_sizes[cluster_index];
    uint32_t num_points_upper_cluster = original_cluster_sizes[cluster_index + 1];
    uint32_t total_num_points_in_current_cluster_subsection = (num_points_lower_cluster + num_points_upper_cluster) / 2;
    float total_height = fabs(combined_cluster_points[first_point_index].pos_coordinates_[1] - combined_cluster_points[last_point_index].pos_coordinates_[1]);

    float base_offset = total_height / total_num_points_in_current_cluster_subsection;
    float height_offset = base_offset;

    for (uint32_t pnt_counter = num_points_lower_cluster/2; pnt_counter < num_points_lower_cluster; ++pnt_counter) {
      combined_cluster_points[first_point_index + pnt_counter].pos_coordinates_[1] += height_offset;
      height_offset += base_offset;
    }

    height_offset = ((num_points_upper_cluster - 1) - (num_points_upper_cluster / 2) + 1) * base_offset;

    for (uint32_t pnt_counter = num_points_upper_cluster - 1; pnt_counter > num_points_upper_cluster / 2; --pnt_counter) {
      //std::cout << height_offset << " index for access: " << last_point_index - pnt_counter << "; size of vecs: " << combined_cluster_points.size() << "\n";
      uint32_t backward_iterating_index = last_point_index - pnt_counter;
      combined_cluster_points[backward_iterating_index].pos_coordinates_[1] -= height_offset;
      height_offset -= base_offset;
    }
    first_point_index += original_cluster_sizes[cluster_index];
    last_point_index += original_cluster_sizes[cluster_index + 1];
  }
}

std::vector<line> 
 generate_lines_from_curve (std::vector<point> const& ordered_points, uint8_t degree, bool is_verbose) {
 
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


   //connect first to last segment
   auto first_point = control_points_vec[0];
   control_points_vec.push_back(first_point);

   //num control points must be >= order (degree + 1)
   /*if (control_points_vec.size() < degree + 1) {
      throw std::runtime_error("Insufficient number of control points");
      degree = control_points_vec.size() - 1;
   }*/
/*
     std::cout << "Control point x_coord: [";
  for(auto cp :  control_points_vec){

      std::cout << "{" << cp.x << ", " << cp.y << ", " << cp.z << "}, ";
    
  }

  std::cout << "]\n\n";
*/
  //generate knot vector
  std::vector<double> knot_vec;
 // auto knot_vec_size = control_points_vec.size() + degree + 1; //knot vector should be of size: num. control points + order
  
  //all leading zeroes except for one
  for(uint32_t i = 0; i < degree; ++i){ 
      knot_vec.push_back(double(0.0));
  }

  //sequence 0 .. #control_points - degree 
  float last_knot_value = 0;
  for(uint32_t knot_counter = 0; knot_counter < control_points_vec.size() + 1 - degree ; ++knot_counter){ //
      knot_vec.push_back(double(knot_counter));
      last_knot_value = knot_counter;
  } 

  //all trailing #control_points except for the first one
  for(uint32_t i = 0; i < degree; ++i){ 
      knot_vec.push_back(double(last_knot_value));
  } 


  //std::cout << "num knots: " << knot_vec.size() << "\n";
  //std::cout << "num control points: " << control_points_vec.size() << "\n"; 

  /*std::cout << "Knot Vec: [";
  for(auto knot :  knot_vec){
     std::cout << knot << " ";
  }

  std::cout << "]\n\n";*/

  //curve fitting
  gpucast::math::nurbscurve3d nurbs_curve;
  nurbs_curve.set_points(control_points_vec.begin(), control_points_vec.end());
  nurbs_curve.degree(degree);
  nurbs_curve.set_knotvector(knot_vec.begin(), knot_vec.end());

  //sample the curve inside the knot span
  std::vector<line> line_segments_vec;
  
  float evaluation_offset = 0.0001;



  #if 1 //dynamic sampling step parameter
 	//intial points
  	float initial_t = evaluation_offset; 
  	float final_t = last_knot_value;
  	
  	auto start_point = nurbs_curve.evaluate(initial_t);
  	auto end_point = nurbs_curve.evaluate(final_t);
  	auto approximated_middle_point = get_midpoint(start_point, end_point); 


  	//uint8_t const inital_vec_size = 2;
  	std::vector<gpucast::math::point3d> evaluated_points_vec(2);
  	evaluated_points_vec.push_back(start_point);
  	evaluated_points_vec.push_back(end_point);

  	std::stack<line_approximation_job> working_stack;

  	//create inital approximation jobs
  	line_approximation_job j_0(initial_t, final_t, approximated_middle_point);
  	working_stack.push(j_0);
 

  	float error_threshold = 0.001;
  	while(!working_stack.empty()){
      //std::cout << "still in the loop!\n";
  		auto current_job = working_stack.top();
  		working_stack.pop();

  		auto start_point = nurbs_curve.evaluate(current_job.start_t_);
  		auto end_point = nurbs_curve.evaluate(current_job.end_t_);
	  	float middle_t = (current_job.start_t_ + current_job.end_t_) / 2.0;
	  	auto middle_point = nurbs_curve.evaluate(middle_t);
	  	auto & approx_middle_point = current_job.approximated_point_;
	  	auto start_coord = lamure::vec3f(middle_point[0], middle_point[1], middle_point[2]);
	  	auto end_coord = lamure::vec3f(approx_middle_point[0], approx_middle_point[1], approx_middle_point[2]);
	  	auto error = utils::compute_distance(start_coord, end_coord);
	  	if (error > error_threshold) {
	  		//prepare 2 new jobs
	  		auto approximated_middle_point_1 = get_midpoint(start_point, nurbs_curve.evaluate(middle_t));
	  		line_approximation_job j_1(current_job.start_t_, middle_t, approximated_middle_point_1);

  			auto approximated_middle_point_2 = get_midpoint(nurbs_curve.evaluate(middle_t), end_point);
  			line_approximation_job j_2(middle_t, current_job.end_t_, approximated_middle_point_2);

  			working_stack.push(j_1);
  			working_stack.push(j_2);

	  	}else{//push data for final storage
	  		//convert to point type used by line struk
	  		float start_coord_arr[3] = {start_point[0], start_point[1], start_point[2]};
	  		float end_coord_arr[3] = {end_point[0], end_point[1], end_point[2]};
	  		point start_point_(start_coord_arr);
	  		point end_point_(end_coord_arr);

	  		line closely_approximated_line(start_point_, end_point_);
	  		line_segments_vec.push_back(closely_approximated_line);
	  	}
  	}

  #else //fixed sampling step parameter

	  float parameter_t = evaluation_offset;
	  float sampling_step = 1.2;
	  
	  while(parameter_t < last_knot_value - sampling_step){

	    auto st_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
	    std::cout << "XXXXXX " << st_sampled_curve_point[0] << ", " << st_sampled_curve_point[1] << ", " << st_sampled_curve_point[2] << "\n";
      parameter_t += sampling_step;
	    auto end_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
	    std::cout << "t" << parameter_t << std::endl;


	    float pos_st_point[3] = {st_sampled_curve_point[0], st_sampled_curve_point[1], st_sampled_curve_point[2]};
	    float pos_end_point[3] = {end_sampled_curve_point[0], end_sampled_curve_point[1], end_sampled_curve_point[2]};

	    point current_start_point(pos_st_point);
	    point current_end_point(pos_end_point);
	    line current_line(current_start_point, current_end_point);
	    line_segments_vec.push_back(current_line);
	  }
  #endif

  if(is_verbose) {
    std::cout << "NUM line segments: " << line_segments_vec.size() << "\n";
  }

  //return sampled points as vector of line segments (lines)
  return line_segments_vec;
}

std::vector<line> 
generate_lines(std::vector<xyzall_surfel_t>& input_data, uint32_t& max_num_line_loops, bool use_nurbs, bool apply_alpha_shapes, bool is_verbose){
	uint32_t current_cluster_id = 0;
	uint8_t degree = 3; //TODO consider changeing this variable to user-defined one

	//parameters used for distance-based sampling; TODO make them use input dependent if still used in the long run
    uint32_t max_num_points = 40;
    bool naive_sampling = false;

    //sort input points according to their y-coordinate 
    std::sort(input_data.begin(), input_data.end(), comparator);
    lamure::vec3f direction_ref_vector (1.0, 0.0, 0.0);


    //inital global computation for whole model 
    //with size to holding appyimately 1000 data points (assuming uniform point distribution)
    uint8_t num_cells_pro_dim =  std::ceil(std::cbrt(input_data.size() / 1000)); 
    float avg_min_distance = utils::compute_avg_min_distance(input_data, num_cells_pro_dim, num_cells_pro_dim, num_cells_pro_dim);
    float distance_threshold =  avg_min_distance * 2.0; 

    std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size());

    std::vector<line> line_data;
 
     std::vector<point> combined_sampled_clusters;
     std::vector<uint32_t> cluster_sizes;
    //adaptive binning (distributes input surfels into descrite num. bins and projects them onto 2d plane)
    auto bins_vec = binning::generate_all_bins(input_data, distance_threshold, max_num_line_loops);

    for (auto const& current_bin_of_surfels : bins_vec){    

        std::vector<clusters_t> all_clusters_per_bin_vector;

        //compute average minimal distance between surfels in a bin
        float avg_min_distance_per_bin = utils::compute_avg_min_distance(current_bin_of_surfels.content_, num_cells_pro_dim, 1, num_cells_pro_dim);

        if(0.0f == avg_min_distance_per_bin) {
          continue;
        }

        //set parameters for DBSCAN
		float eps = avg_min_distance_per_bin * 6.0; // radis of search area
		uint8_t minPoints = 3; //minimal number of data points that should be located inside search ared

		//generate clusters 
		all_clusters_per_bin_vector = clustering::create_DBSCAN_clusters(current_bin_of_surfels.content_, eps, minPoints);

    if(is_verbose) {
		  std::cout << "all_clusters_per_bin: " <<  all_clusters_per_bin_vector.size() << " \n";
    }


        //color clusters
        for( auto& current_cluster : all_clusters_per_bin_vector ) {
          uint32_t current_cluster_color_id = id_to_color_hash(current_cluster_id);
          lamure::vec3b current_cluster_color = color_array[current_cluster_color_id];
          for( auto& point_in_current_cluster : current_cluster) {
            point_in_current_cluster.set_color(current_cluster_color);
          }
          ++current_cluster_id;
        }


       
        //create oulines for underlying shape of a cluster 
        if(all_clusters_per_bin_vector.size() > 0){

          for(auto& current_cluster : all_clusters_per_bin_vector){
            uint32_t cluster_size = current_cluster.size();


            //for nurbs fitting:  num controll points should be at least as much as the order of the curve
            //smaller glusters should be ignored 
            if(cluster_size > degree + 1) { 
              if(is_verbose) {
                std::cout << "Cluster size " << cluster_size << std::endl;
              }
        		
              //no cluster reduction yet => all cluster members are kept;
              auto& sampled_cluster = current_cluster;

        	  //remove unnecessary points in a cluster
              if(naive_sampling){ //use  max-distance approach to select predifinned number of cluster members
                if(is_verbose) {
              	 std::cout << "applying naive sampling with " << max_num_points << "\n";
                }
              	sampled_cluster = sampling::apply_distance_optimization_sampling (current_cluster, max_num_points); 
              }else if(apply_alpha_shapes){ //use alpha_shapes to select only points that make up the concave hull of a cluster
                float cluster_y_coord = current_cluster[0].pos_coordinates_[1];
                std::vector<alpha::cgal_point_2> cgal_points(cluster_size);
                alpha::do_input_conversion(cgal_points, current_cluster);
                auto cgal_line_segments = alpha::generate_alpha_shape(cgal_points);
                sampled_cluster =  alpha::do_output_conversion(cgal_line_segments, cluster_y_coord);
              }

              //sort cluster content (ONLY WHEN NOT USING ALPHA SHAPES)
              #if 0
                auto ordered_cluster = utils::order_points(sampled_cluster, true);
              #else //no sorting (ONLY WHEN USING ALPHA SHAPES)
                auto& ordered_cluster = sampled_cluster;
              #endif

              //recheck cluster size sufficiency after the potential reduction during alpha-shapes detection  
              if( ordered_cluster.size() <= degree + 1 ) {
                  if(is_verbose) {
                    std::cout << "cluster with " << ordered_cluster.size() << " points was skipped \n";
                  }
                  continue;
                }

              #if 1 
                  if(!use_nurbs){ //create straightforward line segments

                    uint32_t num_lines_to_push = (ordered_cluster.size()) - 1;

                    for (uint line_idx = 0; line_idx < num_lines_to_push; ++line_idx) {
                      
                      line_data.emplace_back(ordered_cluster.at(line_idx), ordered_cluster.at(line_idx+1));
                    }
                  }else {//fit NURBS curve and evaluate it 
    /*
                    float sqrt_of_val = std::sqrt(0.5);
                    std::vector<point> fake_cluster = {    
                                                          {1.0, 0, 0}, {sqrt_of_val, 0, sqrt_of_val}, {0, 0, 1}, 
                                                          {- sqrt_of_val, 0, sqrt_of_val}, {-1, 0, 0}, {-sqrt_of_val, 0, -sqrt_of_val},
                                                          {0, 0, -1}, {sqrt_of_val, 0, -sqrt_of_val}
                                                      };


                    std::vector<line> line_data_from_sampled_curve = generate_lines_from_curve(fake_cluster);
    */
                    /*if(line_data_from_sampled_curve.size() > 1) {
                      line_data_from_sampled_curve.push_back( line(line_data_from_sampled_curve[0].end, line_data_from_sampled_curve[line_data_from_sampled_curve.size()-1].start ) );
                    }*/

                    std::vector<line> line_data_from_sampled_curve = generate_lines_from_curve(ordered_cluster, degree);
                    for(auto& current_line : line_data_from_sampled_curve){
                      line_data.emplace_back(current_line);
                    }
                  }
              #endif
                //cluster_sizes.push_back(ordered_cluster.size());
                //combined_sampled_clusters.insert(combined_sampled_clusters.end(), ordered_cluster.begin(), ordered_cluster.end());
            }
          }
        }
        else{
          if(is_verbose) {
            std::cout << "no clusters in the current layer \n";
          }
        }
    }
    return line_data;
    interpolate_cluster_input(combined_sampled_clusters, cluster_sizes);

    if(is_verbose) {
      std::cout << "After interpolating cluster input\n";
      std::cout << "interpolate_cluster_input: " << combined_sampled_clusters.size() << "\n";
    }
    
    std::vector<line> combined_line_data = generate_lines_from_curve(combined_sampled_clusters);
    
    if(is_verbose) {
      std::cout << "After line from curve\n";
    }
   // return combined_line_data; 
}

} //namespace line_gen
} //namespace npr