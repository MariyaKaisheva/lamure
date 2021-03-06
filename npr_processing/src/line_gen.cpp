#include <lamure/npr/line_gen.h>

#include <lamure/npr/binning.h>
#include <lamure/npr/color_hash_map.h>
#include <lamure/npr/sampling.h>
#include <lamure/npr/alpha_shapes_wrapper.h>

#include <chrono>
 
namespace npr {
namespace line_gen{

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
    float total_height = std::fabs(combined_cluster_points[first_point_index].pos_coordinates_[1] - combined_cluster_points[last_point_index].pos_coordinates_[1]);

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

//order-preserving shift of vector elemets (i.e. cluster points) 
//new first element is cluster point with smallest x-coordinate value 
void rotate(clusters_t & point_cluster) {
  uint32_t new_first_index = 0;
  float current_smallest_x = std::numeric_limits<float>::max();
  for(uint32_t point_index = 0; point_index < point_cluster.size(); ++point_index){

    float current_x = point_cluster[point_index].pos_coordinates_[0];
    if(current_smallest_x >= current_x){
      new_first_index = point_index;
      current_smallest_x = current_x;
    }
  }

  std::rotate(point_cluster.begin(), point_cluster.begin() + new_first_index, point_cluster.end());

}


gpucast::math::nurbscurve3d fit_curve(std::vector<point> const& ordered_points, uint8_t degree, bool spiral_look, bool is_verbose){
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
  if(!spiral_look){
    auto first_point = control_points_vec[0];
    control_points_vec.push_back(first_point);
  }


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

  //curve fitting
  gpucast::math::nurbscurve3d nurbs_curve;
  nurbs_curve.set_points(control_points_vec.begin(), control_points_vec.end());
  nurbs_curve.degree(degree);
  nurbs_curve.set_knotvector(knot_vec.begin(), knot_vec.end());

  return nurbs_curve;
}

std::vector<line> evaluate_curve(gpucast::math::nurbscurve3d & nurbs_curve, bool dynamic_sampling_step) {
  
  if(0 == nurbs_curve.knots().size()) {
    return std::vector<line>();
  }
  
  nurbs_curve.normalize_knotvector();
  float last_knot_value = 1.0; 
  //sample the curve inside the knot span
  std::vector<line> line_segments_vec;
  
  float evaluation_offset = 0.00;

 if(dynamic_sampling_step){ //dynamic sampling step parameter
  //intial points
    float initial_t = evaluation_offset;
    float final_t = last_knot_value;
    //float final_t = 1.0;
    
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
 

    float error_threshold = 0.0001;
    while(!working_stack.empty()){
      //std::cout << "still in the loop! stack size: " << working_stack.size() << "\n";
      auto current_job = working_stack.top();
      working_stack.pop();

      auto start_point = nurbs_curve.evaluate(current_job.start_t_);
      auto end_point = nurbs_curve.evaluate(current_job.end_t_);
      float middle_t = (current_job.start_t_ + current_job.end_t_) / 2.0;
      auto middle_point = nurbs_curve.evaluate(middle_t);

      //nurbs_curve.print(std::cout);
      auto & approx_middle_point = current_job.approximated_point_;
      auto start_coord = lamure::vec3f(middle_point[0], middle_point[1], middle_point[2]);
      auto end_coord = lamure::vec3f(approx_middle_point[0], approx_middle_point[1], approx_middle_point[2]);
      auto error = utils::compute_distance(start_coord, end_coord);
      //std::cout << error << " error value \n"; 
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
  }else{ //fixed sampling step parameter

    float parameter_t = evaluation_offset;
    float sampling_step = 1.2;
    
    while(parameter_t < last_knot_value - sampling_step){

      auto st_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
      parameter_t += sampling_step;
      auto end_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
    
      float pos_st_point[3] = {st_sampled_curve_point[0], st_sampled_curve_point[1], st_sampled_curve_point[2]};
      float pos_end_point[3] = {end_sampled_curve_point[0], end_sampled_curve_point[1], end_sampled_curve_point[2]};

      point current_start_point(pos_st_point);
      point current_end_point(pos_end_point);
      line current_line(current_start_point, current_end_point);
      line_segments_vec.push_back(current_line);
    }
  }

  //TODO: 
  /*if(is_verbose) {
    std::cout << "NUM line segments: " << line_segments_vec.size() << "\n";
  }*/

  //return sampled points as vector of line segments (lines)
  return line_segments_vec;
}

void map_to_normalized_range(float& value_to_map, float const input_range_start, float const input_range_end) {
  value_to_map = ( 1.0 / (input_range_end - input_range_start)) * (value_to_map - input_range_start);
}

std::vector<point> blend_between_curves(gpucast::math::nurbscurve3d & top_curve, 
                                        gpucast::math::nurbscurve3d & bottom_curve,
                                        float  max_distance){

  //std::cout << "knot span Before: " << top_curve.get_knotspan() << std::endl;
  top_curve.normalize_knotvector();
  bottom_curve.normalize_knotvector();
  uint32_t num_points_per_winding = 100;
  //float evaluation_offset = 0.00001;

  //both curves have constat y-coordinate value at any sampling position 
  //height (y-coordinate) between the 2 curves gives is used to calculate mun windings 
  //auto top_curve_y = top_curve.evaluate(evaluation_offset)[1];
  //auto bottom_curve_y = bottom_curve.evaluate(evaluation_offset)[1];

  auto top_curve_y = top_curve.evaluate(1.0)[1];
  auto bottom_curve_y = bottom_curve.evaluate(1.0)[1];

  uint32_t windings = std::ceil(std::fabs(top_curve_y - bottom_curve_y) / max_distance);

  /*std::cout << "distance: " << top_curve_y - bottom_curve_y << std::endl;
  std::cout << "max_distance: " << max_distance << std::endl;
  std::cout << "Num windings: " << windings << std::endl;*/

  float sampling_step = (1.0) / num_points_per_winding;
  uint32_t num_blending_points = num_points_per_winding * windings;

  std::vector<point> blended_points(num_blending_points + 1);

  #pragma omp parallel for
  for (uint32_t point_counter = 0; point_counter <= num_blending_points; ++point_counter){

    float accumulated_t = point_counter * sampling_step; // + evaluation_offset;

    float fractional_of_accumulated_t = std::fmod(accumulated_t, 1.0);
    float sampling_parameter_t = fractional_of_accumulated_t; //std::max(evaluation_offset, fractional_of_accumulated_t);

    float blend_parameter = accumulated_t / float(windings);

    //map_to_normalized_range(blend_parameter, evaluation_offset / windings, (windings + evaluation_offset) / windings );

    auto top_curve_sample = top_curve.evaluate(sampling_parameter_t);
    auto bottom_curve_sample = bottom_curve.evaluate(sampling_parameter_t); 

    float blended_x = (1.0 - blend_parameter) * top_curve_sample[0] + blend_parameter * bottom_curve_sample[0];
    float blended_y = (1.0 - blend_parameter) * top_curve_sample[1] + blend_parameter * bottom_curve_sample[1];
    float blended_z = (1.0 - blend_parameter) * top_curve_sample[2] + blend_parameter * bottom_curve_sample[2];
    float pos_coordinates[3] = {blended_x, blended_y, blended_z};
           
    blended_points[point_counter] = point(pos_coordinates);
  }

  return blended_points; 
}

nurbs_vec_t generate_spirals(std::vector<nurbs_vec_t> const& guiding_nurbs_vec,
                             float max_distance){

  //nurbs_vec_t spiral_basis_nurbs_vec;
  std::vector<std::pair<nurbs_t, nurbs_t>> spiral_blending_pairs_vec;
  for(uint32_t bin_index = 0; bin_index < guiding_nurbs_vec.size() - 1; ++bin_index){
    auto & guiding_nurbs_in_current_bin = guiding_nurbs_vec[bin_index];
    auto & guiding_nurbs_in_adjacent_bin = guiding_nurbs_vec[bin_index + 1];
    uint32_t num_cluster_in_current_bin = guiding_nurbs_in_current_bin.size();
    uint32_t num_cluster_in_adjacent_bin = guiding_nurbs_in_adjacent_bin.size();

    std::vector<bool> remaining_available_clusters(num_cluster_in_adjacent_bin,true);



    for(uint32_t cluster_counter = 0; cluster_counter < num_cluster_in_current_bin; ++cluster_counter) {

        auto const& current_cluster_curve = guiding_nurbs_in_current_bin[cluster_counter];

        bool zeros = std::all_of(remaining_available_clusters.begin(), remaining_available_clusters.end(), [](bool i) { return i==false; });

        if(zeros) {
          std::cout << "WAAAAAAAH THIS IS BAD!\n";
          break;
        }
        uint32_t index_of_corresponding_curve = clustering::find_corresponding_cluster_curve(current_cluster_curve, 
                                                                                             guiding_nurbs_in_adjacent_bin,
                                                                                             remaining_available_clusters);
        auto const corresponding_cluster_curve_from_adjacent_bin = guiding_nurbs_in_adjacent_bin[index_of_corresponding_curve];

        if(num_cluster_in_current_bin <= num_cluster_in_adjacent_bin){
          
          spiral_blending_pairs_vec.push_back( std::make_pair(current_cluster_curve, corresponding_cluster_curve_from_adjacent_bin));
        } else {
          bool candidate_is_close = true; //TODO: additional test (e.g overlapping area) can improve connectivity to 
          if(candidate_is_close){
            spiral_blending_pairs_vec.push_back( std::make_pair(current_cluster_curve, corresponding_cluster_curve_from_adjacent_bin));
            remaining_available_clusters[index_of_corresponding_curve] = true;
          }

        }
    } 

  }
    

  bool spiral_look = true;
  uint32_t degree = 3; //TODO make this an input parameter
  nurbs_vec_t final_spiral_segments_vec(spiral_blending_pairs_vec.size());
  //for(uint8_t curve_index = 0; curve_index < spiral_basis_nurbs_vec.size() - 1; ++curve_index){
    
  #pragma omp parallel for
  for(uint8_t curve_index = 0; curve_index < spiral_blending_pairs_vec.size(); ++curve_index){
    std::vector<point> control_points_vec; 

    /*auto new_control_points = blend_between_curves(spiral_basis_nurbs_vec[curve_index],
                                                   spiral_basis_nurbs_vec[curve_index + 1],
                                                   max_distance);*/

    auto new_control_points = blend_between_curves(spiral_blending_pairs_vec[curve_index].first,
                                                   spiral_blending_pairs_vec[curve_index].second,
                                                   max_distance);

    control_points_vec.insert(control_points_vec.end(), new_control_points.begin(), new_control_points.end());


    auto final_spiral_curve = fit_curve(control_points_vec, degree, spiral_look, false);

    //final_spiral_curve.print(std::cout);

    final_spiral_segments_vec[curve_index] = final_spiral_curve;
    //final_spiral_segments_vec.push_back(final_spiral_curve);
  }



  return final_spiral_segments_vec;
}

bool with_detailed_prints = false;  //TODO clean print logic

void prepare_clusters (std::vector<binning::bin> const& bins_vec, 
                       std::vector< std::shared_ptr<std::vector<clusters_t>> > & all_clusters_per_bin_vector_for_all_slices,
                       float eps_factor,
                       uint32_t num_cells_pro_dim,
                       bool is_verbose){

    #pragma omp parallel for
    for(uint32_t bin_index = 0; bin_index < bins_vec.size(); ++bin_index){

         auto const& current_bin_of_surfels = bins_vec[bin_index]; 
        //compute average minimal distance between surfels in a bin
        float avg_min_distance_per_bin = utils::compute_avg_min_distance(current_bin_of_surfels.content_, num_cells_pro_dim, 1, num_cells_pro_dim);

        if(0.0f == avg_min_distance_per_bin) {
          if(is_verbose){std::cout << "Skipping EMPTY BIN \n";}
          continue;
        }

        //set parameters for DBSCAN
        float eps = avg_min_distance_per_bin * eps_factor; // radis of search area
        uint8_t minPoints = 3; //minimal number of data points that should be located inside search ared

        //generate clusters
        auto all_clusters_per_bin_vector = std::make_shared<std::vector<clusters_t>>(clustering::create_DBSCAN_clusters(current_bin_of_surfels.content_, eps, minPoints) );


        if(is_verbose && with_detailed_prints) {
          std::cout << "num clusters in current layer: " <<  all_clusters_per_bin_vector->size() << std::endl;
        }

        all_clusters_per_bin_vector_for_all_slices[bin_index] = all_clusters_per_bin_vector;
    }
}

void clean_clusters_via_alpha_shape_detection(std::vector< std::shared_ptr<std::vector<clusters_t>> > & all_clusters_per_bin_vector_for_all_slices,
                                              std::vector< std::shared_ptr<std::vector<clusters_t>> > & all_alpha_shapes_for_all_bins,
                                              uint32_t degree,
                                              bool radial_binning,
                                              scm::math::vec3f bounding_sphere_center,
                                              bool color,
                                              bool is_verbose){

  uint32_t current_cluster_id = 0;

  #pragma omp parallel for
  for(uint32_t bin_vec_index = 0; bin_vec_index < all_clusters_per_bin_vector_for_all_slices.size(); ++bin_vec_index){
    auto const all_clusters_per_bin_vector = all_clusters_per_bin_vector_for_all_slices[bin_vec_index]; 
    std::shared_ptr<std::vector<clusters_t>> alpha_shapes_in_current_bin = std::make_shared<std::vector<clusters_t>>();
    for(auto& current_cluster : *all_clusters_per_bin_vector){
      //color clusters for visualisation purposes
      if(color) {
        uint32_t current_cluster_color_id = id_to_color_hash(current_cluster_id);
        lamure::vec3b current_cluster_color = color_array[current_cluster_color_id];
        for( auto& point_in_current_cluster : current_cluster) {
          point_in_current_cluster.set_color(current_cluster_color);
        }
        ++current_cluster_id;
      }



      float cluster_y_coord = current_cluster[0].pos_coordinates_[1];
      if(radial_binning){
        cluster_y_coord = bounding_sphere_center.y;
      }
      uint32_t cluster_size = current_cluster.size();

      std::vector<alpha::cgal_point_2> cgal_points(cluster_size);
      alpha::do_input_conversion(cgal_points, current_cluster);
      auto cgal_line_segments = alpha::generate_alpha_shape(cgal_points);
      std::vector<point> ordered_cluster = alpha::do_output_conversion(cgal_line_segments, cluster_y_coord);


      //check if the cleaned and ordereed cluster still containes enough elemets for be later used for nurbs fitting 
      if( ordered_cluster.size() <= degree + 1 ) {
          if(is_verbose && with_detailed_prints) {
            std::cout << "cluster with " << ordered_cluster.size() << " points was skipped \n";
          }
          continue;
      }

      //order-preserving shift of vector elemets (i.e. cluster points) 
      //=> first element has similar position for each cluster to facilitate smooth thansition between clusters of differnt bin-layers
      rotate(ordered_cluster);

      alpha_shapes_in_current_bin->push_back(ordered_cluster);
    }

    all_alpha_shapes_for_all_bins[bin_vec_index] = alpha_shapes_in_current_bin;
  }

}

std::vector<line> 
generate_lines(std::vector<xyzall_surfel_t>& input_data, 
               line_generation_descriptor& line_gen_desc)
{

	uint32_t degree = 3; //TODO consider changeing this variable to user-defined one
  bool color = true; //TODO make dependent on --write_xyz_points

  uint32_t last_el = input_data.size() - 1;

  auto model_height = std::fabs(input_data[0].pos_coordinates[1] - input_data[last_el].pos_coordinates[1]);
  std::cout <<" Min Y from input data " <<input_data[0].pos_coordinates[1] << " Max Y from input data " << input_data[last_el].pos_coordinates[1] << std::endl;
  std::cout <<"---> Initial min distance " << line_gen_desc.min_distance_ << " +  model_height " << model_height << std::endl;

  uint32_t max_num_slicing_planes;
  if(!line_gen_desc.radial_slicing_){
    max_num_slicing_planes = std::floor(model_height / line_gen_desc.min_distance_);
  } else {
    max_num_slicing_planes = std::floor(360.0 / line_gen_desc.min_distance_);
  }
  

  //inital global computation for whole model 
  //with size to holding appyimately 1000 data points (assuming uniform point distribution)
  uint32_t num_cells_pro_dim =  std::ceil(std::cbrt(input_data.size() / 1000)); 


  std::chrono::time_point<std::chrono::system_clock> start_avg_min_distance_est, end_avg_min_distance_est;
  start_avg_min_distance_est = std::chrono::system_clock::now();
  float avg_min_distance = utils::compute_avg_min_distance(input_data, num_cells_pro_dim, num_cells_pro_dim, num_cells_pro_dim);
  end_avg_min_distance_est = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_avg_min_distance_est = end_avg_min_distance_est - start_avg_min_distance_est;

  //write the avg min distance for the calling core function to use it. NOTE: don't use out_avg_min_distance anymore after the next line
  line_gen_desc.out_avg_min_distance_ = avg_min_distance;

  float distance_threshold =  avg_min_distance * 2.0; 

  std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size());
  std::vector<line> line_data;

  //adaptive binning (distributes input surfels into descrite num. bins and projects them onto 2d plane)
  std::chrono::time_point<std::chrono::system_clock> start_binning, end_binning;
  start_binning = std::chrono::system_clock::now();


  scm::math::vec3f bounding_sphere_center = scm::math::vec3f(0.0, 0.0, 0.0); //just inial value which will be recomputed in case of radial binning mode
  auto bins_vec = binning::generate_all_bins(input_data, distance_threshold, 
                                             max_num_slicing_planes, line_gen_desc.radial_slicing_,
                                             line_gen_desc.bounding_sphere_transl_vec_,
                                             bounding_sphere_center,
                                             line_gen_desc.max_distance_, line_gen_desc.is_verbose_);

  end_binning = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_binning = end_binning - start_binning;

  std::vector<nurbs_vec_t> guiding_nurbs_vec;

  
  //clustering 
  std::chrono::time_point<std::chrono::system_clock> start_clustering, end_clustering;
  start_clustering = std::chrono::system_clock::now();

  std::vector< std::shared_ptr<std::vector<clusters_t>> > all_clusters_per_bin_vector_for_all_bins(bins_vec.size());
  prepare_clusters(bins_vec, 
                   all_clusters_per_bin_vector_for_all_bins,
                   line_gen_desc.eps_factor_,
                   num_cells_pro_dim,
                   line_gen_desc.is_verbose_);

  end_clustering = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_clustering_single_bin = end_clustering - start_clustering;

  std::vector<uint32_t> bins_index_to_ignore;

  for(uint32_t bin_index = 0; bin_index < bins_vec.size(); ++bin_index) {
    auto & current_vector_of_clusters_per_bin =  all_clusters_per_bin_vector_for_all_bins.at(bin_index);
    if( current_vector_of_clusters_per_bin == nullptr || current_vector_of_clusters_per_bin->empty() ){
      bins_index_to_ignore.push_back(bin_index);
    }
  }

  std::reverse(bins_index_to_ignore.begin(), bins_index_to_ignore.end()); //swich order to remove bins with higher index first
  for(uint32_t bin_index : bins_index_to_ignore) {
    bins_vec.erase(bins_vec.begin() + bin_index);
    all_clusters_per_bin_vector_for_all_bins.erase(all_clusters_per_bin_vector_for_all_bins.begin() + bin_index);
  }


  if(line_gen_desc.write_intermediate_results_){
    bool use_binning_coloring = true;
    io::write_intermediate_result_out(line_gen_desc.output_base_name_ + "_BINNING.pob",
                                      avg_min_distance,
                                      all_clusters_per_bin_vector_for_all_bins,
                                      line_gen_desc.transformation_mat_,
                                      use_binning_coloring);
    std::cout << "Done with BINNING\n";
  }





  if(line_gen_desc.write_intermediate_results_){ 
    //project points y-value to center of bin
    if(!line_gen_desc.radial_slicing_){
      uint32_t bin_id = 0;
      for( auto& clusters_in_bin_ptr: all_clusters_per_bin_vector_for_all_bins) {

        float center_of_bin = bins_vec[bin_id].pos_along_slicing_axis_;

        auto& clusters_in_bin = *clusters_in_bin_ptr;
        for( auto& cluster : clusters_in_bin) {

          for( auto& point_in_cluster : cluster) {  
            point_in_cluster.pos_coordinates_[1] = center_of_bin;
          }
        }


        ++bin_id;
      }
    } else{
      //TODO

    }
 

    io::write_intermediate_result_out(line_gen_desc.output_base_name_ + "_CLUSTERING.pob",
                                      avg_min_distance,
                                      all_clusters_per_bin_vector_for_all_bins,
                                      line_gen_desc.transformation_mat_);
    std::cout << "Done with CLUSTERING\n";
  }


  //alpha-shapes detection
  std::chrono::time_point<std::chrono::system_clock> start_alpha_shaping, end_alpha_shaping;
  start_alpha_shaping = std::chrono::system_clock::now();//start timing

  std::vector< std::shared_ptr<std::vector<clusters_t>> >  all_alpha_shapes_for_all_bins(all_clusters_per_bin_vector_for_all_bins.size());
  clean_clusters_via_alpha_shape_detection(all_clusters_per_bin_vector_for_all_bins,
                                           all_alpha_shapes_for_all_bins,
                                           degree,
                                           line_gen_desc.radial_slicing_,
                                           bounding_sphere_center,
                                           color,
                                           line_gen_desc.is_verbose_);

  end_alpha_shaping = std::chrono::system_clock::now();//end timing
  std::chrono::duration<double> elapsed_seconds_AS_dection_per_cluster = end_alpha_shaping - start_alpha_shaping; 

  //Thansform data points to the original non AABB space
 if(line_gen_desc.radial_slicing_){
   for (uint32_t bin_index = 0; bin_index < all_alpha_shapes_for_all_bins.size(); ++bin_index) {
    auto const& current_bin = bins_vec[bin_index];
    auto const& transformation_mat = current_bin.radial_rotation_mat_;
    auto & vector_of_alpha_shapes_per_bin = all_alpha_shapes_for_all_bins[bin_index];
    if(vector_of_alpha_shapes_per_bin->size() > 0){
      for(auto & current_alpha_shape : *vector_of_alpha_shapes_per_bin) {
        utils::transform_points_by_matrix(current_alpha_shape, transformation_mat);
      }
    }
  }
 }


  if(line_gen_desc.write_intermediate_results_){
    io::write_intermediate_result_out(line_gen_desc.output_base_name_ + "_ALPHA_SHAPES.pob",
                                avg_min_distance, 
                                all_alpha_shapes_for_all_bins,
                                line_gen_desc.transformation_mat_);
    std::cout << "Done with ALPHA_SHAPES\n";
  }

  std::vector< std::shared_ptr< std::vector<std::vector<line>> >> lines_for_all_bins(all_alpha_shapes_for_all_bins.size());

  //nurbs fitting
  std::chrono::time_point<std::chrono::system_clock> start_nurbs_fitting, end_nurbs_fitting;
  start_nurbs_fitting = std::chrono::system_clock::now();//start timing

  guiding_nurbs_vec.resize(all_alpha_shapes_for_all_bins.size());

  #pragma omp parallel for
  for(uint32_t bin_index = 0; bin_index < all_alpha_shapes_for_all_bins.size(); ++bin_index){
   std::shared_ptr<std::vector<clusters_t>>& all_alpha_shapes_in_current_bin_vec = all_alpha_shapes_for_all_bins.at(bin_index);

    nurbs_vec_t guiding_nurbs_per_layer;
    std::vector<std::vector<line>> lines_per_alpha_shape_in_current_bin;
    for(clusters_t& current_alpha_shaped_cluster : *all_alpha_shapes_in_current_bin_vec ){
      auto cluster_approximating_curve = fit_curve(current_alpha_shaped_cluster, degree, line_gen_desc.spiral_look_, false);
      

      auto evaluated_curve_lines = evaluate_curve(cluster_approximating_curve, true);

      if( 0 != evaluated_curve_lines.size() ) {
        lines_per_alpha_shape_in_current_bin.push_back(evaluated_curve_lines);
      }
      guiding_nurbs_per_layer.push_back(cluster_approximating_curve);
    }

    guiding_nurbs_vec[bin_index] = (guiding_nurbs_per_layer);
    lines_for_all_bins[bin_index] = std::make_shared<std::vector<std::vector<line>>>( lines_per_alpha_shape_in_current_bin ) ;

  }

  end_nurbs_fitting = std::chrono::system_clock::now();//end timing
  std::chrono::duration<double> elapsed_seconds_nurbs_fitting = end_nurbs_fitting - start_nurbs_fitting;

  for(auto& lines_in_current_bin : lines_for_all_bins) {
    for(auto& lines_per_alpha_shape : *lines_in_current_bin) {
      line_data.insert(std::end(line_data), std::begin(lines_per_alpha_shape), std::end(lines_per_alpha_shape));
    }
  }


  //generate spirals timing
  std::chrono::time_point<std::chrono::system_clock> start_spiral_generation, end_spiral_generation;
  start_spiral_generation = std::chrono::system_clock::now();//start timing

  float max_winding_distance = line_gen_desc.max_distance_ / 3.0; //TODO check this for radial slicing mode
  if(line_gen_desc.spiral_look_){
    line_data.clear();
    std::vector<gpucast::math::nurbscurve3d> final_curves_vec = generate_spirals(guiding_nurbs_vec, max_winding_distance);
    bool adaptive = true; 
    for(auto& current_spiral_section : final_curves_vec){
      std::vector<line> line_data_from_sampled_curve = evaluate_curve(current_spiral_section, adaptive);
      line_data.insert(std::end(line_data), std::begin(line_data_from_sampled_curve), std::end(line_data_from_sampled_curve));
    }
  }
  
  end_spiral_generation = std::chrono::system_clock::now();//end timing
  std::chrono::duration<double> elapsed_seconds_spiral_generation = end_spiral_generation - start_spiral_generation;

  if(line_gen_desc.is_verbose_){
    std::cout << "\t --  Time LOG:  -- avg min dist estimation: "<< elapsed_seconds_avg_min_distance_est.count()    << "s\n";
    std::cout << "\t --  Time LOG:  -- binning: "                << elapsed_seconds_binning.count()                 << "s\n";
    std::cout << "\t --  Time LOG:  -- clustering: "             << elapsed_seconds_clustering_single_bin.count()   << "s\n";
    std::cout << "\t --  Time LOG:  -- detecting Alpha-shapes: " << elapsed_seconds_AS_dection_per_cluster.count()  << "s\n";
    std::cout << "\t --  Time LOG:  -- nurbs fitting: "          << elapsed_seconds_nurbs_fitting.count()           << "s\n";
    std::cout << "\t --  Time LOG:  -- spiral generation: "      << (line_gen_desc.spiral_look_ ? elapsed_seconds_spiral_generation.count() : 0) << "s\n";
  }


  return line_data;
}

} //namespace line_gen
} //namespace npr
