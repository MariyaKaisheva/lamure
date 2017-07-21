// Copyright (c) 2014 Bauhaus-Universitaet Weimar
// This Software is distributed under the Modified BSD License, see license.txt.
//
// Virtual Reality and Visualization Research Group 
// Faculty of Media, Bauhaus-Universitaet Weimar
// http://www.uni-weimar.de/medien/vr

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include <limits>
#include <random>

#include <lamure/ren/model_database.h>
#include <lamure/bounding_box.h>
#include <lamure/types.h>
#include <lamure/ren/dataset.h>
#include <lamure/ren/bvh.h>
#include <lamure/ren/lod_stream.h>

#include "clustering.h"
#include "utils.h"
#include "color_hash_map.hpp"
#include "math_wrapper.h"
#include "nurbscurve.hpp"
#include "point.hpp"

#define DEFAULT_PRECISION 15
#define OUTPUT_OBJ 1
#define BIDIRECTIONAL 0

class constrained_polyfit;

char* get_cmd_option(char** begin, char** end, const std::string & option) {
    char** it = std::find(begin, end, option);
    if (it != end && ++it != end)
        return *it;
    return 0;
}

bool cmd_option_exists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}


struct line{
    point start;
    point end;
    float length;
};


//using line_data_t = std::vector<line>;

//sort in descending order based on y coordinate value 
bool comparator (const xyzall_surfel_t& A, const xyzall_surfel_t& B) {
    return A.pos_coordinates[1] < B.pos_coordinates[1];
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
    std::cout << "Total num lines < 1 \n"; 
    return avg_line_length; 
  }
}


uint32_t current_cluster_id = 0;

//using namespace gpucast::math;
std::vector<line> generate_lines_from_curve (std::vector<point> const& cluster_of_points) {

  //std::cout << "incoming cluster size: " << cluster_of_points.size() << std::endl;

  //sort cluster points in correct order
  bool use_euclidian_distance = true;
  auto ordered_points = utils::order_points(cluster_of_points, use_euclidian_distance);
 
  for(uint duplication_idx = 0; duplication_idx < 7; ++duplication_idx) {
    ordered_points.push_back(ordered_points[duplication_idx]);
  }
  //coppy cluster content to vector of control points

  std::vector<gpucast::math::point3d> control_points_vec(ordered_points.size());
  
  #pragma omp parallel for
  for (uint cluster_index = 0; cluster_index < ordered_points.size(); ++cluster_index){
    auto x = ordered_points.at(cluster_index).pos_coordinates_[0];
    auto y = ordered_points.at(cluster_index).pos_coordinates_[1];
    auto z = ordered_points.at(cluster_index).pos_coordinates_[2];
    auto weight = 1.0;
    auto current_control_point = gpucast::math::point3d(x, y, z, weight);
    control_points_vec[cluster_index] = current_control_point;
  }

  //std::cout << "control_points_vec size: " << control_points_vec.size() << std::endl;

   uint degree = 6;

   //num control points must be >= order (degree + 1)
   if (control_points_vec.size() < degree + 1) {
      //throw std::runtime_error("Insufficient number of control points");
      degree = control_points_vec.size() - 1;
   }

  //generate knot vector
  std::vector<double> knot_vec;
 // auto knot_vec_size = control_points_vec.size() + degree + 1; //knot vector should be of size: num. control points + order
  
  float last_knot_value = 0;
  for(uint knot_counter = 0; knot_counter < control_points_vec.size(); ++knot_counter){ //
      knot_vec.push_back(double(knot_counter));
      last_knot_value = knot_counter;
  } 
  for(uint i = 0; i <= degree; ++i){ 
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
  float sampling_step = 0.5;
  //int ppoint_id_counter = 0;
  while(parameter_t < last_knot_value - sampling_step){

    auto st_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
    parameter_t += sampling_step;
    auto end_sampled_curve_point = nurbs_curve.evaluate(parameter_t);
    

    float pos_st_point[3] = {st_sampled_curve_point[0], st_sampled_curve_point[1], st_sampled_curve_point[2]};
    float pos_end_point[3] = {end_sampled_curve_point[0], end_sampled_curve_point[1], end_sampled_curve_point[2]};

    point current_start_point(pos_st_point);
    point current_end_point(pos_end_point);
    line current_line;
    current_line.start = current_start_point;
    current_line.end = current_end_point;
    current_line.length = utils::compute_distance(lamure::vec3f(current_line.start.pos_coordinates_[0], current_line.start.pos_coordinates_[1], current_line.start.pos_coordinates_[2]),
                                           lamure::vec3f(current_line.end.pos_coordinates_[0], current_line.end.pos_coordinates_[1], current_line.end.pos_coordinates_[2])); 
    line_segments_vec.push_back(current_line);
    //parameter_t += sampling_step;
  }

  //return sampled points as vector of line segments (lines)
  return line_segments_vec;
}

std::vector<line> generate_lines(std::vector<xyzall_surfel_t>& input_data, unsigned long number_line_loops, bool use_nurbs, bool apply_naive_clustering){
    
    //sort input points according to their y-coordinate 
    std::sort(input_data.begin(), input_data.end(), comparator);
    lamure::vec3f direction_ref_vector (1.0, 0.0, 0.0);
    float threshold = 0.008; //TODO think of alternative for dynamic calculation of thershold value

    std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size());
    std::vector<line> line_data;

    auto num_elements = input_data.size();
    auto height = input_data.at(num_elements-1).pos_coordinates[1] - input_data.at(0).pos_coordinates[1];
    //std::cout << "height: " << height << std::endl; 
    float const offset = height / number_line_loops;
    int total_num_clusters = 0;

    float current_y_min = input_data.at(0).pos_coordinates[1];
    float current_y_max = current_y_min + offset;

    for(uint i = 0; i < number_line_loops; ++i) {
       
        float current_y_mean = (current_y_min + current_y_max) / 2.0;
        //std::cout << "current_y_mean: " << current_y_mean << std::endl;
        if(threshold > current_y_max - current_y_mean) {
          throw  std::runtime_error("density thershold might be too low");
        } 
        auto copy_lambda = [&]( xyzall_surfel_t const& surfel){return (surfel.pos_coordinates[1] >= current_y_mean - threshold) && (surfel.pos_coordinates[1] <= current_y_mean + threshold);};
        auto it = std::copy_if(input_data.begin(), input_data.end(), current_bin_of_surfels.begin(), copy_lambda);
        current_bin_of_surfels.resize(std::distance(current_bin_of_surfels.begin(), it));
        current_y_min = current_y_max;
        current_y_max += offset;
        
        for (auto& surfel : current_bin_of_surfels) {
          surfel.pos_coordinates[1] = current_y_mean; //project all surfels for a given y_mean value to a single plane
        }
        

        std::vector<clusters_t> all_clusters_per_bin_vector;
        if(apply_naive_clustering){
          all_clusters_per_bin_vector = clustering::create_clusters(current_bin_of_surfels);
        }else{
          float eps = 0.08;
          uint8_t minPots = 3;
          //std::cout << "XX_C_XX STARTING DB SCAN WITH " << current_bin_of_surfels.size() << " surfels\n";
          all_clusters_per_bin_vector = clustering::create_DBSCAN_clusters(current_bin_of_surfels, eps, minPots);
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

        std::cout << "all_clusters_per_bin_vector.size(): " <<  all_clusters_per_bin_vector.size() << " \n";
        line current_line; 


        if(all_clusters_per_bin_vector.size() > 0){
          //
          std::random_device rd;
          std::mt19937 g(rd());

          for(auto& current_cluster : all_clusters_per_bin_vector){
            if(current_cluster.size() > 1) { //at least 2 vertices per cluster are need for one complete line 
             // std::cout << "Cluster size " << current_cluster.size() << std::endl;

              //sampling step (reduction of points to chosen % of original amout of points per cluster)
              float const sampling_rate = 0.7; //persentage points to remain after sampling
              unsigned const min_number_points_in_cell = 2;

              float min_x = std::numeric_limits<float>::max();
              float max_x = std::numeric_limits<float>::lowest();
              float min_z = std::numeric_limits<float>::max();
              float max_z = std::numeric_limits<float>::lowest();
              for (auto & current_point : current_cluster){
                min_x = std::min(min_x, current_point.pos_coordinates_[0]);
                max_x = std::max(max_x, current_point.pos_coordinates_[0]);
                min_z = std::min(min_z, current_point.pos_coordinates_[2]);
                max_z = std::max(max_z, current_point.pos_coordinates_[2]);
              }             
       
              auto num_cells_x_direction = 10;
              auto num_cells_z_direction = 10;
              auto cell_length = (max_z - min_z) / num_cells_z_direction;
              auto cell_width = (max_x - min_x) / num_cells_x_direction;

              std::vector<std::vector<point*>> vector_of_cells (num_cells_x_direction * num_cells_z_direction);

              //split all cluter points into respective grid cells
              for (auto & current_point : current_cluster){ //TODO check constness
                auto x_index = std::min(num_cells_x_direction - 1, std::max(0, int( (current_point.pos_coordinates_[0] - min_x) / cell_width)));
                auto z_index = std::min(num_cells_z_direction - 1, std::max(0, int( (current_point.pos_coordinates_[2]- min_z) / cell_length)));
                int64_t cell_index = z_index * num_cells_x_direction + x_index;
                auto& current_cell = vector_of_cells[cell_index];
                current_cell.push_back(&current_point);
              }
              auto num_points_in_cluster = current_cluster.size();
              int total_num_remaining_points = num_points_in_cluster * sampling_rate;
              std::vector<point> sampled_cluster;
              for (auto& current_cell : vector_of_cells) {

                std::shuffle(current_cell.begin(), current_cell.end(), g);
                uint64_t num_point_to_remain_in_cell = current_cell.size()*sampling_rate;

                std::cout << "Size before: " << current_cell.size() << "\n";

                if( min_number_points_in_cell < current_cell.size() ) {
                  current_cell.resize(num_point_to_remain_in_cell);
                  current_cell.shrink_to_fit(); 
                } else {
                  std::cout << "Did not reduce num points of this cell. " << "\n";
                }
                std::cout << "Size after: " << current_cell.size() << std::endl;
                for(auto& current_point : current_cell){
                  sampled_cluster.push_back(*current_point);
                }
              }

              std::cout << sampled_cluster.size() << " vs " << total_num_remaining_points << "\n";

              auto ordered_sample_cluster = utils::order_points(sampled_cluster, true);
              if(!use_nurbs){ 

                uint32_t num_lines_to_push = (ordered_sample_cluster.size()) - 1;

                for (uint line_idx = 0; line_idx < num_lines_to_push; ++line_idx) { 
                  current_line.start = ordered_sample_cluster.at(line_idx);
                  current_line.end = ordered_sample_cluster.at(line_idx+1);
                  current_line.length = utils::compute_distance(lamure::vec3f(current_line.start.pos_coordinates_[0], current_line.start.pos_coordinates_[1], current_line.start.pos_coordinates_[2]),
                                                          lamure::vec3f(current_line.end.pos_coordinates_[0], current_line.end.pos_coordinates_[1], current_line.end.pos_coordinates_[2]));
                  line_data.push_back(current_line);
                }
              }else {
                std::vector<line> const line_data_from_sampled_curve = generate_lines_from_curve(sampled_cluster);

                for(auto& current_line : line_data_from_sampled_curve){

                  line_data.push_back(current_line);
                } 
              }
              ++total_num_clusters;
            }
          }
        }
        else{
          std::cout << "no clusters in the current layer \n"; 
        }

    }   
    //std::cout << "num Lines: " << line_data.size() << std::endl;
    std::cout << "num line loops: " << total_num_clusters << std::endl;
    return line_data;   
}

int main(int argc, char *argv[]) {
    
    if (argc == 1 || cmd_option_exists(argv, argv+argc, "-h") ||
        !cmd_option_exists(argv, argv+argc, "-f")) {
        
      std::cout << "Usage: " << argv[0] << " -f <input_file>" << std::endl <<
         "Parameters: " << std::endl <<
         "\t-f: specify .bvh input file" << std::endl <<
         "\t-d: (optional) specify depth to extract; default value is the maximal depth, i.e. leaf level" << std::endl <<
         "\t-l: (optional) specify number of slycing layers; default value is 20" << std::endl <<
         "\t--bin_density_threshold: (optional) specify max. distance to mean y_coordinate value; implicitly set surfel selection tolerance; default value is 0.002" << std::endl << 
         "\t--apply_nurbs_fitting: (optional); set flag for curve-fitting to TRUE" << std::endl << 
         "\t--use_dbscan: (optional); set DBSCAN as prefered clustering algorithm" << std::endl << 
         "\t--write_xyz_points: (optional) writes an xyz_point_cloud instead of a *.obj containing line data" << std::endl  <<
         std::endl;
      return 0;
    }


    std::string bvh_filename = std::string(get_cmd_option(argv, argv + argc, "-f"));
    std::string ext = bvh_filename.substr(bvh_filename.size()-3);
    if (ext.compare("bvh") != 0) {
        std::cout << "please specify a .bvh file as input" << std::endl;
        return 0;
    }

    lamure::ren::bvh* bvh = new lamure::ren::bvh(bvh_filename);
    int32_t depth = -1;
    if(cmd_option_exists(argv, argv+argc, "-d")){
      depth = atoi(get_cmd_option(argv, argv+argc, "-d"));
      if(depth > int(bvh->get_depth()) || depth < 0){
        depth = bvh->get_depth();
      }
    }
    else{
      depth = bvh->get_depth();
    }

    bool write_obj_file = !cmd_option_exists(argv, argv + argc, "--write_xyz_points");

    unsigned long number_line_loops = 25; //TODO make number dependent on model size??
    if(cmd_option_exists(argv, argv+argc, "-l")){
     number_line_loops = atoi(get_cmd_option(argv, argv+argc, "-l")); //user input
    }
    
    std::string obj_filename = bvh_filename.substr(0, bvh_filename.size()-4)+ "_d" + std::to_string(depth) + "_l" + std::to_string(number_line_loops) + ".obj";
    std::string xyz_all_filename = bvh_filename.substr(0, bvh_filename.size()-4)+ "_d" + std::to_string(depth) + "_l" + std::to_string(number_line_loops) + ".xyz_all";
    std::cout << "input: " << bvh_filename << std::endl;
   
    std::cout << "output: ";
    if(write_obj_file) {
      std::cout << obj_filename; 
    }
    else {
      std::cout << xyz_all_filename;
    }
    std::cout << std::endl;

    std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
    lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
    in_access->open(lod_filename);

    size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);

    std::cout << "working with surfels at depth " << depth << "; Max depth for this moidel is "<<  bvh->get_depth() <<std::endl;
    lamure::node_t num_leafs = bvh->get_length_of_depth(depth);
 

    auto total_num_surfels = num_leafs*bvh->get_primitives_per_node();
    std::vector<xyzall_surfel_t> surfels_vector(total_num_surfels);

    auto num_offset_nodes = std::pow(2, depth) - 1; //only valid for fanfactor 2
    auto length_in_bytes = num_leafs * size_of_node;
    in_access->read((char*) &surfels_vector[0], num_offset_nodes * size_of_node, length_in_bytes);
  
  
    //clean data from degenerated padding surfels
    surfels_vector.erase(std::remove_if(surfels_vector.begin(), 
                                        surfels_vector.end(),
                                        [](xyzall_surfel_t s){return s.radius_ <= 0.0;}),
                         surfels_vector.end());
    bool use_nurbs = cmd_option_exists(argv, argv + argc, "--apply_nurbs_fitting");
    bool apply_naive_clustering = !cmd_option_exists(argv, argv + argc, "--use_dbscan");
    auto line_data = generate_lines(surfels_vector, number_line_loops, use_nurbs, apply_naive_clustering);
    
    #if 0
    std::cout << "Num lines BEFORE clean up: " << line_data.size() << std::endl;
    //clean data
    auto avg_line_length = compute_global_average_line_length(line_data); 
    line_data.erase(std::remove_if(line_data.begin(),
                                   line_data.end(),
                                   [&](line l){return l.length >= 4.8 * avg_line_length;}),
                    line_data.end());
    std::cout << "Num lines AFTER clean up: " << line_data.size() << std::endl;
    #endif

    if(write_obj_file) {
      std::ofstream output_file(obj_filename);
      unsigned long vert_counter = 1;

      if (output_file.is_open()){
          for (uint i = 0; i < line_data.size(); ++i){

           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).start.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).start.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION)<< line_data.at(i).start.pos_coordinates_[2] << "\n";
           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[2] << "\n";
           //vertex duplication to emulate triangle
           auto x_offset =  line_data.at(i).end.pos_coordinates_[0] / 1000000.0;
           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.pos_coordinates_[0] + x_offset<< " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[2] << "\n";
           output_file << "f " << vert_counter << " " << (vert_counter + 1) << " " << (vert_counter + 2) << "\n";
           vert_counter += 3;
          }

          //TODO connect first and last line segment of the loop
        
          output_file.close();
      }
      else{
        std::cout << "Cannot open output file to write to! \n";
      }

    }
    else {
      //consider hidden translation
      const scm::math::vec3f& translation = bvh->get_translation();
      std::ofstream output_file(xyz_all_filename);
      lamure::vec3f const fixed_upward_normal(0.0, 1.0, 0.0);
      lamure::vec3f const fixed_forward_normal(0.0, 0.0, 1.0);
      float const fixed_radius(0.02);
      if (output_file.is_open()){
          for (uint i = 0; i < line_data.size(); ++i){
           output_file << std::setprecision(DEFAULT_PRECISION) << translation.x + line_data.at(i).start.pos_coordinates_[0] << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << translation.y + line_data.at(i).start.pos_coordinates_[1] << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << translation.z + line_data.at(i).start.pos_coordinates_[2] << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_upward_normal.x << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_upward_normal.y << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_upward_normal.z << " ";
           output_file << (int) line_data.at(i).start.r_  << " ";
           output_file << (int) line_data.at(i).start.g_ << " ";
           output_file << (int) line_data.at(i).start.b_ << " ";  
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_radius << std::endl;
           #if BIDIRECTIONAL
           output_file << std::setprecision(DEFAULT_PRECISION) << translation.x + line_data.at(i).start.pos_coordinates_[0] << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << translation.y + line_data.at(i).start.pos_coordinates_[1] << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << translation.z + line_data.at(i).start.pos_coordinates_[2] << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_forward_normal.x << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_forward_normal.y << " ";
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_forward_normal.z << " ";
           output_file << (int) line_data.at(i).start.r_  << " ";
           output_file << (int) line_data.at(i).start.g_ << " ";
           output_file << (int) line_data.at(i).start.b_ << " ";  
           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_radius << std::endl;
           #endif
          }
        
          output_file.close();
      }
      else{
        std::cout << "Cannot open output file to write to! \n";
      }
    }
   
    std::cout << "ok \n";

    delete in_access;
    delete bvh;

    return 0;
}



