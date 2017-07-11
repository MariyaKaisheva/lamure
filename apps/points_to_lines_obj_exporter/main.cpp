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

#include <lamure/ren/model_database.h>
#include <lamure/bounding_box.h>
#include <lamure/types.h>
#include <lamure/ren/dataset.h>
#include <lamure/ren/bvh.h>
#include <lamure/ren/lod_stream.h>

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

using bins_t = std::vector<xyzall_surfel_t>;
using clusters_t = std::vector<point>;
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

lamure::vec3f normalize (lamure::vec3f const& centroid_surfel_vec) {
  auto vector_length = sqrt(centroid_surfel_vec.x*centroid_surfel_vec.x + 
                            centroid_surfel_vec.y*centroid_surfel_vec.y + 
                            centroid_surfel_vec.z*centroid_surfel_vec.z);
  lamure::vec3f result = centroid_surfel_vec / vector_length;
  return result;
}

//computes dot product of 2 vectors
float dot(lamure::vec3f const& vec_A, lamure::vec3f const& vec_B) {
  float product = vec_A.x*vec_B.x + vec_A.y*vec_B.y + vec_A.z*vec_B.z;
  return product;
}

float compute_distance(lamure::vec3f const& pos1, lamure::vec3f const& pos2) {

  lamure::vec3f distance_vector((pos1.x - pos2.x), (pos1.y - pos2.y), (pos1.z - pos2.z));
  float result = sqrt(distance_vector.x*distance_vector.x +
                      distance_vector.y*distance_vector.y +
                      distance_vector.z*distance_vector.z);
  return result;
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

std::pair<float, point> find_nearest_neighbour (point const& start_point, std::vector<point> const& all_points) {
  float current_min = std::numeric_limits<float>::max();
  point current_nearest_neighbour;

  bool is_valid = false;

  for (auto& p : all_points){
    if((!p.is_used_)){
      float distance = compute_distance(lamure::vec3f(p.pos_coordinates_[0], p.pos_coordinates_[1], p.pos_coordinates_[2]), 
                                        lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]));
      if( distance > 0.0 && distance <= current_min ){
        current_min = distance;
        current_nearest_neighbour = p;

        is_valid = true;
      }
    }
  }

  if( !is_valid ) {
    current_min = std::numeric_limits<float>::max();
  }

  std::pair<float, point> result(current_min, current_nearest_neighbour);
  return result;
}


std::vector<point*> find_near_neighbours(float max_distance, point const& start_point, std::vector<point>& all_points){
  std::vector<point*> neighbourhood;
  for(auto& p : all_points){

    float distance = compute_distance(lamure::vec3f(p.pos_coordinates_[0], p.pos_coordinates_[1], p.pos_coordinates_[2]), 
                                      lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]));
    if(distance <= max_distance ){
      neighbourhood.push_back(&p);
    }
  }
  
  return neighbourhood; 
}

std::vector<clusters_t> create_clusters (bins_t& all_surfels_per_layer){
   
  //create vector of points, containg all surfels per layer in point format
  std::vector<point> all_points_per_layer;
  int32_t id_counter = 0; 
  for (auto& s: all_surfels_per_layer){
    point current_point(s.pos_coordinates, id_counter, s.r_, s.g_, s.b_, false);
    all_points_per_layer.push_back(current_point);
    ++id_counter;
  }

  //auto centroid = compute_average_position_per_layer(all_surfels_per_layer);
  std::vector<point> point_cluster;
  std::vector<clusters_t> clusters_vector;
  float distance_threshold = 60.0; 
  
  float last_measured_distance = 0;  //TODO think about this value

  //int cluster_counter  = 0; 

  //distribute points to clusters based on distance to closest available neigbour point
  for(uint point_iterator = 0; point_iterator < all_points_per_layer.size(); ++point_iterator){
      int point_index = point_iterator;
     // bool next_point_is_used = false;
      do {
        point* current_point = &all_points_per_layer[point_index];
        if(!current_point->is_used_){
            point_cluster.push_back(*current_point);
            current_point->is_used_ = true;

            //reset distance record for each new cluster; no distance has been stored so far
            if(point_cluster.size() < 2){
              last_measured_distance = std::numeric_limits<float>::max() / (2.0 * distance_threshold); 
            }
            
            auto search_result =  find_nearest_neighbour(*current_point, all_points_per_layer);
            auto distance_to_nearest_unused_point = search_result.first; 
            auto next_point = search_result.second; //function find_nearest_neighbour() makes sure next point is not used so far
            
            
            if(distance_to_nearest_unused_point <= (last_measured_distance * distance_threshold)){  
              point_cluster.push_back(next_point);
             // next_point.is_used_ = true;
              last_measured_distance  = distance_to_nearest_unused_point; 
              point_index = next_point.id_;
            }
            else{
             // next_point.is_used_ = false;
              //push current cluster to common container and start a new
              clusters_vector.push_back(point_cluster);
             // ++cluster_counter; 
              point_cluster.clear();
            }
            //next_point_is_used = next_point.is_used_;
        } else {
          break;
        }
      } while(true); //!next_point_is_used is always true 
  }  

  //std::cout << cluster_counter << std::endl;
  return clusters_vector; 
}


void expand_cluster(point const& current_point, std::vector<point*>& neighbourhood, std::vector<point>& all_points_per_layer, clusters_t& current_cluster, float eps, uint8_t min_points) {
  current_cluster.push_back(current_point);
  //for(auto& neighbour_point : neighbourhood) {
  auto size = neighbourhood.size();

  for(uint32_t point_itr = 0; point_itr < size; ++point_itr) {
    auto& neighbour_point = neighbourhood.at(point_itr);  
    if(!neighbour_point->is_used_){ //if P' is not visited
      neighbour_point->is_used_ = true; //mark P' as visited
      auto next_neighbourhood = find_near_neighbours(eps, *neighbour_point, all_points_per_layer); //region Query (P', eps)
      if(next_neighbourhood.size() >= min_points){ 
        std::cout << "Size of neighbourhood vec 1: " << size << std::endl;
        for(auto& next_neighbour_point : next_neighbourhood){
         // std::cout << "Point id: " <<  next_neighbour_point->id_ << " next_neighbourhood.size(): " << next_neighbourhood.size() << std::endl;
          if(!next_neighbour_point->is_used_) {
            if( std::find(neighbourhood.begin(), neighbourhood.end(), next_neighbour_point) == neighbourhood.end() ){
              std::cout << "Adding a point!" <<  std::endl;
              neighbourhood.push_back(next_neighbour_point);
              //std::cout 
            }
          }
        }

      }
      std::cout << "Size of neighbourhood vec 2: " << size << std::endl;
      size = neighbourhood.size();

      if( !neighbour_point->is_member_of_cluster_ ) {
        neighbour_point->is_member_of_cluster_ = true;
        current_cluster.push_back(*neighbour_point);
         std::cout << "Adding a point to cluster!" <<  std::endl;
      }
    }
  }
}

std::vector<clusters_t> create_DBSCAN_clusters (bins_t& all_surfels_per_layer, float eps, uint8_t min_points){
  //create vector of points, containg all surfels per layer in point format
  std::vector<point> all_points_per_layer;
  uint32_t id_counter = 0; 
  for (auto& s: all_surfels_per_layer){
    point current_point(s.pos_coordinates, id_counter, s.r_, s.g_, s.b_);
    all_points_per_layer.push_back(current_point);
    ++id_counter;
  }

  std::vector<clusters_t> clusters_vec;


  uint32_t current_cluster_id = 0;

  //iterate over all points 
  for (auto& current_point : all_points_per_layer) {
    if(current_point.is_used_){
      continue; //point has been already evaluated
    } 
    else{
      current_point.is_used_ = true;
      auto neighbourhood = find_near_neighbours(eps, current_point, all_points_per_layer);

      neighbourhood.reserve(1000000);
        //std::cout << "Num points in eps range: " << neighbourhood.size() << std::endl;
      if(neighbourhood.size() < min_points){
        continue; //point is considered noise and is not added to any cluster
      }
      else{
        clusters_t current_cluster;
        expand_cluster(current_point, neighbourhood, all_points_per_layer, current_cluster, eps, min_points);

        uint32_t current_cluster_color_id = id_to_color_hash(current_cluster_id);
        lamure::vec3b current_cluster_color = color_array[current_cluster_color_id];
        for( auto& point_in_current_cluster : current_cluster) {
          point_in_current_cluster.set_color(current_cluster_color);
        }
        ++current_cluster_id;

        clusters_vec.push_back(current_cluster);

      }
    }
  }

  return clusters_vec;
}

//using namespace gpucast::math;
std::vector<line> generate_lines_from_curve (std::vector<point>& cluster_of_points) {

  std::cout << "incoming cluster size: " << cluster_of_points.size() << std::endl;
 
  //coppy cluster content to vector of control points

  std::vector<gpucast::math::point3d> control_points_vec(cluster_of_points.size());
  
  #pragma omp parallel for
  for (uint cluster_index = 0; cluster_index < cluster_of_points.size(); ++cluster_index){
    auto x = cluster_of_points.at(cluster_index).pos_coordinates_[0];
    auto y = cluster_of_points.at(cluster_index).pos_coordinates_[1];
    auto z = cluster_of_points.at(cluster_index).pos_coordinates_[2];
    auto weight = 1.0;
    auto current_control_point = gpucast::math::point3d(x, y, z, weight);
    control_points_vec[cluster_index] = current_control_point;
  }

  std::cout << "control_points_vec size: " << control_points_vec.size() << std::endl;

   uint degree = 5;

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
    current_line.length = compute_distance(lamure::vec3f(current_line.start.pos_coordinates_[0], current_line.start.pos_coordinates_[1], current_line.start.pos_coordinates_[2]),
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
    float threshold = 0.01; //TODO think of alternative for dynamic calculation of thershold value

    std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size()); 
    std::vector<line> line_data;
    std::vector<line> line_data_from_sampled_curve;
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
          all_clusters_per_bin_vector = create_clusters(current_bin_of_surfels);
        }else{
          float eps = 0.05;
          uint8_t minPots = 3;
          //std::cout << "XX_C_XX STARTING DB SCAN WITH " << current_bin_of_surfels.size() << " surfels\n";
          all_clusters_per_bin_vector = create_DBSCAN_clusters(current_bin_of_surfels, eps, minPots);
        }

        std::cout << "all_clusters_per_bin_vector.size(): " <<  all_clusters_per_bin_vector.size() << " \n";
        line current_line; 

        /*-std::cout << "XX_C_XX: Num Points per Cluster:\n";
        for(auto const& points_per_cluster : all_clusters_per_bin_vector) {
          std::cout << "XX_C_XX: " << points_per_cluster.size() << "\n";
        } */

        if(all_clusters_per_bin_vector.size() > 0){

          for(auto& current_cluster : all_clusters_per_bin_vector){
            if(current_cluster.size() > 1) { //at least 2 vertices per cluster are need for one complete line 
             // std::cout << "Cluster size " << current_cluster.size() << std::endl;


              lamure::vec3f centroid_pos = compute_cluster_centroid_position(current_cluster);
              auto angle_sorting_lambda = [&](point const& surfel_A,
                                              point const& surfel_B){
                                                   // std::cout << "Sorting, sorting cluster\n";
                                                    lamure::vec3f surfel_position_A (surfel_A.pos_coordinates_[0], surfel_A.pos_coordinates_[1], surfel_A.pos_coordinates_[2]);
                                                    lamure::vec3f surfel_position_B (surfel_B.pos_coordinates_[0], surfel_B.pos_coordinates_[1], surfel_B.pos_coordinates_[2]);
                                                    lamure::vec3f centroid_surfel_vec_A = normalize(surfel_position_A - centroid_pos);
                                                    lamure::vec3f centroid_surfel_vec_B = normalize(surfel_position_B - centroid_pos);
                                                    auto plane_rotation_angle_A = dot(centroid_surfel_vec_A, direction_ref_vector);
                                                    auto plane_rotation_angle_B = dot(centroid_surfel_vec_B, direction_ref_vector);
                                                    //both vectors are in the lower half of the unit cirle => 
                                                    //sort in descending order of angle b/n centroid_serfel_vector and reference vector
                                                    if((surfel_position_A.z <= centroid_pos.z) && (surfel_position_B.z <= centroid_pos.z)){
                                                      return plane_rotation_angle_A >= plane_rotation_angle_B; 
                                                    }
                                                    //both vectors are in the upper half of the unit cirle => 
                                                    //sort in ascending order of angle b/n centroid_serfel_vector and reference vector
                                                    else if((surfel_position_A.z >= centroid_pos.z) && (surfel_position_B.z >= centroid_pos.z)) {

                                                      return plane_rotation_angle_A <= plane_rotation_angle_B;
                                                    } 
                                                    //vectors are in opposite halfs of the unit cirle => 
                                                    //sort in descending order of z coordinate
                                                    else {
                                                      return surfel_position_A.z <= surfel_position_B.z;
                                                    }
                                              };

              std::sort(current_cluster.begin(), current_cluster.end(), angle_sorting_lambda);
              
              if(!use_nurbs){ 
                for (uint j = 0; j < (current_cluster.size()) - 1; ++j){ 
                  current_line.start = current_cluster.at(j);
                  current_line.end = current_cluster.at(j+1);
                  current_line.length = compute_distance(lamure::vec3f(current_line.start.pos_coordinates_[0], current_line.start.pos_coordinates_[1], current_line.start.pos_coordinates_[2]),
                                                          lamure::vec3f(current_line.end.pos_coordinates_[0], current_line.end.pos_coordinates_[1], current_line.end.pos_coordinates_[2]));
                  line_data.push_back(current_line);
                }
              }else {
                line_data_from_sampled_curve = generate_lines_from_curve(current_cluster);
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
           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[2] << "\n";
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
           /*output_file << 80 << " ";
           output_file << 190 << " ";
           output_file << 230 << " ";*/
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
           /*output_file << 30 << " ";
           output_file << 150 << " ";
           output_file << 100 << " ";
           */
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



