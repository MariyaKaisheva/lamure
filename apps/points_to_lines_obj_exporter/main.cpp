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
#include <algorithm>
#include <memory>
#include <limits>


#include <lamure/ren/model_database.h>
#include <lamure/bounding_box.h>
#include <lamure/types.h>
#include <lamure/ren/dataset.h>
#include <lamure/ren/bvh.h>
#include <lamure/ren/lod_stream.h>


#define DEFAULT_PRECISION 15

/*char* get_cmd_option(char** begin, char** end, const std::string & option) {
    char** it = std::find(begin, end, option);
    if (it != end && ++it != end)
        return *it;
    return 0;
}*/

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
point() : pos_coordinates_{0.0, 0.0, 0.0}, id_(0), neighbour_id_(-1), is_used_(false) {}
point(float* pos, int32_t id, int32_t neighbour_id,  bool is_used ) : 
                          pos_coordinates_{pos[0], pos[1], pos[2]},
                          id_(id),
                          neighbour_id_(neighbour_id),
                          is_used_(is_used) {}
float pos_coordinates_[3];
int32_t id_;
int32_t neighbour_id_;
bool is_used_;
};

//std::vector<point> my_vec(10000, point({0.0,0.0,0.0}, -1, false));
using bins_t = std::vector<xyzall_surfel_t>;
using clusters_t = std::vector<point>;
struct line{
    point start;
    point end;
};

//sort in descending order based on y coordinate value 
bool comparator (const xyzall_surfel_t& A, const xyzall_surfel_t& B) {
    return A.pos_coordinates[1] > B.pos_coordinates[1];
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
  std::cout << "centroid_pos: " << average_position << std::endl;
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

std::vector<clusters_t> create_clusters (bins_t& all_surfels_per_layer){
   
  //create vector of points, containg all surfels per layer in point format
  std::vector<point> all_points_per_layer;
  int32_t id_counter = 0; 
  for (auto& s: all_surfels_per_layer){
    point current_point(s.pos_coordinates, id_counter, -1, false);
    all_points_per_layer.push_back(current_point);
    ++ id_counter;
  }

  //auto centroid = compute_average_position_per_layer(all_surfels_per_layer);
  std::vector<point> point_cluster;
  std::vector<clusters_t> clusters_vector;
  float distance_threshold = 8.4; 
  
  float last_measured_distance = 0;  //TODO think about this value

  int cluster_counter  = 0; 

  //distribute points to clusters based on distance to closest available neigbour point
  for(uint point_iterator = 0; point_iterator < all_points_per_layer.size(); ++point_iterator){
    //if(point_iterator >= 0 && point_iterator < all_points_per_layer.size()){
      int point_index = point_iterator;
      bool next_point_is_used = false;
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
              last_measured_distance  = distance_to_nearest_unused_point; 
              current_point->neighbour_id_ = next_point.id_;
              point_index = next_point.id_;
            }
            else{
              //push current cluster to common container and start a new
              //next_point.is_used_ = false; 
              clusters_vector.push_back(point_cluster);
              ++cluster_counter; 
              point_cluster.clear();
            }

            next_point_is_used = next_point.is_used_;
        } else {
          break;
        }
      } while(!next_point_is_used);
/*    }
    else{
      std::cout << "Invalid point_index";
    }
*/
  }  

  std::cout << cluster_counter << std::endl;
  return clusters_vector; 
}

std::vector<line> generate_lines(std::vector<xyzall_surfel_t>& input_data, unsigned long number_line_loops){
    
    //sort input points according to their y-coordinate 
    std::sort(input_data.begin(), input_data.end(), comparator);
    lamure::vec3f direction_ref_vector (1.0, 0.0, 0.0);
    float threshold = 0.002; //TODO think of alternative for dynamic calculation of thershold value

    std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size()); 
    std::vector<line> line_data;
    unsigned long counter = 0; 
    unsigned long offset = floor(input_data.size() / number_line_loops);
    for(uint i = 0; i < number_line_loops; ++i) {
        float current_y_min = input_data.at(counter).pos_coordinates[1];
        float current_y_max = input_data.at(counter + offset).pos_coordinates[1];
        counter += offset;

        float current_y_mean = (current_y_min + current_y_max) / 2.0;

        auto copy_lambda = [&]( xyzall_surfel_t const& surfel){return (surfel.pos_coordinates[1] >= current_y_mean - threshold) && (surfel.pos_coordinates[1] <= current_y_mean + threshold);};
        auto it = std::copy_if(input_data.begin(), input_data.end(), current_bin_of_surfels.begin(), copy_lambda);
        current_bin_of_surfels.resize(std::distance(current_bin_of_surfels.begin(), it));
        
        for (auto& surfel : current_bin_of_surfels) {
          surfel.pos_coordinates[1] = current_y_mean; //project all surfels for a given y_mean value to a single plane
        }
        //std::sort(current_bin_of_surfels.begin(), current_bin_of_surfels.end(), comparator_x_dim);
        auto all_clsters_per_bin_vector = create_clusters(current_bin_of_surfels);

        std::cout << "all_clsters_per_bin_vector.size(): " <<  all_clsters_per_bin_vector.size() << " \n";
        line current_line; 
        if(all_clsters_per_bin_vector.size() > 0){

          for(auto& current_cluster : all_clsters_per_bin_vector){
            if(current_cluster.size() > 1) {
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

              //std::sort(current_cluster.begin(), current_cluster.end(), angle_sorting_lambda);
              for (uint j = 0; j < (current_cluster.size()) - 1; ++j){ 
                current_line.start = current_cluster.at(j);
                current_line.end = current_cluster.at(j+1);
                line_data.push_back(current_line);
              }
            }
          }
        }
        else{
          std::cout << "no clusters in the current layer \n"; 
        }

    }   
    std::cout << "num Limes (should be >= num layers??) " << line_data.size() << std::endl; 
    return line_data;   
}
int main(int argc, char *argv[]) {
    
    if (argc == 1 ) {
        
      std::cout << "Usage: " << argv[0] << "<flags> -f <input_file>" << std::endl <<
         "\t-f: selects .bvh input file" << std::endl <<
         std::endl;
      return 0;
    }

   // std::string bvh_filename = std::string(get_cmd_option(argv, argv + argc, "-f"));
    std::string bvh_filename = argv[1];
    std::string ext = bvh_filename.substr(bvh_filename.size()-3);
    if (ext.compare("bvh") != 0) {
        std::cout << "please specify a .bvh file as input" << std::endl;
        return 0;
    }

    lamure::ren::bvh* bvh = new lamure::ren::bvh(bvh_filename);
    int32_t depth = -1;
    if(argc > 2){
      depth = atoi(argv[2]);
      if(depth > int(bvh->get_depth()) || depth < 0){
        depth = bvh->get_depth();
      }
    }
    else{
      depth = bvh->get_depth();
    }
    unsigned long number_line_loops = 20; //TODO make number dependent on model size??
    if(argc > 3){
     number_line_loops = atoi(argv[3]); //user input
    }
    
    std::string obj_filename = bvh_filename.substr(0, bvh_filename.size()-4)+ "_d" + std::to_string(depth) + "_l" + std::to_string(number_line_loops) + ".obj";
    std::cout << "input: " << bvh_filename << std::endl;
    std::cout << "output: " << obj_filename << std::endl;
   
    std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
    lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
    in_access->open(lod_filename);

    size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);

    std::cout << "working with surfels at depth " << depth << std::endl;
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

    auto line_data = generate_lines(surfels_vector, number_line_loops);

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
   
    std::cout << "ok \n";

    delete in_access;
    delete bvh;

    return 0;
}



