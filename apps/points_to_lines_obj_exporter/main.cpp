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
#include <cmath>


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
  float x_, y_, z_;
  uint8_t r_, g_, b_, fake_;
  float radius_;
  float nx_, ny_, nz_;
};

struct line{
    xyzall_surfel_t start;
    xyzall_surfel_t end;
};

//sort in descending order based on y coordinate value 
bool comparator (const xyzall_surfel_t& A, const xyzall_surfel_t& B) {
    return A.y_ > B.y_;
};

lamure::vec3f compute_average_position_per_layer (std::vector<xyzall_surfel_t> const& surfels_vector) {
  float average_x = 0.0;
  float average_y = 0.0;
  float average_z = 0.0;
  for(auto const& surfel : surfels_vector){
    average_x += surfel.x_;
    average_y += surfel.y_;
    average_z += surfel.z_;
  }
  auto number_of_surfels_per_layer = surfels_vector.size();
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

std::vector<line> generate_lines(std::vector<xyzall_surfel_t>& input_data){
    
    //sort input points according to their y-coordinate 
    std::sort(input_data.begin(), input_data.end(), comparator);
    lamure::vec3f direction_ref_vector (1.0, 0.0, 0.0);
    unsigned long number_line_loops = 20; //TODO make number dependent on model size
    float threshold = 0.002; //TODO think of alternative for dynamic calculation of thershold value
    std::vector< std::vector <xyzall_surfel_t>> bins; 
    std::vector<xyzall_surfel_t> current_bin_of_surfels(input_data.size()); 
    std::vector<line> line_data;
    unsigned long counter = 0; 
    unsigned long offset = floor(input_data.size() / number_line_loops);
    for(uint i = 0; i < number_line_loops; ++i) {
        float current_y_min = input_data.at(counter).y_;
        float current_y_max = input_data.at(counter + offset).y_;
        counter += offset;

        float current_y_mean = (current_y_min + current_y_max) / 2.0;

        auto copy_lambda = [&]( xyzall_surfel_t const& surfel){return (surfel.y_ >= current_y_mean - threshold) && (surfel.y_ <= current_y_mean + threshold);};
        auto it = std::copy_if(input_data.begin(), input_data.end(), current_bin_of_surfels.begin(), copy_lambda);
        current_bin_of_surfels.resize(std::distance(current_bin_of_surfels.begin(), it));
        
    
        for (auto& surfel : current_bin_of_surfels) {
          surfel.y_ = current_y_mean; //project all surfels for a given y_mean value to a single plane
        }

        lamure::vec3f centroid = compute_average_position_per_layer(current_bin_of_surfels);

        auto angle_sorting_lambda = [&](xyzall_surfel_t const& surfel_A,
                                        xyzall_surfel_t const& surfel_B){

                                              lamure::vec3f surfel_position_A (surfel_A.x_, surfel_A.y_, surfel_A.z_);
                                              lamure::vec3f surfel_position_B (surfel_B.x_, surfel_B.y_, surfel_B.z_);
                                              lamure::vec3f centroid_surfel_vec_A = normalize(surfel_position_A - centroid);
                                              lamure::vec3f centroid_surfel_vec_B = normalize(surfel_position_B - centroid);
                                              auto plane_rotation_angle_A = dot(centroid_surfel_vec_A, direction_ref_vector);
                                              auto plane_rotation_angle_B = dot(centroid_surfel_vec_B, direction_ref_vector);
                                              //both vectors are in the lower half of the unit cirle => 
                                              //sort in descending order of angle b/n centroid_serfel_vector and reference vector
                                              if((surfel_position_A.z <= centroid.z) && (surfel_position_B.z <= centroid.z)){
                                                return plane_rotation_angle_A >= plane_rotation_angle_B; 
                                              }
                                              //both vectors are in the upper half of the unit cirle => 
                                              //sort in ascending order of angle b/n centroid_serfel_vector and reference vector
                                              else if((surfel_position_A.z >= centroid.z) && (surfel_position_B.z >= centroid.z)) {

                                                return plane_rotation_angle_A <= plane_rotation_angle_B;
                                              } 
                                              //vectors are in opposite halfs of the unit cirle => 
                                              //sort in descending order of z coordinate
                                              else {
                                                return surfel_position_A.z <= surfel_position_B.z;
                                              }
                                        };

        std::sort(current_bin_of_surfels.begin(), current_bin_of_surfels.end(), angle_sorting_lambda);
        line current_line; 
        if(current_bin_of_surfels.size() > 2) {
          for (uint j = 0; j < (current_bin_of_surfels.size()) - 1; ++j){ 
            current_line.start = current_bin_of_surfels.at(j);
            current_line.end = current_bin_of_surfels.at(j+1);
            line_data.push_back(current_line);
          }
        }
    }   
    
    return line_data;   
};

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

    std::string obj_filename = bvh_filename.substr(0, bvh_filename.size()-4)+ "_l" + std::to_string(depth) + ".obj";
    std::cout << "input: " << bvh_filename << std::endl;
    std::cout << "output: " << obj_filename << std::endl;
   
    std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
    lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
    in_access->open(lod_filename);

    size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);

    std::cout << "working with surfels at depth " << depth << std::endl;
    lamure::node_t num_leafs = bvh->get_length_of_depth(depth);

    std::ofstream out_stream;
    out_stream.open(obj_filename, std::ios::out | std::ios::trunc); // Any current content is discarded
    out_stream.close();

    //uint64_t num_surfels_excluded = 0;

    auto total_num_surfels = num_leafs*bvh->get_primitives_per_node();
    std::vector<xyzall_surfel_t> surfels_vector(total_num_surfels);

    auto num_offset_nodes = std::pow(2, depth) - 1; //only valid for fanfactor 2
    auto length_in_bytes = num_leafs * size_of_node;
    in_access->read((char*) &surfels_vector[0], num_offset_nodes * size_of_node, length_in_bytes);
  
    //std::ios::openmode mode = std::ios::out | std::ios::app;
    //out_stream.open(obj_filename, mode);

    //std::string filestr;
    //std::stringstream ss(filestr);
    
    //clean data from degenerated padding surfels
    surfels_vector.erase(std::remove_if(surfels_vector.begin(), 
                                        surfels_vector.end(),
                                        [](xyzall_surfel_t s){return s.radius_ <= 0.0;}),
                         surfels_vector.end());

    /*for (unsigned int i = 0; i < surfels_vector.size() ; ++i) {

        xyzall_surfel_t const& s = surfels_vector.at(i); 
        
        if(s.radius_ <= 0.0f){
          ++num_surfels_excluded;
          continue;
        } 
        ss << "v " <<  s.x_ << " " << s.y_ << " " << s.z_ << "\n"; //TODO << std::setprecision(DEFAULT_PRECISION)
    }*/

    auto line_data = generate_lines(surfels_vector);

    std::ofstream output_file(obj_filename);
 
    unsigned long vert_counter = 1;
    //uint i = 0;

    if (output_file.is_open()){
        std::cout << "ok \n";
        for (uint i = 0; i < line_data.size(); ++i){
        //while(i < line_data.size()){
         output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).start.x_ << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).start.y_ << " " << std::setprecision(DEFAULT_PRECISION)<< line_data.at(i).start.z_ << "\n";
         output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.x_ << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.y_ << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.z_ << "\n";
         //vertex duplication to emulate triangle
         output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.x_ << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.y_ << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.z_ << "\n";
         output_file << "f " << vert_counter << " " << (vert_counter + 1) << " " << (vert_counter + 2) << "\n";
         vert_counter += 3;
         //i += 250;
        }
        output_file.close();
    }
    else{
      std::cout << "Cannot open output file to write to! \n";
    }
    //out_stream << std::setprecision(DEFAULT_PRECISION) << ss.rdbuf();
    //out_stream.close();
    //std::cout << "done. (" << num_surfels_excluded << " surfels out of " << surfels_vector.size()<< " excluded)" << std::endl;

    delete in_access;
    delete bvh;

    return 0;
}


