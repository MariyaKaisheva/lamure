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

#include <scm/gl_core/math.h>

#include <lamure/ren/model_database.h>
#include <lamure/bounding_box.h>
#include <lamure/types.h>
#include <lamure/ren/dataset.h>
#include <lamure/ren/bvh.h>
#include <lamure/ren/lod_stream.h>

#include "input_output.h"
#include "binning.h"
#include "utils.h"
#include "sampling.h"
#include "color_hash_map.hpp"
#include "math_wrapper.h"
#include "point.hpp"
#include "alpha-shapes_wrapper.h"
#include "line_gen.h"


int main(int argc, char** argv)
{
  
    if (argc == 1 || io::cmd_option_exists(argv, argv+argc, "-h") ||
        !io::cmd_option_exists(argv, argv+argc, "-f") || 
        !io::cmd_option_exists(argv, argv+argc, "-t")) {
        io::print_help_message(argv);
      return 0;
    }

    std::string txt_filename = std::string(io::get_cmd_option(argv, argv + argc, "-t"));
    auto user_defined_rot_mat = io::read_in_transformation_file(txt_filename);

    std::string bvh_filename = std::string(io::get_cmd_option(argv, argv + argc, "-f"));
    std::string ext = bvh_filename.substr(bvh_filename.size()-3);
    if (ext.compare("bvh") != 0) {
        std::cout << "please specify a .bvh file as input" << std::endl;
        return 0;
    }

    lamure::ren::bvh* bvh = new lamure::ren::bvh(bvh_filename);
    int32_t depth = -1;
    if(io::cmd_option_exists(argv, argv+argc, "-d")){
      depth = atoi(io::get_cmd_option(argv, argv+argc, "-d"));
      if(depth > int(bvh->get_depth()) || depth < 0){
        depth = bvh->get_depth();
      }
    }
    else{
      depth = bvh->get_depth();
    }

    bool write_obj_file = !io::cmd_option_exists(argv, argv + argc, "--write_xyz_points");
    uint32_t max_number_line_loops = 85;
    if(io::cmd_option_exists(argv, argv+argc, "-l")){
     max_number_line_loops = atoi(io::get_cmd_option(argv, argv+argc, "-l")); //user input
    }
    

    std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
    lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
    in_access->open(lod_filename);

    size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);

    std::cout << "working with surfels at depth " << depth << "; Max depth for this model is "<<  bvh->get_depth() <<std::endl;
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

    //apply inverser transformation on the local coordinate system of the model
    //such that y_axis will allign with user-defined direction 
    auto inverse_rot_mat = scm::math::inverse(user_defined_rot_mat);
    for(auto& surfel : surfels_vector){
        scm::math::vec3f original_coord = scm::math::vec4f(surfel.pos_coordinates[0], surfel.pos_coordinates[1], surfel.pos_coordinates[2], 1.0f);
        auto transformed_coord = inverse_rot_mat * original_coord;
        surfel.pos_coordinates[0] = transformed_coord.x;
        surfel.pos_coordinates[1] = transformed_coord.y;
        surfel.pos_coordinates[2] = transformed_coord.z; 
    }

    //check user preferences
    bool use_nurbs = io::cmd_option_exists(argv, argv + argc, "--apply_nurbs_fitting");
   //bool apply_naive_clustering = !io::cmd_option_exists(argv, argv + argc, "--use_dbscan");
    bool apply_alpha_shapes = io::cmd_option_exists(argv, argv + argc, "--apply_alpha_shapes");


    //create line representation of original input data
    auto line_data = generate_lines(surfels_vector, max_number_line_loops, use_nurbs, apply_alpha_shapes);

    //transform data again to return to the original model orientation 
    utils::transform(line_data, user_defined_rot_mat);
    
    #if 0 //remove potential oulier line segments; 
    std::cout << "Num lines BEFORE clean up: " << line_data.size() << std::endl;
    auto avg_line_length = utils::compute_global_average_line_length(line_data); 
    line_data.erase(std::remove_if(line_data.begin(),
                                   line_data.end(),
                                   [&](line l){return l.length >= 10 * avg_line_length;}),
                    line_data.end());
    std::cout << "Num lines AFTER clean up: " << line_data.size() << std::endl;
    #endif

    std::string obj_filename = bvh_filename.substr(0, bvh_filename.size()-4)+ "_d" + std::to_string(depth) + "_l" + std::to_string(max_number_line_loops) + ".obj";
    std::string xyz_all_filename = bvh_filename.substr(0, bvh_filename.size()-4)+ "_d" + std::to_string(depth) + "_l" + std::to_string(max_number_line_loops) + ".xyz_all";
   

    if(write_obj_file){
      io::write_output(write_obj_file, obj_filename, line_data, bvh);
    }else{
      io::write_output(write_obj_file, xyz_all_filename, line_data, bvh);
    }
        
    std::cout << "ok \n";

    delete in_access;
    delete bvh;

    return 0;
}