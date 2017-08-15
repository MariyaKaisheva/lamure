#include "core.h"


void core::
 generate_line_art(scm::math::mat4f user_defined_rot_mat, 
                std::string bvh_filename,
                int32_t depth, 
                bool write_obj_file,
                bool use_nurbs,
                bool apply_alpha_shapes,
                uint32_t max_number_line_loops)
    {

        lamure::ren::bvh* bvh = new lamure::ren::bvh(bvh_filename);

        std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
        lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
        in_access->open(lod_filename);

        size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);

        if(depth > int(bvh->get_depth()) || depth < 0){
        	depth = bvh->get_depth();
      	}

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

        size_t num_surfels_in_vector = surfels_vector.size();
        #pragma omp parallel for
        for(size_t surfel_idx = 0; surfel_idx < num_surfels_in_vector; ++surfel_idx) {
        //for(auto& surfel : surfels_vector){
            auto& surfel = surfels_vector[surfel_idx];
            scm::math::vec3f original_coord = scm::math::vec4f(surfel.pos_coordinates[0], surfel.pos_coordinates[1], surfel.pos_coordinates[2], 1.0f);
            auto transformed_coord = inverse_rot_mat * original_coord;
            surfel.pos_coordinates[0] = transformed_coord.x;
            surfel.pos_coordinates[1] = transformed_coord.y;
            surfel.pos_coordinates[2] = transformed_coord.z; 
        }

        //create line representation of original input data
        auto line_data = line_gen::generate_lines(surfels_vector, max_number_line_loops, use_nurbs, apply_alpha_shapes);

        std::cout << "Num generated lines: " << line_data.size() << "\n";
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

        std::string bvh_filename_without_path = bvh_filename.substr(bvh_filename.find_last_of("/\\") + 1); 
        std::string bvh_filename_without_path_and_extension = bvh_filename_without_path.substr(0, bvh_filename_without_path.size() - 4 );

        //std::string rot_substring_without_path = rot_filename.substr(rot_filename.find_last_of("/\\") + 1); 
        //std::string rot_substring_without_path_and_extension = rot_substring_without_path.substr(0, rot_substring_without_path.size() - 4 );

        //std::string rot_substring = rot_filename.substr(0, rot_filename.size()-4);

        std::string output_files_base_name =   bvh_filename_without_path_and_extension + "_LA_preview";
                                            /* + "_d" + std::to_string(depth) 
                                             + "_l" + std::to_string(max_number_line_loops) 
                                             + "_"  + /*+ rot_substring_without_path_and_extension;*/

        std::string obj_filename = output_files_base_name + ".obj";
        std::string xyz_all_filename = output_files_base_name + ".xyz_all";
       

        if(write_obj_file){
          io::write_output(write_obj_file, obj_filename, line_data, bvh);
        }else{
          io::write_output(write_obj_file, xyz_all_filename, line_data, bvh);
        }
         
        std::cout << "NURBS usage: " <<  use_nurbs << std::endl;
        std::cout << "Alpha-shapes usage: " <<  apply_alpha_shapes << std::endl;  
        std::cout << "--------------- ok ----------------\n";

        delete in_access;
        delete bvh;
    }  