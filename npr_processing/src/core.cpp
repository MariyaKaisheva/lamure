#include <lamure/npr/core.h>

#include <chrono>

namespace npr {
namespace core {

void
 generate_line_art(scm::math::mat4f const& user_defined_rot_mat, 
                   std::string const& bvh_filename,
                   int32_t depth, 
                   bool write_intermediate_results,
                   bool spiral_look,
                   std::string output_base_name,
                   float min_distance,
                   float max_distance,
                   bool use_nurbs,
                   bool apply_alpha_shapes,
                   bool without_lod_adjustment,
                   bool is_verbose)
    {

        lamure::ren::bvh* bvh = new lamure::ren::bvh(bvh_filename);

        std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
        lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
        in_access->open(lod_filename);

        size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);

        if(depth > int(bvh->get_depth()) || depth < 0){
        	depth = bvh->get_depth();
      	}

        if(is_verbose) {
            std::cout << "------------------------------------------------------------------------------\n" << std::endl;
            std::cout << "------------------------------------- L O G ----------------------------------\n" << std::endl;
            std::cout << "------------------------------------------------------------------------------\n" << std::endl;

            std::cout << "Working with surfels at depth " << depth << "; Max depth for this model is "<<  bvh->get_depth() 
                      << std::endl <<std::endl;
        }

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
            auto& surfel = surfels_vector[surfel_idx];
            scm::math::vec3f original_coord = scm::math::vec4f(surfel.pos_coordinates[0], surfel.pos_coordinates[1], surfel.pos_coordinates[2], 1.0f);
            auto transformed_coord = inverse_rot_mat * original_coord;
            surfel.pos_coordinates[0] = transformed_coord.x;
            surfel.pos_coordinates[1] = transformed_coord.y;
            surfel.pos_coordinates[2] = transformed_coord.z; 
        }

        io::stage_content_storage intermediate_visalization_struct;


        //sort in descending order based on y coordinate value 
        auto comparator = [] (const xyzall_surfel_t& A, const xyzall_surfel_t& B) {
            return A.pos_coordinates[1] < B.pos_coordinates[1];
        };

        //sort input points according to their y-coordinate 
        std::sort(surfels_vector.begin(), surfels_vector.end(), comparator);

        
        //validate input 
        if(min_distance < 0 || max_distance < 0){
            auto const& suggested_distance_thresholds = utils::estimate_binning_densities(surfels_vector, is_verbose);
            min_distance = suggested_distance_thresholds.first;
            max_distance = suggested_distance_thresholds.second;
            if(is_verbose){
                std::cout << "\t min_distance: " << min_distance << "\n";
                std::cout << "\t max_distance: " << max_distance << "\n";
                std::cout << "----------------------------\n";
            }
        }


        auto line_data =  line_gen::generate_lines(surfels_vector,
                                                   min_distance, max_distance,
                                                   output_base_name, //used to write out intermediate stages
                                                   write_intermediate_results,
                                                   use_nurbs, apply_alpha_shapes,
                                                   spiral_look, is_verbose);

        if(is_verbose) {
            std::cout << "Num generated lines: " << line_data.size() << "\n";
        }

        //transform data again to return to the original model orientation
        std::chrono::time_point<std::chrono::system_clock> start_inverse_rotation, end_inverse_rotation;
        start_inverse_rotation = std::chrono::system_clock::now();
        utils::transform(line_data, user_defined_rot_mat);
        end_inverse_rotation = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_inverse_rotation = end_inverse_rotation - start_inverse_rotation;


        std::string obj_filename = output_base_name + ".obj";
        std::string xyz_all_filename = output_base_name + ".xyz_all";


        std::chrono::time_point<std::chrono::system_clock> start_writing_output, end_writing_output;
        start_writing_output = std::chrono::system_clock::now();

        io::write_output( obj_filename, line_data, bvh);

        end_writing_output = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_writing_output = end_writing_output - start_writing_output;
        
        if(is_verbose) {

            std::cout << "\t --  Time LOG:  -- rotating model to original position: " << elapsed_seconds_inverse_rotation.count() << "s\n";
            std::cout << "\t --  Time LOG:  -- writing output: " << elapsed_seconds_writing_output.count() << "s\n";
            std::cout << "-----------------------------------\n" << std::endl;
            std::cout << "NURBS usage: " <<  use_nurbs << std::endl;
            std::cout << "Alpha-shapes usage: " <<  apply_alpha_shapes << std::endl;
            std::cout << "--------------- ok ----------------\n";
        }

        delete in_access;
        delete bvh;
    }

} //namespace core
} //namespace npr