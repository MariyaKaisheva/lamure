// Copyright (c) 2014 Bauhaus-Universitaet Weimar
// This Software is distributed under the Modified BSD License, see license.txt.
//
// Virtual Reality and Visualization Research Group 
// Faculty of Media, Bauhaus-Universitaet Weimar
// http://www.uni-weimar.de/medien/vr
#include <lamure/npr/core.h>

using namespace npr;

int main(int argc, char** argv)
{
    //input validation
    bool valid_input = io::check_user_input(argv, argc);
    if (argc == 1 || io::cmd_option_exists(argv, argv+argc, "-h") ||
        !io::cmd_option_exists(argv, argv+argc, "-f") || !valid_input) {

        io::print_help_message(argv);
        return 0;
    }

    //get input for point cloud model and corresponding rotaion file
    std::string bvh_filename = std::string(io::get_cmd_option(argv, argv + argc, "-f"));
    std::string ext_bvh = bvh_filename.substr(bvh_filename.size()-3);
    if ((ext_bvh.compare("bvh") != 0) ){ //validate input paramenter
        std::cout << "please specify a .bvh file as input" << std::endl;
        return 0;
    }

    //get input for rotaion file (it must contain 4 values that describe the ange and axis of rotation of slicing plane)
    //if no such input provided, use defaut plane orientation 
    scm::math::mat4f user_defined_rot_mat = scm::math::make_rotation(0.0f, 0.0f, 0.0f, 1.0f);
    std::string rot_filename = "";
    if(io::cmd_option_exists(argv, argv+argc, "-t")){
        rot_filename = std::string(io::get_cmd_option(argv, argv + argc, "-t"));
        std::string ext_rot = rot_filename.substr(rot_filename.size()-3);
        if(ext_rot.compare("rot") != 0){ //validate input paramenter
            std::cout << "Please specify correct .rot file after -t " << std::endl;
            return 0;
        }
        user_defined_rot_mat = io::read_in_transformation_file(rot_filename);
    }


    //check user preferences
    int32_t depth = -1;
    if(io::cmd_option_exists(argv, argv+argc, "-d")){
      depth = atoi(io::get_cmd_option(argv, argv+argc, "-d"));
    }

    bool write_obj_file = !io::cmd_option_exists(argv, argv + argc, "--write_xyz_points");
    uint32_t max_number_line_loops = 95;

    bool is_verbose_option_1 = io::cmd_option_exists(argv, argv + argc, "--verbose");
    bool is_verbose_option_2 = io::cmd_option_exists(argv, argv + argc, "-v");
    bool is_verbose = is_verbose_option_1 | is_verbose_option_2;

    bool without_lod_adjustment = io::cmd_option_exists(argv, argv + argc, "--no_reduction");

    if(io::cmd_option_exists(argv, argv+argc, "-l")){
     max_number_line_loops = atoi(io::get_cmd_option(argv, argv+argc, "-l")); //user input
    }

    bool use_nurbs = !io::cmd_option_exists(argv, argv + argc, "--no_nurbs_fitting");
    bool apply_alpha_shapes = !io::cmd_option_exists(argv, argv + argc, "--no_alpha_shapes");
    bool spiral_look = io::cmd_option_exists(argv, argv + argc, "--generate_spirals");
    
    std::string bvh_filename_without_path = bvh_filename.substr(bvh_filename.find_last_of("/\\") + 1); 
    std::string bvh_filename_without_path_and_extension = bvh_filename_without_path.substr(0, bvh_filename_without_path.size() - 4 );

    scm::math::quat<double> output_quat;
    output_quat = scm::math::quat<double>::from_matrix(scm::math::mat4d(user_defined_rot_mat));
    double angle; 
    scm::math::vec3d axis; 
    output_quat.retrieve_axis_angle(angle, axis);

    std::string output_base_name = bvh_filename_without_path_and_extension
                                         + "_d" + std::to_string(depth)
                                         + "_angle_"  + std::to_string(angle);

    core::generate_line_art(user_defined_rot_mat, bvh_filename, depth, 
                            write_obj_file, use_nurbs,
                            apply_alpha_shapes, spiral_look,
                            output_base_name, max_number_line_loops, 
                            without_lod_adjustment, is_verbose);

    return 0;
}