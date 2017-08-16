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
    if (argc == 1 || io::cmd_option_exists(argv, argv+argc, "-h") ||
        !io::cmd_option_exists(argv, argv+argc, "-f") || 
        !io::cmd_option_exists(argv, argv+argc, "-t")) {
        io::print_help_message(argv);
      return 0;
    }

    std::string rot_filename = std::string(io::get_cmd_option(argv, argv + argc, "-t"));
    std::string bvh_filename = std::string(io::get_cmd_option(argv, argv + argc, "-f"));
    std::string ext_bvh = bvh_filename.substr(bvh_filename.size()-3);
    std::string ext_rot = rot_filename.substr(rot_filename.size()-3);
    if ((ext_bvh.compare("bvh") != 0) && (ext_rot.compare("rot") != 0)){
        std::cout << "please specify a .bvh file as input" << std::endl;
        return 0;
    }

    int32_t depth = -1;
    if(io::cmd_option_exists(argv, argv+argc, "-d")){
      depth = atoi(io::get_cmd_option(argv, argv+argc, "-d"));
    }

    bool write_obj_file = !io::cmd_option_exists(argv, argv + argc, "--write_xyz_points");
    uint32_t max_number_line_loops = 85;

    bool is_verbose_option_1 = io::cmd_option_exists(argv, argv + argc, "--verbose");
    bool is_verbose_option_2 = io::cmd_option_exists(argv, argv + argc, "-v");
    bool is_verbose = is_verbose_option_1 | is_verbose_option_2;

    if(io::cmd_option_exists(argv, argv+argc, "-l")){
     max_number_line_loops = atoi(io::get_cmd_option(argv, argv+argc, "-l")); //user input
    }
    //check user preferences
    bool use_nurbs = io::cmd_option_exists(argv, argv + argc, "--apply_nurbs_fitting");
    bool apply_alpha_shapes = io::cmd_option_exists(argv, argv + argc, "--apply_alpha_shapes");
    auto user_defined_rot_mat = io::read_in_transformation_file(rot_filename);


    std::string bvh_filename_without_path = bvh_filename.substr(bvh_filename.find_last_of("/\\") + 1); 
    std::string bvh_filename_without_path_and_extension = bvh_filename_without_path.substr(0, bvh_filename_without_path.size() - 4 );

    //std::string rot_substring_without_path = rot_filename.substr(rot_filename.find_last_of("/\\") + 1); 
    //std::string rot_substring_without_path_and_extension = rot_substring_without_path.substr(0, rot_substring_without_path.size() - 4 );

    //std::string rot_substring = rot_filename.substr(0, rot_filename.size()-4);
    scm::math::quat<double> output_quat;
    output_quat = scm::math::quat<double>::from_matrix(scm::math::mat4d(user_defined_rot_mat));
    double angle; 
    scm::math::vec3d axis; 
    output_quat.retrieve_axis_angle(angle, axis);

    std::string output_base_name = bvh_filename_without_path_and_extension
                                         + "_d" + std::to_string(depth)
                                         + "_angle_"  + std::to_string(angle);

    core::generate_line_art(user_defined_rot_mat, bvh_filename, depth, 
                        write_obj_file, use_nurbs, apply_alpha_shapes,
                        output_base_name, max_number_line_loops, is_verbose);

    return 0;
}