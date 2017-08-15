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

    core::generate_line_art(user_defined_rot_mat, bvh_filename, depth, 
                        write_obj_file, use_nurbs, apply_alpha_shapes,
                        max_number_line_loops, is_verbose);

    return 0;
}