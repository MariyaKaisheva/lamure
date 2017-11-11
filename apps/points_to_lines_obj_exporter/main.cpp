// Copyright (c) 2014 Bauhaus-Universitaet Weimar
// This Software is distributed under the Modified BSD License, see license.txt.
//
// Virtual Reality and Visualization Research Group 
// Faculty of Media, Bauhaus-Universitaet Weimar
// http://www.uni-weimar.de/medien/vr
#include <lamure/npr/core.h>
#include <random>

using namespace npr;



void generate_random_orthogonal_90_degree_orientations() {

   std::random_device rd;
   std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
   
   std::uniform_real_distribution<> distribution_between_0_and_1(0.0, 1.0);
   
   double rand_x_axis_component = distribution_between_0_and_1(gen);
   double rand_y_axis_component = distribution_between_0_and_1(gen);
   double rand_z_axis_component = distribution_between_0_and_1(gen);

    scm::math::vec3d random_rotation_axis{rand_x_axis_component,
                                          rand_y_axis_component,
                                          rand_z_axis_component};

    //handle degenerated case
    if(random_rotation_axis == scm::math::vec3d(0.0, 0.0, 0.0)) {
        random_rotation_axis = scm::math::vec3d(1.0, 0.0, 0.0); 
    }

    random_rotation_axis = scm::math::normalize(random_rotation_axis);

    // compute arbitrary tangent vectors (see surfel rendering code)
    scm::math::vec3d ms_n = random_rotation_axis;
    scm::math::vec3d tmp_ms_u = scm::math::vec3d(0.0, 0.0, 0.0);
    if(random_rotation_axis[2] != 0.0) {
      tmp_ms_u = scm::math::vec3d( 1.0, 1.0, (-random_rotation_axis[0] -random_rotation_axis[1])/random_rotation_axis[2]);
    } else if (ms_n.y != 0.0) {
      tmp_ms_u = scm::math::vec3d( 1.0, (-random_rotation_axis[0] -random_rotation_axis[2])/random_rotation_axis[1], 1.0);
    } else {
      tmp_ms_u = scm::math::vec3d( (-random_rotation_axis[1] -random_rotation_axis[2])/random_rotation_axis[0], 1.0, 1.0);
    }


  // assign tangent vectors
  scm::math::vec3d ms_u = scm::math::normalize(tmp_ms_u);
  scm::math::vec3d ms_v = scm::math::normalize(scm::math::cross(ms_n, ms_u));

  std::ofstream rand_arbitrary_axis_file_0("rand_perp_axis_0.rot");
  rand_arbitrary_axis_file_0 << "90.0 " << ms_n[0] << " " << ms_n[1] << " " << ms_n[2];
  rand_arbitrary_axis_file_0.close();

  std::ofstream rand_tangent_axis_file_1("rand_perp_axis_1.rot");
  rand_tangent_axis_file_1 << "90.0 " << ms_u[0] << " " << ms_u[1] << " " << ms_u[2];
  rand_tangent_axis_file_1.close();

  std::ofstream rand_tangent_axis_file_2("rand_perp_axis_2.rot");
  rand_tangent_axis_file_2 << "90.0 " << ms_v[0] << " " << ms_v[1] << " " << ms_v[2];
  rand_tangent_axis_file_2.close();
}

int main(int argc, char** argv)
{

    bool generate_random_perpendicular_transformations = io::cmd_option_exists(argv, argv + argc, "--create_random_axes");
    
    if(generate_random_perpendicular_transformations) {
      generate_random_orthogonal_90_degree_orientations();
      std::cout << "\nSuccessfully created three random perpendicular transformation files: [rand_perp_axis_0.rot, rand_perp_axis_1.rot, rand_perp_axis_2.rot]\n";
      std::cout << "Exiting.\n\n";
      return 0;
    }

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
        std::cout << "Please, specify a .bvh file as input" << std::endl;
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
            std::cout << "Please, specify correct .rot file after -t " << std::endl;
            return 0;
        }
        user_defined_rot_mat = io::read_in_transformation_file(rot_filename);
    }


    //check user preferences
    int32_t depth = -1;
    if(io::cmd_option_exists(argv, argv+argc, "-d")){
      depth = atoi(io::get_cmd_option(argv, argv+argc, "-d"));
    }

    //parse intermediate stage enabling
    bool write_intermediate_results = io::cmd_option_exists(argv, argv + argc, "--write_stages");

    //parse verbosity parameters
    bool is_verbose_option_1 = io::cmd_option_exists(argv, argv + argc, "--verbose");
    bool is_verbose_option_2 = io::cmd_option_exists(argv, argv + argc, "-v");
    bool is_verbose = is_verbose_option_1 | is_verbose_option_2;

    //parse slicing mode
    bool use_radial_slicing = io::cmd_option_exists(argv, argv + argc, "-r");

    //parse bin-distance related parameters (default: estimated later)
    float min_distance = -1.0;
    float max_distance = -1.0;
    io::parse_float_parameter(argc, argv, min_distance, "--min");
    io::parse_float_parameter(argc, argv, max_distance, "--max");
    

    //parse DBSCAN eps_factor (default: 10)
    float eps_factor = 10.0;
    io::parse_float_parameter(argc, argv, eps_factor, "--eps");

    //parse line-color related parameters (default: black)
    float red_channel_value = 0.0;
    float green_channel_value = 0.0;
    float blue_channel_value = 0.0;

    io::parse_float_parameter(argc, argv, red_channel_value, "--red");
    io::parse_float_parameter(argc, argv, green_channel_value, "--green");
    io::parse_float_parameter(argc, argv, blue_channel_value, "--blue");

    bool spiral_look = io::cmd_option_exists(argv, argv + argc, "--generate_spirals");
    
    //retrieve rotation axis and angle
    scm::math::quat<double> output_quat;
    output_quat = scm::math::quat<double>::from_matrix(scm::math::mat4d(user_defined_rot_mat));
    double angle; 
    scm::math::vec3d axis; 
    output_quat.retrieve_axis_angle(angle, axis);

    //get model related path names
    std::string bvh_filename_without_path = bvh_filename.substr(bvh_filename.find_last_of("/\\") + 1); 
    std::string bvh_filename_without_path_and_extension = bvh_filename_without_path.substr(0, bvh_filename_without_path.size() - 4 );

    //retrieve the name for the model creation based on the parameters
    std::string output_base_name = npr::io::create_output_base_name(bvh_filename_without_path_and_extension, depth, angle, axis, spiral_look, use_radial_slicing, min_distance, max_distance, eps_factor );

    scm::math::vec3f translation_components_vec(0.0, 0.0, 0.0);

    //call to the actual line art creation npr-library function
    core::generate_line_art(user_defined_rot_mat,
                            translation_components_vec,
                            bvh_filename, depth, 
                            write_intermediate_results, spiral_look,
                            output_base_name, min_distance, max_distance, use_radial_slicing,
                            red_channel_value, green_channel_value, blue_channel_value,
                            eps_factor, is_verbose);

    return 0;
}