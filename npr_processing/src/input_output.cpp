#include <lamure/npr/input_output.h>


namespace npr {
namespace io {

 bool check_user_input(char** argv, int argc) {
	std::map<std::string, std::string> valid_options; 
	prepare_options_with_descriptions(valid_options);

	bool valid_input = true;
	for (int i = 0; i < argc; ++i){
		std::string current_option = argv[i];

		if(current_option.find_first_of("-") == 0){
			auto option_it = valid_options.find(current_option);
			if(option_it == std::end(valid_options) ){
				std::cout << "INVALID OPTION: " << current_option << "\n";
				valid_input = false;
			}
		}
	}

	return valid_input; 
}

 bool cmd_option_exists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

std::string create_output_base_name(std::string const& base_name_without_path_and_extension, int32_t bvh_depth, 
                                    	   float rot_angle, scm::math::vec3d const& rot_axis,
                                    	   bool spiral_look, bool radial_slicing,
                                    	   float min_bin_distance, float max_bin_distance, 
                                    	   float DBSCAN_epsilon )  {

    auto reset_stringstream = [] (std::stringstream& strstr_to_reset) {
        strstr_to_reset.str( std::string() );
        strstr_to_reset.clear();
    };

    std::string output_base_name = base_name_without_path_and_extension;

    std::stringstream conversion_strstr;
    conversion_strstr << std::fixed << std::setprecision(2) << rot_angle;
    std::string fixed_precision_angle_string;
    conversion_strstr >> fixed_precision_angle_string;

    reset_stringstream(conversion_strstr);
    conversion_strstr << std::fixed << std::setprecision(2) << rot_axis[0];
    std::string fixed_precision_rot_axis_x_string;
    conversion_strstr >> fixed_precision_rot_axis_x_string;

    reset_stringstream(conversion_strstr);
    conversion_strstr.clear();
    conversion_strstr << std::fixed << std::setprecision(2) << rot_axis[1];
    std::string fixed_precision_rot_axis_y_string;
    conversion_strstr >> fixed_precision_rot_axis_y_string;

    reset_stringstream(conversion_strstr);
    conversion_strstr << std::fixed << std::setprecision(2) << rot_axis[2];
    std::string fixed_precision_rot_axis_z_string;
    conversion_strstr >> fixed_precision_rot_axis_z_string;

    reset_stringstream(conversion_strstr);
    conversion_strstr << std::fixed << std::setprecision(2) << DBSCAN_epsilon;
    std::string fixed_precision_epsilon_string;
    conversion_strstr >> fixed_precision_epsilon_string;


    output_base_name +=   "_d"         + std::to_string(bvh_depth) //TODO fix name for leaf level
                        + "_angle_"    + fixed_precision_angle_string
                        + "_Ax_"       + fixed_precision_rot_axis_x_string
                        + "_Ay_"       + fixed_precision_rot_axis_y_string
                        + "_Az_"       + fixed_precision_rot_axis_z_string
                        + "_eps_"      + fixed_precision_epsilon_string;


    if( !(min_bin_distance < 0) ) {
    reset_stringstream(conversion_strstr);
        conversion_strstr << std::fixed << std::setprecision(2) << min_bin_distance;
        std::string fixed_precision_min_distance_string;
        conversion_strstr >> fixed_precision_min_distance_string;
        output_base_name += "_min_"      + fixed_precision_min_distance_string;
  
    }
    
    if( !(max_bin_distance < 0) ) {
    reset_stringstream(conversion_strstr);
        conversion_strstr << std::fixed << std::setprecision(2) << max_bin_distance;
        std::string fixed_precision_max_distance_string;
        conversion_strstr >> fixed_precision_max_distance_string;
        output_base_name += "_max_"      + fixed_precision_max_distance_string;
    }

    output_base_name += "_r"  + std::to_string(radial_slicing);

    output_base_name += "_s"  + std::to_string(spiral_look);


    return output_base_name;
}


char* get_cmd_option(char** begin, char** end, const std::string & option) {
    char** it = std::find(begin, end, option);
    if (it != end && ++it != end)
        return *it;
    return 0;
}


void parse_float_parameter(int argc, char** argv, float& float_parameter, std::string const& parameter_name ) {
	if(io::cmd_option_exists(argv, argv+argc, parameter_name)){
 	float_parameter = atof(io::get_cmd_option(argv, argv+argc, parameter_name)); //user input
	}
}

 void parse_int_parameter(int argc, char** argv, int& int_parameter, std::string const& parameter_name ) {
    if(io::cmd_option_exists(argv, argv+argc, parameter_name)){
     int_parameter = atof(io::get_cmd_option(argv, argv+argc, parameter_name)); //user input
    }
}

 void prepare_options_with_descriptions(std::map<std::string, std::string>& options_with_descriptions_vec){
	//TODO ubdate the --no... flags
	options_with_descriptions_vec.emplace("-j",                                ": (optional) num threads to be used");
	options_with_descriptions_vec.emplace("-f", 				   ": (REQUIRED) specify .bvh input file");
	options_with_descriptions_vec.emplace("-t", 				   ": (optional) specify .rot input file that contains rotation for slicing plane");
	options_with_descriptions_vec.emplace("-d", 				   ": (optional) specify depth to extract; default value is the maximal depth, i.e. leaf level");
	options_with_descriptions_vec.emplace("-r", 				   ": (optional) set slicing mode to radial");
	options_with_descriptions_vec.emplace("-v", 				   ": (optional) set flag for print-outs to TRUE");
	options_with_descriptions_vec.emplace("--red", 				   ": (optional) set color value (float: 0.0 - 1.0) for red channel of lines");
	options_with_descriptions_vec.emplace("--green", 			   ": (optional) set color value (float: 0.0 - 1.0) for green channel of lines");
	options_with_descriptions_vec.emplace("--blue", 			   ": (optional) set color value (float: 0.0 - 1.0) for blue channel of lines");
	options_with_descriptions_vec.emplace("--eps",                 ": (optional) set value for DBSCAN epsylon paramerter");
	options_with_descriptions_vec.emplace("--min",                 ": (optional) set value for the minimal distance between 2 layers");
	options_with_descriptions_vec.emplace("--max",                 ": (optional) set value for the maximal distance between 2 layers");
	options_with_descriptions_vec.emplace("--verbose",             ": (optional) set flag for print-outs to TRUE");
	options_with_descriptions_vec.emplace("--write_stages",        ": (optional) set flag writing out intermediate results to TRUE");
	options_with_descriptions_vec.emplace("--generate_spirals",    ": (optional) set flag for spiral look to TRUE");

	options_with_descriptions_vec.emplace("--create_random_axes",  ": (optional) creates three perpendicular axes based on a random drawing and exits the program");
}


 void print_help_message(char** argv){
	std::map<std::string, std::string> options_with_descriptions; 
	prepare_options_with_descriptions(options_with_descriptions);
	std::cout << "Usage of " << argv[0] << std::endl <<
	"Parameters: " << std::endl; 

	for(auto const& option_description_pair : options_with_descriptions){
		std::cout << "\t" << option_description_pair.first
				  <<  option_description_pair.second
				  << std::endl;
	}

	std::cout << std::endl;	
}


scm::math::mat4f read_in_transformation_file(std::string input_filename){
	std::ifstream in_stream(input_filename);
    std::string line_buffer;
    float angle, axis_x, axis_y, axis_z; 

	while(std::getline(in_stream, line_buffer)) {
		std::stringstream my_ss(line_buffer);		
		my_ss >> angle;
		my_ss >> axis_x;
		my_ss >> axis_y;
		my_ss >> axis_z;
	}

	in_stream.close();
	
	return scm::math::make_rotation(angle, axis_x, axis_y, axis_z);
}

 void write_intermediate_result_out(std::string output_filename,
										  float avg_min_distance,
										  std::vector< std::shared_ptr<std::vector<clusters_t>> > const& all_clusters_per_bin_vector_for_all_slices,
                      scm::math::mat4f const& transformation_mat,
										  bool binning){

  uint32_t const outer_vector_size = all_clusters_per_bin_vector_for_all_slices.size();

  std::vector< std::shared_ptr<std::vector<clusters_t>> > deep_copied_all_clusters_per_bin_vector_for_all_slices(outer_vector_size);

  for( uint32_t bin_vector_id = 0; bin_vector_id < outer_vector_size; ++bin_vector_id) {

    deep_copied_all_clusters_per_bin_vector_for_all_slices[bin_vector_id] = std::make_shared< std::vector<clusters_t> >( *(all_clusters_per_bin_vector_for_all_slices[bin_vector_id])  );
    //*(deep_copied_all_clusters_per_bin_vector_for_all_slices[bin_vector_id]) = *(all_clusters_per_bin_vector_for_all_slices[bin_vector_id]);
  }

  //auto all_clusters_per_bin_vector_for_all_slices_transformed = all_clusters_per_bin_vector_for_all_slices;
  utils::transform_intermediate_points_to_original_model_orientation(deep_copied_all_clusters_per_bin_vector_for_all_slices,
                                                                    transformation_mat);


	std::ofstream output_file(output_filename);

	if(output_file.is_open()){
		//output_file <<"o testPOB" << std::endl;
		uint32_t current_cluster_id = 0;
		uint32_t current_bin_id = 0;
		int32_t current_cluster_color_id = -1;
		float thickness = avg_min_distance * 0.25f;
		for (uint32_t bin_index = 0; bin_index < deep_copied_all_clusters_per_bin_vector_for_all_slices.size(); ++bin_index){
			auto const& all_clusters_per_bin_vector = deep_copied_all_clusters_per_bin_vector_for_all_slices[bin_index];
			++current_bin_id;
			for(uint32_t cluster_index = 0; cluster_index < all_clusters_per_bin_vector->size(); ++cluster_index){
				output_file << "o line_object_" << current_cluster_id << "\n";
				if(binning){
					current_cluster_color_id = id_to_color_hash(current_bin_id % 2);
				}else{
					current_cluster_color_id = id_to_color_hash(current_cluster_id);
				}

				lamure::vec3b current_cluster_color = color_array[current_cluster_color_id];
				++current_cluster_id;

				auto const& cluster = all_clusters_per_bin_vector->at(cluster_index);

				for(uint32_t point_index = 0; point_index < cluster.size(); ++point_index){
					auto const& point = cluster[point_index];
					float x = point.pos_coordinates_[0];
					float y = point.pos_coordinates_[1];
					float z = point.pos_coordinates_[2];
					output_file <<"v " << x << " " << y << " " << z << "\n";
					output_file <<"c " << current_cluster_color[0] / 255.0 << " " << current_cluster_color[1] / 255.0 << " " << current_cluster_color[2] / 255.0 << "\n";
					output_file <<"t " << thickness << " \n"; 
				}

			}
          
		}
	}
	output_file.close();
}


 void write_output_obj(std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh){

      std::ofstream output_file(output_filename);
      unsigned long vert_counter = 1;
      //consider hidden translation
      const scm::math::vec3f& translation = bvh->get_translation();
      if (output_file.is_open()){
          for (uint i = 0; i < line_data.size(); ++i){

           output_file << "v " << std::setprecision(DEFAULT_PRECISION) << translation.x + line_data.at(i).start.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << translation.y + line_data.at(i).start.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION)<< translation.z + line_data.at(i).start.pos_coordinates_[2] << "\n";
           output_file << "v " << std::setprecision(DEFAULT_PRECISION) << translation.x +line_data.at(i).end.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << translation.y + line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << translation.z + line_data.at(i).end.pos_coordinates_[2] << "\n";
           //vertex duplication to emulate triangle
           auto x_offset =  (translation.x + line_data.at(i).end.pos_coordinates_[0]) / 1000000.0;
           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  translation.x + line_data.at(i).end.pos_coordinates_[0] + x_offset << " " << std::setprecision(DEFAULT_PRECISION) << translation.y + line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << translation.z + line_data.at(i).end.pos_coordinates_[2] << "\n";
           output_file << "f " << vert_counter << " " << (vert_counter + 1) << " " << (vert_counter + 2) << "\n";
           vert_counter += 3;
          }

          output_file.close();
      }
      else{
        std::cout << "<LAMURE_NPR_PROCESSING>: Cannot open output file to write to! \n";
      }

    std::cout << "OUTPUT FILE: " << output_filename << std::endl;
}


 void write_output_lob(std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh, float avg_min_distance, float r, float g, float b){

      std::ofstream output_file(output_filename);

      float thickness = avg_min_distance ;
      //consider hidden translation
	std:: cout << "THICKNESS: " << thickness << "\n";
      float float_max{std::numeric_limits<float>::max()};
      lamure::vec3f previous_end_vertex_position{float_max, float_max, float_max};
      const scm::math::vec3f& translation = bvh->get_translation();
      std::cout << "HIDDEN TRANSLATION " << translation << "\n"; 

      uint32_t total_line_object_counter = 0;
      if (output_file.is_open()){



          for (uint line_idx = 0; line_idx < line_data.size(); ++line_idx){
      	    auto const& current_start_vertex_position = line_data.at(line_idx).end.pos_coordinates_;
      	    auto const& current_end_vertex_position   = line_data.at(line_idx).start.pos_coordinates_;

            if( (previous_end_vertex_position[0] != current_start_vertex_position[0] ) ||
            	(previous_end_vertex_position[1] != current_start_vertex_position[1] ) ||
            	(previous_end_vertex_position[2] != current_start_vertex_position[2] )
              ) {

              if( 0 != line_idx) {
		        output_file << "v " << std::setprecision(DEFAULT_PRECISION) << translation.x + previous_end_vertex_position[0] << " " << std::setprecision(DEFAULT_PRECISION) << translation.y + previous_end_vertex_position[1] << " " << std::setprecision(DEFAULT_PRECISION)<< translation.z + previous_end_vertex_position[2] << "\n";
		        output_file << "c " << r << " " << g << " " << b << "\n";
		        output_file << "t " << thickness << "\n";	              	
              }

              output_file << "o line_object_" << total_line_object_counter++ << "\n";
            }

           output_file << "v " << std::setprecision(DEFAULT_PRECISION) << translation.x + current_start_vertex_position[0] << " " << std::setprecision(DEFAULT_PRECISION) << translation.y + current_start_vertex_position[1] << " " << std::setprecision(DEFAULT_PRECISION)<< translation.z + current_start_vertex_position[2] << "\n";
		   output_file << "c " << r << " " << g << " " << b << "\n";
           output_file << "t " << thickness << "\n";


           if( (line_data.size() - 1) == line_idx) {
		     output_file << "v " << std::setprecision(DEFAULT_PRECISION) << translation.x + current_end_vertex_position[0] << " " << std::setprecision(DEFAULT_PRECISION) << translation.y + current_end_vertex_position[1] << " " << std::setprecision(DEFAULT_PRECISION)<< translation.z + current_end_vertex_position[2] << "\n";
		     output_file << "c " << r << " " << g << " " << b << "\n";
		     output_file << "t " << thickness << "\n";	
           }

           //set current end to previous end vertex position before starting next loop iteration
           for(int dim_idx = 0; dim_idx < 3; ++dim_idx) {
           	previous_end_vertex_position[dim_idx] = current_end_vertex_position[dim_idx];
           }
          }

          output_file.close();
      }
      else{
        std::cout << "<LAMURE_NPR_PROCESSING>: Cannot open output file to write to! \n";
      }

    std::cout << "OUTPUT FILE: " << output_filename << std::endl;
  }

  std::vector<xyzall_surfel_t> extract_data_from_fixed_lod(std::string bvh_filename,
                                                           int32_t depth){
    lamure::ren::bvh* bvh = new lamure::ren::bvh(bvh_filename);
    std::string lod_filename = bvh_filename.substr(0, bvh_filename.size()-3) + "lod";
    lamure::ren::lod_stream* in_access = new lamure::ren::lod_stream();
    in_access->open(lod_filename);

    size_t size_of_node = (uint64_t)bvh->get_primitives_per_node() * sizeof(lamure::ren::dataset::serialized_surfel);
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
    delete in_access;
    delete bvh; 
    return surfels_vector;
  }

}
}
