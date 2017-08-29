#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <map>

#include <lamure/ren/bvh.h>

#include "utils.h"

#define DEFAULT_PRECISION 15

namespace npr {
namespace io {

	inline char* get_cmd_option(char** begin, char** end, const std::string & option) {
	    char** it = std::find(begin, end, option);
	    if (it != end && ++it != end)
	        return *it;
	    return 0;
	}

	inline bool cmd_option_exists(char** begin, char** end, const std::string& option) {
	    return std::find(begin, end, option) != end;
	}

	inline void prepare_options_with_descriptions(std::map<std::string, std::string>& options_with_descriptions_vec){

		options_with_descriptions_vec.emplace("-f", ": (REQUIRED) specify .bvh input file");
		options_with_descriptions_vec.emplace("-t", ": (optional) specify .rot input file that contains rotation for slicing plane");
		options_with_descriptions_vec.emplace("-d", ": (optional) specify depth to extract; default value is the maximal depth, i.e. leaf level");
		options_with_descriptions_vec.emplace("-l", ": (optional) specify max number of slicing layers; vaule should be more than 5");
		options_with_descriptions_vec.emplace("-v", ": (optional); set flag for print-outs to TRUE");
		options_with_descriptions_vec.emplace("--no_reduction", ": (optional) set flag reduce num slicing layers proporional to selected LoD to FALSE");
		options_with_descriptions_vec.emplace("--apply_nurbs_fitting", ": (optional) set flag for curve-fitting to TRUE");
		options_with_descriptions_vec.emplace("--verbose", ": (optional); set flag for print-outs to TRUE");
		options_with_descriptions_vec.emplace("--apply_alpha_shapes", ": (optional) set flag for alpha-shaped to TRUE");
		options_with_descriptions_vec.emplace("--write_xyz_points", ": (optional) writes an xyz_point_cloud instead of a *.obj containing line data");
		options_with_descriptions_vec.emplace("--generate_spirals", ": (optional) set flag for spiral look to TRUE");
	}


	inline void print_help_message(char** argv){
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


	inline bool check_user_input(char** argv, int argc) {
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

	inline scm::math::mat4f read_in_transformation_file(std::string input_filename){
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

	inline void write_output(bool write_obj_file, std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh){
		if(write_obj_file) {

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

	    }
	    else {
	      //consider hidden translation
	      const scm::math::vec3f& translation = bvh->get_translation();
	      std::ofstream output_file(output_filename);
	      lamure::vec3f const fixed_upward_normal(0.0, 1.0, 0.0);
	      lamure::vec3f const fixed_forward_normal(0.0, 0.0, 1.0);
	      float const fixed_radius(0.02);
	      if (output_file.is_open()){
	          for (uint i = 0; i < line_data.size(); ++i){
	           output_file << std::setprecision(DEFAULT_PRECISION) << translation.x + line_data.at(i).start.pos_coordinates_[0] << " ";
	           output_file << std::setprecision(DEFAULT_PRECISION) << translation.y + line_data.at(i).start.pos_coordinates_[1] << " ";
	           output_file << std::setprecision(DEFAULT_PRECISION) << translation.z + line_data.at(i).start.pos_coordinates_[2] << " ";
	           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_upward_normal.x << " ";
	           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_upward_normal.y << " ";
	           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_upward_normal.z << " ";
	           output_file << (int) line_data.at(i).start.r_  << " ";
	           output_file << (int) line_data.at(i).start.g_ << " ";
	           output_file << (int) line_data.at(i).start.b_ << " ";  
	           output_file << std::setprecision(DEFAULT_PRECISION) << fixed_radius << std::endl;
	          }
	        
	          output_file.close();
	      }
	      else{
	        std::cout << "<LAMURE_NPR_PROCESSING>: Cannot open output file to write to! \n";
	      }
	    }

	    std::cout << "<LAMURE_NPR_PROCESSING>: Output: " << output_filename << std::endl;
	}

} //namespace io
} //namespace npr


#endif //INPUT_OUTPUT_H