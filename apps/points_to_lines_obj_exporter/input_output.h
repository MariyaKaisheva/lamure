#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <lamure/ren/bvh.h>

#include "utils.h"

#define DEFAULT_PRECISION 15

namespace io{

	inline char* get_cmd_option(char** begin, char** end, const std::string & option) {
	    char** it = std::find(begin, end, option);
	    if (it != end && ++it != end)
	        return *it;
	    return 0;
	}

	inline bool cmd_option_exists(char** begin, char** end, const std::string& option) {
	    return std::find(begin, end, option) != end;
	}

	inline void write_output(bool write_obj_file, std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh){
		if(write_obj_file) {
	      std::ofstream output_file(output_filename);
	      unsigned long vert_counter = 1;

	      if (output_file.is_open()){
	          for (uint i = 0; i < line_data.size(); ++i){

	           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).start.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).start.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION)<< line_data.at(i).start.pos_coordinates_[2] << "\n";
	           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.pos_coordinates_[0] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[2] << "\n";
	           //vertex duplication to emulate triangle
	           auto x_offset =  line_data.at(i).end.pos_coordinates_[0] / 1000000.0;
	           output_file << "v " << std::setprecision(DEFAULT_PRECISION) <<  line_data.at(i).end.pos_coordinates_[0] + x_offset<< " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[1] << " " << std::setprecision(DEFAULT_PRECISION) << line_data.at(i).end.pos_coordinates_[2] << "\n";
	           output_file << "f " << vert_counter << " " << (vert_counter + 1) << " " << (vert_counter + 2) << "\n";
	           vert_counter += 3;
	          }

	          output_file.close();
	      }
	      else{
	        std::cout << "Cannot open output file to write to! \n";
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
	        std::cout << "Cannot open output file to write to! \n";
	      }
	    }
	}
}



#endif //INPUT_OUTPUT_H