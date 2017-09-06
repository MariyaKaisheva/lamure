#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <map>
#include <memory>

#include <lamure/ren/bvh.h>
#include <lamure/npr/binning.h>
#include <lamure/npr/clustering.h>
#include <lamure/npr/color_hash_map.h>

#include <scm/gl_core/math/mat4_gl.h>

#include "utils.h"

#define DEFAULT_PRECISION 15

namespace npr {
namespace io {

	struct stage_content_storage{

		std::vector<binning::bin> binns_;
		std::vector< std::shared_ptr<std::vector<clusters_t>> > clusters_;
		std::vector<std::vector< std::shared_ptr<std::vector<point> > > > alpha_shapes_;
	};

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

		options_with_descriptions_vec.emplace("-f", 				   ": (REQUIRED) specify .bvh input file");
		options_with_descriptions_vec.emplace("-t", 				   ": (optional) specify .rot input file that contains rotation for slicing plane");
		options_with_descriptions_vec.emplace("-d", 				   ": (optional) specify depth to extract; default value is the maximal depth, i.e. leaf level");
		options_with_descriptions_vec.emplace("-s", 				   ": (optional) specify output stage");
		options_with_descriptions_vec.emplace("-v", 				   ": (optional) set flag for print-outs to TRUE");
		options_with_descriptions_vec.emplace("--min",                 ": (optional) set value for the minimal distance between 2 layers");
		options_with_descriptions_vec.emplace("--max",                 ": (optional) set value for the maximal distance between 2 layers");
		options_with_descriptions_vec.emplace("--no_reduction",        ": (optional) set flag reduce num slicing layers proporional to selected LoD to FALSE");
		options_with_descriptions_vec.emplace("--apply_nurbs_fitting", ": (optional) set flag for curve-fitting to TRUE");
		options_with_descriptions_vec.emplace("--verbose",             ": (optional) set flag for print-outs to TRUE");
		options_with_descriptions_vec.emplace("--apply_alpha_shapes",  ": (optional) set flag for alpha-shaped to TRUE");
		options_with_descriptions_vec.emplace("--write_stages",        ": (optional) set flag writing out intermediate results to TRUE");
		options_with_descriptions_vec.emplace("--generate_spirals",    ": (optional) set flag for spiral look to TRUE");
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

	inline void write_intermediate_result_out(std::string output_filename,
											  float avg_min_distance,
											  std::vector< std::shared_ptr<std::vector<clusters_t>> > const& all_clusters_per_bin_vector_for_all_slices,
											  bool binning = false){

		std::ofstream output_file(output_filename);

		if(output_file.is_open()){
			output_file <<"o testPOB" << std::endl;
			uint32_t current_cluster_id = 0;
			uint32_t current_bin_id = 0;
			int32_t current_cluster_color_id = -1;
			float thickness = avg_min_distance * 0.25f;
			for (uint32_t bin_index = 0; bin_index < all_clusters_per_bin_vector_for_all_slices.size(); ++bin_index){
				auto const& all_clusters_per_bin_vector = all_clusters_per_bin_vector_for_all_slices[bin_index];
				++current_bin_id;
				for(uint32_t cluster_index = 0; cluster_index < all_clusters_per_bin_vector->size(); ++cluster_index){
					
					if(binning){
						current_cluster_color_id = id_to_color_hash(current_bin_id);
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


	inline void write_output(std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh){

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

} //namespace io
} //namespace npr


#endif //INPUT_OUTPUT_H