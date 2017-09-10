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

	//intermediate stage output types
	struct stage_content_storage{

		std::vector<binning::bin> binns_;
		std::vector< std::shared_ptr<std::vector<clusters_t>> > clusters_;
		std::vector<std::vector< std::shared_ptr<std::vector<point> > > > alpha_shapes_;
	};

	//verifies that no options were set that can not be used in the npr creation
	bool check_user_input(char** argv, int argc);

	//creates formatted output name based on the input parameters
	std::string create_output_base_name(std::string const& base_name_without_path_and_extension, int32_t bvh_depth, 
	                                    float rot_angle, scm::math::vec3d const& rot_axis,
	                                    bool spiral_look, 
	                                    float min_bin_distance, float max_bin_distance, 
	                                    float DBSCAN_epsilon);

	//checks if an argument is set (for bool parameters)
	bool cmd_option_exists(char** begin, char** end, std::string const& option);

	//returns the parsed options
	char* get_cmd_option(char** begin, char** end, std::string const& option);

	//helper function to facilitate float parameter parsing
	void parse_float_parameter(int argc, char** argv, float& float_parameter, std::string const& parameter_name );

	//helper function to facilitate int parameter parsing
	void parse_int_parameter(int argc, char** argv, int& int_parameter, std::string const& parameter_name );

	void prepare_options_with_descriptions(std::map<std::string, std::string>& options_with_descriptions_vec);

	void print_help_message(char** argv);


	scm::math::mat4f read_in_transformation_file(std::string input_filename);

	//writes out the intermediate stages as colored *.pob objects to render points in guacamole
	void write_intermediate_result_out(std::string output_filename,
										float avg_min_distance,
										std::vector< std::shared_ptr<std::vector<clusters_t>> > const& all_clusters_per_bin_vector_for_all_slices,
										bool binning = false);

	//writes an *.obj file with the final result that contains degenerated triangles to be seen as lines
	void write_output_obj(std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh);

	//writes a *.lob file with the final result that adheres to the *.lob format to render line objects in guacamole
	void write_output_lob(std::string output_filename, std::vector<line> const& line_data, lamure::ren::bvh* bvh, float avg_min_distance, float r = 1.0, float g = 1.0, float b = 1.0);

} //namespace io
} //namespace npr


#endif //INPUT_OUTPUT_H