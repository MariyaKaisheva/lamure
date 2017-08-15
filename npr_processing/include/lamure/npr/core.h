#ifndef CORE_H
#define CORE_H


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include <limits>
#include <random>

#include <scm/gl_core/math.h>

#include <lamure/bounding_box.h>
#include <lamure/types.h>

#include <lamure/ren/bvh.h>
#include <lamure/ren/lod_stream.h>
#include <lamure/ren/model_database.h>
#include <lamure/ren/dataset.h>

#include <lamure/npr/input_output.h>
#include <lamure/npr/binning.h>
#include <lamure/npr/utils.h>
#include <lamure/npr/sampling.h>
#include <lamure/npr/point.h>
#include <lamure/npr/line_gen.h>

namespace npr {
namespace core {
    void generate_line_art(scm::math::mat4f const& user_defined_rot_mat, 
                            std::string const& bvh_filename, 
                            int32_t depth, 
                            bool write_obj_file,
                            bool use_nurbs,
                            bool apply_alpha_shapes,
                            uint32_t max_number_line_loops,
                            bool is_verbose = false);
} //namespace core
} //namespace npr
#endif //CORE_H