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

#include <lamure/ren/model_database.h>
#include <lamure/bounding_box.h>
#include <lamure/types.h>
#include <lamure/ren/dataset.h>
#include <lamure/ren/bvh.h>
#include <lamure/ren/lod_stream.h>

#include "input_output.h"
#include "binning.h"
#include "utils.h"
#include "sampling.h"
#include "math_wrapper.h"
#include "point.hpp"
#include "line_gen.h"

namespace core{
    void generate_line_art(scm::math::mat4f user_defined_rot_mat, 
                            std::string bvh_filename, 
                            int32_t depth, 
                            bool write_obj_file,
                            bool use_nurbs,
                            bool apply_alpha_shapes,
                            uint32_t max_number_line_loops);
}

#endif //CORE_H