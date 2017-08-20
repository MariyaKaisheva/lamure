#ifndef LINE_GEN_H
#define LINE_GEN_H


#include "nurbscurve.h"
#include "clustering.h"

#include <stack>


namespace npr {
namespace line_gen {
    struct line_approximation_job{
        float start_t_;
        float end_t_;
        gpucast::math::point3d approximated_point_;

        line_approximation_job (float start_t, float end_t, gpucast::math::point3d middle_point) : start_t_(start_t),
                                                                                                 end_t_(end_t), 
                                                                                                 approximated_point_(middle_point) {}
    };
    //sort in descending order based on y coordinate value 
    bool comparator (const xyzall_surfel_t& A, const xyzall_surfel_t& B);

    gpucast::math::point3d get_midpoint(gpucast::math::point3d const& start_point, 
                                        gpucast::math::point3d const& end_point);

    void interpolate_cluster_input(std::vector<point> & combined_cluster_points, 
                                   std::vector<uint32_t> const& original_cluster_sizes);

    std::vector<line> generate_lines_from_curve (std::vector<point> const& ordered_points, uint8_t degree = 3, bool is_verbose = false);

    std::vector<line> generate_lines(std::vector<xyzall_surfel_t>& input_data, 
                                     uint32_t& max_num_line_loops, 
                                     bool use_nurbs, 
                                     bool apply_alpha_shapes,
                                     bool is_verbose = false);
} //namespace line_gen
} //namespace npr

#endif //LINE_GEN_H