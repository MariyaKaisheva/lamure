#ifndef LINE_GEN_H
#define LINE_GEN_H


#include "nurbscurve.h"
#include <lamure/npr/clustering.h>
#include <lamure/npr/binning.h>
#include <lamure/npr/input_output.h>

#include <stack>
#include <memory>


namespace npr {
namespace line_gen {

using nurbs_vec_t = std::vector<gpucast::math::nurbscurve3d>;

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

    void rotate(clusters_t & point_cluster);

    gpucast::math::nurbscurve3d fit_curve(std::vector<point> const& ordered_points, uint8_t degree,  bool spiral_look, bool is_verbose);

    //std::vector<line> generate_lines_from_curve (std::vector<point> const& ordered_points, uint8_t degree = 3, bool is_verbose = false);

    std::vector<line> evaluate_curve(gpucast::math::nurbscurve3d& nurbs_curve, bool dynamic_sampling_step);

    std::vector<point> blend_between_curves(gpucast::math::nurbscurve3d& top_curve, 
                                            gpucast::math::nurbscurve3d& bottom_curve,
                                            float max_distance);

    nurbs_vec_t generate_spirals(std::vector<nurbs_vec_t> const& guiding_nurbs_vec, float max_distance);

    void prepare_clusters (std::vector<binning::bin> & bin_vec, 
                           std::vector< std::shared_ptr<std::vector<clusters_t>> > & all_clusters_per_bin_vector_for_all_slices,
                           uint32_t num_cells_pro_dim,
                           bool is_verbose);

    void clean_clusters_via_alpha_shape_detection(std::vector< std::shared_ptr<std::vector<clusters_t>> > & all_clusters_per_bin_vector_for_all_slices,
                                                  //std::vector<std::vector< std::shared_ptr<std::vector<point> > > > & all_alpha_shapes_for_all_bins,
                                                  std::vector< std::shared_ptr<std::vector<clusters_t>> >  & all_alpha_shapes_for_all_bins,
                                                  uint32_t degree,
                                                  bool color,
                                                  bool is_verbose);

    std::vector<line> generate_lines(std::vector<xyzall_surfel_t>& input_data, 
                                     float min_distance,
                                     float max_distance,
                                     float& out_avg_min_distance,
                                     std::string output_base_name,
                                     bool write_intermediate_results,
                                     float eps_factor = 10.0,
                                     bool use_nurbs = true, 
                                     bool apply_alpha_shapes = true,
                                     bool spiral_look = false,
                                     bool is_verbose = false);



} //namespace line_gen
} //namespace npr

#endif //LINE_GEN_H