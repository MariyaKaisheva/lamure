#ifndef ALPHA_SHAPES_H_
#define ALPHA_SHAPES_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Cartesian.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <list>

namespace alpha{


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_2  cgal_point_2;
typedef K::Segment_2  Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K>  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

template <class OutputIterator>
void
alpha_edges( const Alpha_shape_2&  A,
         OutputIterator out)
{
  for(Alpha_shape_edges_iterator it =  A.alpha_shape_edges_begin();
      it != A.alpha_shape_edges_end();
      ++it){
    *out++ = A.segment(*it);
  }
}

void do_input_conversion(std::vector<cgal_point_2> & cgal_points_vec, std::vector<point> const& points_vec){
	
	for(uint point_itr = 0; point_itr < points_vec.size(); ++point_itr){
		auto const& current_point = points_vec[point_itr];
		auto hx = current_point.pos_coordinates_[0];
		auto hz = current_point.pos_coordinates_[2];
		cgal_point_2 current_cgal_point(hx, hz);
		cgal_points_vec[point_itr] = current_cgal_point;
	}
}

std::vector<Segment> generate_alpha_shape(std::vector<cgal_point_2> const& points){
	Alpha_shape_2 A(points.begin(), points.end(), FT(10), Alpha_shape_2::REGULARIZED);
	Alpha_shape_2 A_second_run(points.begin(), points.end(), FT(*A.find_optimal_alpha(1)), Alpha_shape_2::REGULARIZED);
 	std::vector<Segment> segments;
  	alpha_edges( A_second_run, std::back_inserter(segments));
  	//std::cout << segments.size() << " alpha shape edges" << std::endl;
  	std::cout << "Optimal alpha 1st run: " << *A.find_optimal_alpha(1)<<std::endl;
  	//std::cout << "Optimal alpha 2nd run: " << *A_second_run.find_optimal_alpha(1)<<std::endl;

  	return segments;
}

std::vector<line> do_output_conversion_to_line_segments(std::vector<Segment> const& segments, float prjected_coord){
	std::vector<line> line_data;
	for(uint i = 0; i < segments.size() - 1; ++i){
		auto const& current_segment = segments[i];
		auto const& current_cgal_start_point = current_segment.point(0);
		auto const& current_cgal_end_point = current_segment.point(1);

		auto const& next_cgal_segment = segments[i + 1];
		auto const& next_cgal_start_point = next_cgal_segment.point(0);
		auto const& next_cgal_end_point = next_cgal_segment.point(1);

		/* if( ( (current_cgal_end_point.hx() == next_cgal_start_point.hx()) && (current_cgal_end_point.hy() == next_cgal_start_point.hy()) ) ||
			( (current_cgal_start_point.hx() == next_cgal_end_point.hx()) && (current_cgal_start_point.hy() == next_cgal_end_point.hy()) )

		  ){*/
				float start_point_coordinates[3] = {float(current_cgal_start_point.hx()), prjected_coord, float(current_cgal_start_point.hy())};
				float end_point_coordinates[3] = {float(current_cgal_end_point.hx()), prjected_coord, float(current_cgal_end_point.hy())};

				point start_point(start_point_coordinates);
				point end_point(end_point_coordinates);
				auto start_coord = lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]);
				auto end_coord = lamure::vec3f(end_point.pos_coordinates_[0], end_point.pos_coordinates_[1], end_point.pos_coordinates_[2]);
				auto length = utils::compute_distance(start_coord, end_coord);

				line current_line(start_point, end_point, length);
				line_data.push_back(current_line);
		/*}else{
			break;
		}*/

	}
	return line_data;

}

std::vector<point> do_output_conversion(std::vector<Segment> const& segments, float prjected_coord){
	std::vector<point> clened_sorted_points_vec;
	for(uint i = 0; i < segments.size(); ++i){
		auto const& current_segment = segments[i];
		auto cgal_point = current_segment.point(0);
		float point_coord[3] = {float(cgal_point.hx()), prjected_coord, float(cgal_point.hy())}; 
		point current_point(point_coord);
		clened_sorted_points_vec.push_back(current_point);
	}
	return clened_sorted_points_vec;
}

}


#endif //ALPHA_SHAPES_H_