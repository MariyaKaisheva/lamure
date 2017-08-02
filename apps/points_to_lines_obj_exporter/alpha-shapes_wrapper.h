#ifndef ALPHA_SHAPES_H_
#define ALPHA_SHAPES_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/algorithm.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

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
typedef CGAL::Quotient<CGAL::MP_Float> Number_type;
typedef CGAL::Cartesian<Number_type> Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel> Traits_2;
typedef Traits_2::Point_2 Monotone_Curve_Vertex;
typedef Traits_2::X_monotone_curve_2  Monotone_Curve_Segment;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;



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
	Alpha_shape_2 A(points.begin(), points.end(), FT(0.5), Alpha_shape_2::GENERAL);
	Alpha_shape_2 A_second_run(points.begin(), points.end(), FT(*A.find_optimal_alpha(1)), Alpha_shape_2::GENERAL);
 	std::vector<Segment> segments;
  	//alpha_edges( A_second_run, std::back_inserter(segments));
  	alpha_edges( A_second_run, std::back_inserter(segments));
 	std::vector<Monotone_Curve_Segment> segments_as_monotone_curves;

 	for( auto const& segment_to_convert : segments) {

 		Monotone_Curve_Vertex start_segment_point(segment_to_convert.point(0).hx(), segment_to_convert.point(0).hy());
  		Monotone_Curve_Vertex end_segment_point  (segment_to_convert.point(1).hx(), segment_to_convert.point(1).hy());
 		segments_as_monotone_curves.emplace_back(start_segment_point, end_segment_point);
 	}

 	segments.clear();

  	Arrangement_2 arrangement;
  	CGAL::insert( arrangement, &segments_as_monotone_curves[0], &segments_as_monotone_curves[segments_as_monotone_curves.size()]);
  	Arrangement_2::Face_const_handle unbounded_face = arrangement.unbounded_face();
  	  // Print the boundary of each of the holes.
	 Arrangement_2::Hole_const_iterator hole_iterator;
	  int index = 1;

	  for (hole_iterator = unbounded_face->holes_begin(); hole_iterator != unbounded_face->holes_end(); ++hole_iterator, ++index) {
	    std::cout << " Hole #" << index << ": \n";

		Arrangement_2::Ccb_halfedge_const_circulator start_of_hole_circulator = *hole_iterator;
	  	Arrangement_2::Ccb_halfedge_const_circulator halfedge_circulator = start_of_hole_circulator;

	  	uint32_t num_edges_in_hole = 0;
	  	do{
	  		++num_edges_in_hole;
	  		auto point_0 = halfedge_circulator->source()->point();
	  		auto point_1 = halfedge_circulator->target()->point();

	  		cgal_point_2 converted_point_0(CGAL::to_double(point_0.hx()), CGAL::to_double(point_0.hy()) );
	  		cgal_point_2 converted_point_1(CGAL::to_double(point_1.hx()), CGAL::to_double(point_1.hy()) );

	  		segments.emplace_back(converted_point_0, converted_point_1);
	  	} while(++halfedge_circulator != start_of_hole_circulator);

	  	std::cout <<"num edges in hole: " << num_edges_in_hole << "\n";
	 }



  	return segments;
}

std::vector<line> do_output_conversion_to_line_segments(std::vector<Segment> const& segments, float prjected_coord){
	std::vector<line> line_data;
	for(uint i = 0; i < segments.size(); ++i){
		auto const& current_segment = segments[i];
		auto const& current_cgal_start_point = current_segment.point(0);
		auto const& current_cgal_end_point = current_segment.point(1);

		float start_point_coordinates[3] = {float(current_cgal_start_point.hx()), prjected_coord, float(current_cgal_start_point.hy())};
		float end_point_coordinates[3] = {float(current_cgal_end_point.hx()), prjected_coord, float(current_cgal_end_point.hy())};

		point start_point(start_point_coordinates);
		point end_point(end_point_coordinates);
		auto start_coord = lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]);
		auto end_coord = lamure::vec3f(end_point.pos_coordinates_[0], end_point.pos_coordinates_[1], end_point.pos_coordinates_[2]);
		auto length = utils::compute_distance(start_coord, end_coord);

		line current_line(start_point, end_point, length);
		line_data.push_back(current_line);

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