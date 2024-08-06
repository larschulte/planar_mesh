// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Alpha_shape_2.h>
// #include <CGAL/Alpha_shape_face_base_2.h>
// #include <CGAL/Alpha_shape_vertex_base_2.h>
// #include <CGAL/Regular_triangulation_2.h>
// #include <fstream>
// #include <iostream>
// #include <list>
// #include <vector>
// #include <CGAL/draw_polygon_2.h>
// #include <CGAL/Delaunay_triangulation_2.h>

// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef K::FT FT;
// typedef K::Weighted_point_2 Weighted_point;
// typedef K::Point_2 Point;
// typedef K::Segment_2 Segment;
// typedef CGAL::Regular_triangulation_vertex_base_2<K> Rvb;
// typedef CGAL::Alpha_shape_vertex_base_2<K, Rvb> Vb;
// typedef CGAL::Regular_triangulation_face_base_2<K> Rf;
// typedef CGAL::Alpha_shape_face_base_2<K, Rf> Fb;
// typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
// typedef CGAL::Regular_triangulation_2<K, Tds> Triangulation_2;
// typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
// typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/algorithm.h>

#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene_options.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

typedef K::FT                                                FT;
typedef K::Point_2                                           Point;
typedef K::Segment_2                                         Segment;

typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>          Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;

typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;

int main()
{
	// data
	// std::list<Weighted_point> wpoints = 
	// {
	// 	Weighted_point(K::Point_2(0, 0), 1), 
	// 	Weighted_point(K::Point_2(1, 0), 1), 
	// 	Weighted_point(K::Point_2(3, 0), 1), 
	// 	Weighted_point(K::Point_2(6, 0), 1),
	// 	Weighted_point(K::Point_2(10, 0), 1),
	// 	Weighted_point(K::Point_2(0, 1), 1), 
	// 	Weighted_point(K::Point_2(1, 1), 1), 
	// 	Weighted_point(K::Point_2(3, 1), 1), 
	// 	Weighted_point(K::Point_2(6, 1), 1),
	// 	Weighted_point(K::Point_2(10, 1), 1),
	// };

	std::list<Point> points = 
	{
		Point(0, 0),
		Point(1, 0),
		Point(3, 0),
		Point(6, 0),
		Point(10, 0),
		Point(0, 1),
		Point(1, 1),
		Point(3, 1),
		Point(6, 1),
		Point(10, 1),
	};
	
	Alpha_shape_2 A(points.begin(), points.end(), 2, Alpha_shape_2::GENERAL);
	std::cout << "Alpha value: " << A.get_alpha() << std::endl;

	// // compute alpha shape
	// Alpha_shape_2 A(wpoints.begin(), wpoints.end(), 2, Alpha_shape_2::GENERAL);
	// std::cout << "Alpha value: " << A.get_alpha() << std::endl;
	// std::cout << "Optimal alpha: " << *A.find_optimal_alpha(1) << std::endl;

	// extract alpha shape edges
	std::vector<Segment> segments;
	for (Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(); it != A.alpha_shape_edges_end(); ++it)
	{
		segments.push_back(A.segment(*it));
	}
	std::cout << segments.size() << " alpha shape edges" << std::endl;

	// extract alpha shape vertices
	std::vector<Point> vertices;
	for (Alpha_shape_2::Alpha_shape_vertices_iterator it = A.alpha_shape_vertices_begin(); it != A.alpha_shape_vertices_end(); ++it)
	{
		vertices.push_back(A.point(*it));
	}
	std::cout << vertices.size() << " alpha shape vertices" << std::endl;

	// draw the segments
	CGAL::Graphics_scene scene;
	for (Segment s : segments)
	{
		// random color
		CGAL::IO::Color random_color = CGAL::IO::Color(rand() % 256, rand() % 256, rand() % 256);
		scene.add_segment(s.source(), s.target(), random_color);
	}
	for (Point p : vertices)
	{
		scene.add_point(p, CGAL::IO::Color(255, 0, 0));
	}
    CGAL::draw_graphics_scene(scene);


	// return 
	return 0;
}