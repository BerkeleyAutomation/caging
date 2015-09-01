// Copyright (c) 2000, 2001, 2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.5-branch/Geomview/demo/Geomview/gv_terrain.cpp $
// $Id: gv_terrain.cpp 47395 2008-12-12 08:37:38Z afabri $
//
//
// Author(s)     : Sylvain Pion

#include <CGAL/Cartesian.h>
#include <iostream>
#include <Alpha_shapes_disconnection.h>    
#include <CGAL/Timer.h>

#ifndef CGAL_USE_GEOMVIEW
int main()
{
  std::cout << "Geomview doesn't work on Windows, so..." << std::endl;
  return 0;
}
#else

#include <fstream>
#include <unistd.h> // for sleep()

#include <CGAL/Triangulation_euclidean_traits_xy_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_2.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <CGAL/intersections.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>
#include <list>
#include <math.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Gt;

typedef CGAL::Timer                         Timer;

typedef CGAL::Alpha_shape_vertex_base_3<Gt>         Vb;
typedef CGAL::Alpha_shape_cell_base_3<Gt>           Fb;
typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_3<Gt,Tds>       Triangulation_3;
typedef CGAL::Alpha_shape_3<Triangulation_3>        Alpha_shape_3;

typedef Alpha_shape_3::Cell_handle          Cell_handle;
typedef Alpha_shape_3::Vertex_handle        Vertex_handle;
typedef Alpha_shape_3::Facet                Facet;
typedef Alpha_shape_3::Edge                 Edge;
typedef Gt::Weighted_point                  Weighted_point;
typedef Gt::Bare_point                      Bare_point;



int main()
{


  Bare_point qs(-.01,.01,.01);
  Bare_point qg(.01,.01,-.01);

  std::list<Weighted_point> lwp, lwpint;


  Timer t1;
  Timer t2;

  //input : a small molecule

  /*lwp.push_back(Weighted_point(Bare_point( 1, -1, -1), 4));
  lwp.push_back(Weighted_point(Bare_point(-1,  1, -1), 4));
  lwp.push_back(Weighted_point(Bare_point(-1, -1,  1), 4));
  lwp.push_back(Weighted_point(Bare_point( 1,  1,  1), 4));
  lwp.push_back(Weighted_point(Bare_point( 2,  2,  2), 5));
  lwp.push_back(Weighted_point(Bare_point( 2,  2,  3), 5));
  lwp.push_back(Weighted_point(Bare_point( 3,  2,  1), 5));
  lwp.push_back(Weighted_point(Bare_point( 10,  10,  10), 1));*/
  
  //input : a sphere of spheres
  t1.start();

  t2.start();
  lwp = generate_sphere_of_spheres(-3.0, -3.0, -3.0, 1.5 , 0.1 , 30, 40);
  lwpint = generate_sphere_of_spheres(3.0, 3.0, 3.0, 1.5 , 0.1 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(-3.0, -3.0, 3.0, 2 , 0.1 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(-3.0, 3.0, -3.0, 1.5 , 0.1 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(3.0, -3.0, -3.0, 1.5 , 0.075 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(3.0, 3.0, -3.0, 1.5 , 0.075 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(3.0, -3.0, 3.0, 1.5 , 0.1 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(-3.0, 3.0, 3.0, 1.5 , 0.1 , 30, 40);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  lwpint = generate_sphere_of_spheres(0.0, 0.0, 0.0, 2.5 , 0.15 , 35, 45);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());
  t2.stop();


  lwpint = generate_random_spheres(-5,-5,-5,5,5,5, .2, 1.2, 10000);
  lwp.insert(lwp.end(),lwpint.begin(),lwpint.end());

  std::cout << "Sphere gen time is: " << t2.time() << std::endl;
  std::cout << "Num_spheres is: " << lwp.size() << std::endl;  
  t2.reset();
  

  //build alpha_shape  in GENERAL mode and set alpha=0

  t2.start();
  Triangulation_3 triangulation(lwp.begin(), lwp.end());
  //Alpha_shapes_disconnection  asd(lwp.begin(), lwp.end());
  t2.stop();
  
  std::cout << "Regular Triangulation gen time is: " << t2.time() << std::endl;
  std::cout << "Number of cells in Triangulation: " << triangulation.number_of_finite_cells() << std::endl;
  t2.reset();

  t2.start();
  Alpha_shapes_disconnection  asd(triangulation);
  t2.stop();
  
  std::cout << "Alpha Shape and struct gen time is: " << t2.time() << std::endl;
  t2.reset();

  std::cout << "Number of connected components: " << asd.Num_sets() << std::endl;

  

/*  std::cout << "Truth or dare: ";
  t2.start();
  if(asd.Query_disconnection(qs,qg)){  //test disconnection
    t2.stop();
    std::cout << "Connected";}
  else{
    t2.stop();
    std::cout << "Disconnected";}
  std::cout << "\n"; 

  std::cout << "Disconnection query time is: " << t2.time() << std::endl;
  t2.reset();
  
  t1.stop();*/

  std::cout << "Total disconnection runtime is: " << t1.time() << std::endl;
  t1.reset();
  
  // use different colors, and put a few sleeps/clear.

//display segment
  CGAL::Geomview_stream gv(CGAL::Bbox_3(-6, -6, -6, 6, 6, 6));
  gv.set_line_width(2);
  // gv.set_trace(true);
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();
  
  gv << asd.as;
  std::cout << "Drawing 3d Alpha Shapes.\n";
  gv.set_wired(false);


  while(true)
  {
    std::cout << "Enter 0 to finish or anything else to query" << std::endl;
    char ch;
    std::cin >> ch;
    if(ch == '0') break;


    std::cout << "Enter two points to query disconnection" << std::endl;

    double x1,y1,z1, x2,y2,z2;

    std::cout << "First point x: " << std::endl;
    std::cin >> x1;

    std::cout << "First point y: " << std::endl;
    std::cin >> y1;

    std::cout << "First point z: " << std::endl;
    std::cin >> z1;

    std::cout << "Second point x: " << std::endl;
    std::cin >> x2;

    std::cout << "Second point y: " << std::endl;
    std::cin >> y2;

    std::cout << "Second point z: " << std::endl;
    std::cin >> z2;


    qs = Bare_point(x1,y1,z1);
    qg = Bare_point(x2,y2,z2);


    std::cout << "Truth or dare: ";
    t2.start();
    if(asd.Query_disconnection(qs,qg)){  //test disconnection
      t2.stop();
      std::cout << "Connected";}
    else{
      t2.stop();
      std::cout << "Disconnected";}
    std::cout << "\n"; 

    std::cout << "Disconnection query time is: " << t2.time() << std::endl;
    t2.reset();
    

  }


  

  return 0;
}
#endif
