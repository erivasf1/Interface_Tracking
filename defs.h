#ifndef _DEFS_H_
#define _DEFS_H_
// Header file containing declarations
#include <map>
#include <Utils.h>
#include <cassert>
#include <deque>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <Vector3D.h>
#include <Vector2D.h>
using namespace std;

struct Mesh { // Mesh struct declarations -- Base Class

  Mesh(int x,int y); // constructor

  void print_mesh(); // functions`
  
  int x_nodes,y_nodes; // variables
  vector<int> x_list;
  vector<int> y_list;
};

struct Point { // Point struct declarations
  
  Point(int x,int y); //constructor
  
  
  void print_point(); // functions

  int x_coord,y_coord; // variables

};

struct Triangle { // Triangle struct declarations - derived from Base class - Point
  Triangle(Point pt1,Point pt2,Point pt3) :p1{pt1},p2{pt2},p3{pt3} { //constructor definition  
   tri.push_back(p1); 
   tri.push_back(p2); 
   tri.push_back(p3); 
  }
  void print_triangle(); //functions
  
private: 
  int x,y;
  Point p1{x,y};Point p2{x,y};Point p3{x,y};
  vector<Point> tri;
};

struct Tools { //Tools used to compare different user-defined types
  
  // functions for Tool
  void point_mesh_valid(Mesh m,Point pt);  //functions  is_valid def. - checks for validity of point
  void extract_coords(vector<int>& nlist,vector<int>& xlist,vector<int>& ylist,vector<int>& zlist); 
  void ReadMeshFileInTopFormat(const char *filename, vector<Vec3D> &Xs, vector<Int2> &Es);
  double min_list(vector<double> &list); // returns the smallest number of a given list
  double point_distance(Vec3D &pt1,Vec3D &pt2); // returns the distance between 2 points
  Vec3D closest_point(Vec3D &pt1,vector<Vec3D> &pts); // provides the closest Grid point for each embedded surface node
  void max_min_coords(double &xmax,double &ymax,double &xmin,double &ymin,vector<Vec3D> &es); //provides the greatest x & y value of a given vector <Vec3D>
  int getIndex(vector<Vec3D> &v,Vec3D k); //returns the index value of a Vec3D in a Vec3D list
  vector<Vec3D> remove(vector<Vec3D> &v,Vec3D &N); //returns the vector of the desired deleted Vec3D element
  void grid_nodes_solid_domain(double &xmax,double &ymax,double &xmin,double &ymin,double &delta_x,double &delta_y,vector<double> &xcoord,vector<double> &ycoord,vector<Vec3D> &es); //provides the grid nodes that are contained in the solid domain, including the "wetted" surface -----> ONLY WORKS FOR RECTANGLES AS OF NOW!!! 
  void grid_element_solid_domain(); //!Still need to figure out arguments for this
  void flood_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color); //performs the flood fill algorithm
  void intersect_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color,int &color_val,vector<Vec3D> &surface_nodes,vector<Int2> &surface_connectivities,vector<double> &xcoords,vector<double> &ycoords,vector<Vec3D> &intersecting_nodes,vector<Int2> &intersecting_edges);
  bool intersect(int grid_node1i,int grid_node1j,int grid_node2i,int grid_node2j,int &imax,int &jmax,int &imin,int &jmin,vector<Vec3D> &surface_nodes,vector<Int2> &surface_connectivities,vector<double> &xcoords,vector<double> &ycoords); //returns true if an interesection is detected

  double point_of_intersection(Vec3D &grid_node1,Vec3D &grid_node2,Vec3D &surface_node1,Vec3D &surface_node2); //bezier-parameter calc. to determine intersection
  void SpaceVariable3D_print(double*** color,int &imax,int &jmax); //prints the color values of each point
};

class Topology:Tools{ //Topology is used to determine the overlapping area of 2 shapes -- inherits the functions of the Tool class
  
  //connectivities & nodes of Shape 1
  vector<Int2> shape1_elements; vector<Vec3D> shape1_nodes;
  //connectivities & nodes of Shape 2 
  vector<Int2> shape2_elements; vector<Vec3D> shape2_nodes;

public:
  //constructor
  Topology(vector<Int2> connectivities1,vector<Vec3D> nodes1,vector<Int2> connectivities2,vector<Vec3D> nodes2){
    shape1_elements = connectivities1; shape1_nodes = nodes1;  
    shape2_elements = connectivities2; shape2_nodes = nodes2;  
  }
  //functions

  bool is_inside(Vec3D &point,vector<Int2> &elements,vector<Vec3D> &nodes);// returns true if point is inside shape
  Vec3D pt_of_intersection(Vec3D &pt1_A,Vec3D &pt2_A,Vec3D &pt1_B,Vec3D &pt2_B); // returns the intersecting point between 2 line segments
  vector<Vec3D> inside_points(); // returns the inside points of a shape into another shape , also vice versa
  vector<Vec3D> intersecting_points(); // returns the coordinates of the intersecting nodes of the overlaped shape
  void shape_construct(vector<Int2> &overlap_connectivities,vector<Vec3D> &overlap_nodes,vector<Vec3D> &intersecting_pts,vector<Vec3D> &inside_pts); // generates the connectivites and nodes of the overlap shape
  void shape_fill(vector<Int2> &overlap_connectivities,vector<Vec3D> &overlap_nodes,vector<double> &xcoords,vector<double> &ycoords,double*** &color, int color_val); // fills all the grid nodes inside the overlap shape with a color
  void print_overlap_shape(); //prints the connectivities and nodes of the overlapped shape
};
#endif
