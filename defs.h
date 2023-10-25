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
  
  void point_mesh_valid(Mesh m,Point pt);  //functions  is_valid def. - checks for validity of point
  void extract_coords(vector<int>& nlist,vector<int>& xlist,vector<int>& ylist,vector<int>& zlist); 
  void ReadMeshFileInTopFormat(const char *filename, vector<Vec3D> &Xs, vector<Int2> &Es);
  double min_list(vector<double> &list); // returns the smallest number of a given list
  double point_distance(Vec3D &p1,double &p2x,double &p2y); // returns the distance between 2 points
  void closest_point(vector<Vec3D> &d,vector<Vec3D> &es,vector<double> &Gx,vector<double> &Gy); // provides the closest Grid point for each embedded surface node
  void max_min_coords(double &xmax,double &ymax,double &xmin,double &ymin,vector<Vec3D> &es); //provides the greatest x & y value of a given vector <Vec3D>
  void grid_nodes_solid_domain(double &xmax,double &ymax,double &xmin,double &ymin,double &delta_x,double &delta_y,vector<double> &xcoord,vector<double> &ycoord,vector<Vec3D> &es); //provides the grid nodes that are contained in the solid domain, including the "wetted" surface -----> ONLY WORKS FOR RECTANGLES AS OF NOW!!! 
  void grid_element_solid_domain(); //!Still need to figure out arguments for this
  void flood_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color); //performs the flood fill algorithm
  void SpaceVariable3D_print(double*** color,int &imax,int &jmax); //prints the color values of each point
  vector<Vec3D> grid_nodes(int &imax,int &jmax); //returns a vector of all the grid nodes
};
#endif
