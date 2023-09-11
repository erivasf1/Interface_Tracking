#ifndef _DEFS_H_
#define _DEFS_H_
// Header file containing declarations
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>
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
  
  int x,y; 
  Point p{x,y};
  Mesh mt{x,y};
};
#endif
