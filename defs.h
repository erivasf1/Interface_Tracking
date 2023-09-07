// Header file containing declarations
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <algorithm>
using namespace std;

void error(string s) {
  throw runtime_error(s);
}

struct Mesh { // Mesh struct declarations -- Base Class

  Mesh(int x,int y){ // constructor
  x_nodes = x;y_nodes = y;
  }

  void create_mesh(); // functions
  void print_mesh();
  
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
  Triangle(Point pt1,Point pt2,Point pt3) { //constructor definition  
   p1=pt1;p2=pt2;p3=pt3;
   tri.push_back(pt1); 
   tri.push_back(pt2); 
   tri.push_back(pt3); 
  }
  void print_triangle(); //functions
  
private: 
  int x,y;
  Point p1{x,y};Point p2{x,y};Point p3{x,y};
  vector<Point> tri;
};


