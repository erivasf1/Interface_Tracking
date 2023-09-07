// Function definitions program
#include "defs.h"
 
// Mesh class function definitions

void Mesh::create_mesh() { //create_mesh def.
  cout<<"x_nodes: "<<x_nodes<<"\t"<<"y_nodes: "<<y_nodes<<" in create_mesh location before list creations\n";
  for (int i=0;i<x_nodes;++i) {
    x_list.push_back(i);
  }
  for (int i=0;i<y_nodes;++i) {
    y_list.push_back(i);
  }
  cout<<"x_nodes: "<<x_nodes<<"\t"<<"y_nodes: "<<y_nodes<<" in create_mesh location after list creations\n";
}
void Mesh::print_mesh() { //print_mesh def.
  cout << "x_nodes: \n";
  for (int i=0;i<x_list.size();++i) {
    cout << x_list[i] << "\n";
  }

  cout << "y_nodes: \n";
  for (int i=0;i<y_list.size();++i) {
    cout << y_list[i] << "\n";
  }
}

// Point class function definitions
Point::Point(int x,int y) //constructor
  :x_coord{x},y_coord{y} {
}

void Point::print_point() {
  cout<<"["<<x_coord<<","<<y_coord<<"]\n";
}

//Triangle class function definitions
void Triangle::print_triangle() {
  cout<<"["<<p1.print_point()<<","<<p2.print_point()<<","<<p3.print_point()<<"]\n";
}

// Global functions
void point_mesh_valid(Mesh m,Point pt) { //is_valid def. - checks for validity of point
  bool validity = (m.x_nodes>=pt.x_coord && m.y_nodes>=pt.y_coord)? true:false;
  if (validity==false) {
    cout <<"Point is not in bounds: ";
    pt.print_point();
    error("Invalid Point");
  }
}
