// Function definitions program
#include "defs.h"
 
// Mesh class function definitions

Mesh::Mesh(int x,int y) { // constructor def. - creates the vector of x_nodes and y_nodes
  x_nodes = x;y_nodes = y;
  cout<<"x_nodes: "<<x_nodes<<"\t"<<"y_nodes: "<<y_nodes<<" in create_mesh location before list creations\n";

  for (int i=0;i<x_nodes;++i) x_list.push_back(i);

  for (int i=0;i<y_nodes;++i) y_list.push_back(i);

  cout<<"x_nodes: "<<x_nodes<<"\t"<<"y_nodes: "<<y_nodes<<" in create_mesh location after list creations\n";
}

void Mesh::print_mesh() { //print_mesh def.
  cout << "x_nodes: \n";
  for (unsigned i=0;i<x_list.size();++i) {
    cout << x_list[i] << "\n";
  }

  cout << "y_nodes: \n";
  for (unsigned i=0;i<y_list.size();++i) {
    cout << y_list[i] << "\n";
  }
}

// Point class function definitions
Point::Point(int x,int y) //constructor
  :x_coord{x},y_coord{y} {
}

void Point::print_point() {
  cout<<"["<<x_coord<<","<<y_coord<<"]";
}


// Tools class function definitions 
void Tools::point_mesh_valid(Mesh m,Point pt) { //is_valid def. - checks for validity of point
  bool validity = (m.x_nodes>=pt.x_coord && m.y_nodes>=pt.y_coord)? true:false;
  if (validity==false) {
    cout <<"Point is not in bounds: ";
    p.print_point();
    //error("Invalid Point");
  }
}
