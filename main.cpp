// Main file of boundary tracking program - Erick Rivas
#include "defs.h"


void point_mesh_valid(Mesh m,Point pt); //is_valid def. - checks for validity of point
void error(string s);

int main() {
try{
  Mesh grid{3,3}; // defines object Mesh type called grid
  grid.print_mesh();

  //creation of all points in mesh
  vector<Point>mesh_pts;
  for (int i=0;i<grid.x_nodes;++i) {
    for (int j=0;j<grid.y_nodes;++j) {
      Point pt{i,j};
      pt.print_point();
      cout<<endl;
      mesh_pts.push_back(pt);
    }
  }

  Point pt1{1,2}; //Point object creations
  Point pt2{2,2};
  Point pt3{1,1};
  
  
  
  point_mesh_valid(grid,pt1); // checks if points are valid in mesh grid
  point_mesh_valid(grid,pt2);
  point_mesh_valid(grid,pt3); 
  
  cout<<endl;

  

/*
  Triangle t1{pt1,pt2,pt3};

  t1.print_triangle();*/
}
catch (runtime_error& e) {
  cerr << "Runtime Error: " << e.what() << "\n";
}
return 0;
}


// Global functions

void error(string s) {
  throw runtime_error(s);
}

