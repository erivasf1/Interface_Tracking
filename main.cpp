// Main file of boundary tracking program - Erick Rivas
#include "defs.cpp"



int main() {
try{
  Mesh grid{3,3}; // defines object Mesh type called grid

  Point pt1{1,2}; //Point object creations
  Point pt2{2,2};
  Point pt3{1,1};
  
  
  point_mesh_valid(grid,pt1); // checks if points are valid in mesh grid
  point_mesh_valid(grid,pt2);
  point_mesh_valid(grid,pt3); 
  
  pt1.print_point();

  Triangle t1{pt1,pt2,pt3};

}
catch (runtime_error& e) {
  cerr << "Runtime Error: " << e.what() << "\n";
}
return 0;
}
