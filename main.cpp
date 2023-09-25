// Main file of boundary tracking program - Erick Rivas
#include<defs.h>
#include<SpaceOperator.h>
#include<MeshGenerator.h>
#include<IoData.h>
#include<time.h>

using std::vector;

MPI_Comm m2c_comm;

//void point_mesh_valid(Mesh m,Point pt); //is_valid def. - checks for validity of point
void error(string s);

int main(int argc, char* argv[]) {

  // -----------------------------------------------------------------
  // Initialize the environment
  [[maybe_unused]] clock_t start_time = clock(); //for timing purpose only
  //! Initialize MPI 
  MPI_Init(NULL,NULL); //called together with all concurrent programs -> MPI_COMM_WORLD
  //! Print header (global proc #0, assumed to be a M2C proc)
  m2c_comm = MPI_COMM_WORLD; //temporary, just for the next few lines of code
  //! Read user's input file (read the parameters)
  IoData iod(argc, argv);
  //! Partition MPI, if there are concurrent programs
  MPI_Comm comm = m2c_comm; //this is going to be the M2C communicator
  //! Finalize IoData (read additional files and check for errors)
  iod.finalize();
  //! Initialize PETSc
  PETSC_COMM_WORLD = comm;
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0);
  // -----------------------------------------------------------------

  // -----------------------------------------------------------------
  vector<double> xcoords, dx, ycoords, dy, zcoords, dz;
  MeshGenerator meshgen;
  meshgen.ComputeMeshCoordinatesAndDeltas(iod.mesh, xcoords, ycoords, zcoords, dx, dy, dz);
  // -----------------------------------------------------------------

  //! Setup PETSc data array (da) structure for nodal variables
  DataManagers3D dms(comm, xcoords.size(), ycoords.size(), zcoords.size());


  SpaceOperator spo(comm, dms, iod, xcoords, ycoords, zcoords, dx, dy, dz);



  // ----------------------------------------------
  //  START
  // ----------------------------------------------

  SpaceVariable3D V(comm, &(dms.ghosted1_1dof));
  SpaceVariable3D& coordinates(spo.GetMeshCoordinates());

  double*** v = V.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  
  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);
  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++) {

      v[0][j][i] = coords[0][j][i][0] + coords[0][j][i][1];

    }






  V.RestoreDataPointerAndInsert();
  coordinates.RestoreDataPointerToLocalVector();



  V.StoreMeshCoordinates(coordinates);
  V.WriteToVTRFile("Grid.vtr", "MySol");






try{
  
  cout << "X coords. of grid:\n"; //prints out the x coords. of the grid
  for (unsigned int v=0;v<xcoords.size();v++) {
    cout <<"x = "<<xcoords[v]<<endl;
  }
  cout << "Y coords. of grid:\n"; //prints out the y coords. of the grid
  for (unsigned int v=0;v<ycoords.size();v++) {
    cout <<"y = "<<ycoords[v]<<endl;
  }
  cout << "Z coords. of grid:\n"; //prints out the z coords. of the grid
  for (unsigned int v=0;v<zcoords.size();v++) {
    cout <<"x = "<<zcoords[v]<<endl;
  }

  cout << "Node #'s & coordinates of embedded surface:\n";
  Tools tool;  
  vector<int> node_list;
  vector<int> x_list;
  vector<int> y_list;
  vector<int> z_list;
  tool.extract_coords(node_list,x_list,y_list,z_list); // function call to extract coords. from .top file
  cout<<"Node #s\n";
  for (unsigned int i=0;i<node_list.size();i++){
    cout<<node_list[i]<<endl;
  }
  cout<<"X coords\n";
  for (unsigned int i=0;i<x_list.size();i++){
    cout<<x_list[i]<<endl;
  }
  cout<<"Y coords\n";
  for (unsigned int i=0;i<y_list.size();i++){
    cout<<y_list[i]<<endl;
  }
  cout<<"Z coords\n";
  for (unsigned int i=0;i<z_list.size();i++){
    cout<<z_list[i]<<endl;
  }  
  
  vector<Vec3D> Nodes;vector<Int3> Elements;
  tool.ReadMeshFileInTopFormat("embedded_surface1.top",Nodes,Elements);//function to extract coords. from m2c
  cout<<"Nodes of embedded surface:\n";
  for (unsigned int i=0;i<Nodes.size();i++){
    cout<<"Point: ";
    for (unsigned int j=0;j<3;j++){
      cout<<Nodes[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<"Elements of embedded surface:\n";
  for (unsigned int i=0;i<Elements.size();i++){
    cout<<"Element: ";
    for (unsigned int j=0;j<3;j++){
      cout<<Elements[i].v[j]<<" ";
    }
    cout<<endl;
  }

/*
  Point pt1{1,2}; //Point object creations
  Point pt2{2,2};
  Point pt3{1,1};
  
  
  Tools tools; 
  tools.point_mesh_valid(grid,pt1); // checks if points are valid in mesh grid
  tools.point_mesh_valid(grid,pt2);
  tools.point_mesh_valid(grid,pt3); 
  
  cout<<endl;

  


  Triangle t1{pt1,pt2,pt3};

  t1.print_triangle();*/
}
catch (runtime_error& e) {
  cerr << "Runtime Error: " << e.what() << "\n";
}







  V.Destroy();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();
  MPI_Finalize();

  return 0;
}


// Global functions

void error(string s) {
  throw runtime_error(s);
}

