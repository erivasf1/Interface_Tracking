// Main file of boundary tracking program - Erick Rivas
#include<defs.h>
#include<SpaceOperator.h>
#include<MeshGenerator.h>
#include<IoData.h>
#include<time.h>

using std::vector;

MPI_Comm m2c_comm;


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
  vector<Vec3D> Nodes;vector<Int2> Elements;
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
    for (unsigned int j=0;j<2;j++){
      cout<<Elements[i].v[j]<<" ";
    }
    cout<<endl;
  }

  //distance between 2 points 
  double dist=tool.point_distance(Nodes[0],xcoords[0],ycoords[0]);
  cout<<"Node: "<<Nodes[0].v[0]<<","<<Nodes[0].v[1]<<endl;
  cout<<"Grid Point: "<<xcoords[0]<<","<<ycoords[0]<<endl;
  cout<<"Distance: "<<dist<<endl;

  //closest grid points to each embedded surface node
  vector<Vec3D> close_grid_points;
  tool.closest_point(close_grid_points,Nodes,xcoords,ycoords);
  cout<<"Closest Grid points to each embedded surface\n";
  for (unsigned int i=0;i<close_grid_points.size();i++){
    cout<<"Node"<<i+1<<": "<<close_grid_points[i].v[0]<<","<<close_grid_points[i].v[1]<<","<<close_grid_points[i].v[2]<<endl;
  }

/**  //max coordinats of embedded surface
  //tool.max_coords(xmax,ymax,xmin,ymin,Nodes); //TODO
  cout<<"Embedded surface x max: "<<xmax<<endl;
  cout<<"Embedded surface y max: "<<ymax<<endl;
  cout<<"Embedded surface x min: "<<xmin<<endl;
  cout<<"Embedded surface y min: "<<ymin<<endl;
**/

  // grid info.
  cout<<"imax: "<<imax<<endl<< "jmax: "<<jmax<<endl;
  cout<<"i0: "<<i0<<endl<< "j0: "<<j0<<endl;

 
  /*coords = (Vec3D***)coordinates.GetDataPointer();

  double xmax,ymax,xmin,ymin;
  tool.max_min_coords(xmax,ymax,xmin,ymin,Nodes);
  for(int j=j0;j<jmax;j++){
    for(int i=i0;i<imax;i++){
      if(xcoords[j]>=xmin && xcoords[j]<=xmax && ycoords[i]>=ymin && ycoords[i]<=ymax) //condition if grid coord. is inside of embedded surface
        color[0][j][i] = 1;
      else
        color[0][j][i] = 0;
    }
  }*/

  SpaceVariable3D Color(comm, &(dms.ghosted1_1dof));
  double*** color = Color.GetDataPointer();

  tool.flood_fill(color,imax,jmax);

  Color.RestoreDataPointerAndInsert();
  //coordinates.RestoreDataPointerToLocalVector();
  Color.WriteToVTRFile("Color.vtr", "Color");




  Color.Destroy();
  V.Destroy();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();
  MPI_Finalize();

  return 0;

}


