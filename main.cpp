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
  cout<<"imax: "<<imax<<endl;
  for(int j=j0; j<jmax; j++)
    for(int i=i0; i<imax; i++) {

      v[0][j][i] = coords[0][j][i][0] + coords[0][j][i][1];

    }


  //cout<<"node(0,0):\n"<<"x-coord: "<<coords[0][0][0][0]<<endl
  //<<"y-coord: "<<coords[0][0][0][1]<<endl;  



  V.RestoreDataPointerAndInsert();
  coordinates.RestoreDataPointerToLocalVector();



  V.StoreMeshCoordinates(coordinates);
  V.WriteToVTRFile("Grid.vtr", "MySol");




  cout<<"dx = "<<dx[0]<<endl;

  //modifying xcoords and ycoords
  //xcoords.insert(xcoords.begin(),xcoords[0]-dx[0]);
  xcoords.push_back(xcoords.back()+dx[0]);

  //ycoords.insert(ycoords.begin(),ycoords[0]-dy[0]);
  ycoords.push_back(ycoords.back()+dy[0]);
  
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
    cout <<"z = "<<zcoords[v]<<endl;
  }

  cout << "Node #'s & coordinates of embedded surface:\n";
  Tools tool; 
  vector<Vec3D> Nodes;vector<Int2> Elements;
  tool.ReadMeshFileInTopFormat("embedded_arbitrary.top",Nodes,Elements);//function to extract coords. from m2c
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
/*
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

**  //max coordinats of embedded surface
  //tool.max_coords(xmax,ymax,xmin,ymin,Nodes); //TODO
  cout<<"Embedded surface x max: "<<xmax<<endl;
  cout<<"Embedded surface y max: "<<ymax<<endl;
  cout<<"Embedded surface x min: "<<xmin<<endl;
  cout<<"Embedded surface y min: "<<ymin<<endl;
**/

  // Color grid info.
  cout<<"imax: "<<imax<<endl<< "jmax: "<<jmax<<endl;
  cout<<"i0: "<<i0<<endl<< "j0: "<<j0<<endl;


  SpaceVariable3D Color(comm, &(dms.ghosted1_1dof));
  double*** color = Color.GetDataPointer();

  int start_nodex = 0; //starting at node(i=0,j=0)
  int start_nodey = 0;
  cout<<"After start nodes"<<endl;

/*hard coding a square boundary
  vector<int> i_coords{2,3,4,5,6,7}; //for top and bottom side of square
  vector<int> j_coords{2,7};
  for (auto j=0;j<j_coords.size();j++){
    for (auto i=0;i<i_coords.size();i++){
      color[0][j_coords[j]][i_coords[i]]=2;
    }
  }
//neWidth",2 <---- Future note: find out what this means!
  vector<int> j2_coords{2,3,4,5,6,7}; //for top and bottom side of square
  vector<int> i2_coords{2,7};
  for (auto i=0;i<i2_coords.size();i++){
    for (auto j=0;j<j2_coords.size();j++){
      color[0][j2_coords[j]][i2_coords[i]]=2;
    }
  }
*/
 

  cout<<"Before Flood fill."<<endl; 
  tool.flood_fill(start_nodex,start_nodey,imax,jmax,i0,j0,color,Nodes,Elements,xcoords,ycoords);
  //color[0][11][15]=3;
  cout<<"After Flood fill & before print_color"<<endl;
 // tool.SpaceVariable3D_print(color,imax,jmax); //printing the color values of each point
  cout<<"After print_color \n";
  cout<<"imax: "<<imax<<endl<< "jmax: "<<jmax<<endl;
  cout<<"i0: "<<i0<<endl<< "j0: "<<j0<<endl;
  cout<<"dx = "<<dx[0]<<endl;
  double x1=0;double y1=1.14286;double x2=0;double y2=0.914286;double x3=-1;double y3=0.9;double x4=1.1;double y4=1.1;
  double t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //1st bezier parameter
  double u = ((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //2nd bezier parameter
  cout<<"t = "<<t<<"\t"<<"u = "<<u<<endl;

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


