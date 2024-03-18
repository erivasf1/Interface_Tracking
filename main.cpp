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
  PetscInitialize(&argc, &argv, argc>=3 ? argv[2] : (char*)0, (char*)0); // -----------------------------------------------------------------

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
  
  /*cout << "X coords. of grid:\n"; //prints out the x coords. of the grid
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
*/
  Tools tool; 
  vector<Vec3D> Nodes_rect1;vector<Int2> Elements_rect1; //nodes & connectivities of sub
  vector<Vec3D> Nodes_rect2;vector<Int2> Elements_rect2; //nodes & connectivities of arb. shape
  
  tool.ReadMeshFileInTopFormat("rect1.top",Nodes_rect1,Elements_rect1);//function to extract coords. from m2c
  tool.ReadMeshFileInTopFormat("rect2.top",Nodes_rect2,Elements_rect2);//function to extract coords. from m2c
  // Color grid info.
  cout<<"imax: "<<imax<<endl<< "jmax: "<<jmax<<endl;
  cout<<"i0: "<<i0<<endl<< "j0: "<<j0<<endl;


  SpaceVariable3D Color(comm, &(dms.ghosted1_1dof));
  double*** color = Color.GetDataPointer();

  vector<Vec3D> intersecting_nodes; //for recording intersecting connectivities and nodes
  vector<Int2> intersecting_edges;

  int start_nodex = 0; //starting at node(i=0,j=0)
  int start_nodey = 0;

  int shape1_color = 2; int shape2_color = 3; int overlap_color = 4;

  tool.intersect_fill(start_nodex,start_nodey,imax,jmax,i0,j0,color,shape1_color,Nodes_rect1,Elements_rect1,xcoords,ycoords,intersecting_nodes,intersecting_edges); //intersect fill of sub

  tool.intersect_fill(start_nodex,start_nodey,imax,jmax,i0,j0,color,shape2_color,Nodes_rect2,Elements_rect2,xcoords,ycoords,intersecting_nodes,intersecting_edges); //intersect fill of arb.shape

 //----------------------------TOPOLOGICAL CHANGE TEST----------------------------------------
 
  Topology top(Elements_rect1,Nodes_rect1,Elements_rect2,Nodes_rect2); 

  //INSIDE POINTS
  vector<Vec3D> pts_inside;
  pts_inside = top.inside_points(); //vector of inside points
  cout << "inside points: "<<pts_inside.size()<<endl; //would be a good idea to have this as a function -- e.g. disp_points()
  for (int i=0;i<pts_inside.size();i++){
    cout<<"Point "<<i+1<<": "<<pts_inside[i].v[0]<<"\t"<<pts_inside[i].v[1]<<endl;
  }

  //INTERSECTING POINTS
  vector<Vec3D> intersecting_pts;
  intersecting_pts = top.intersecting_points(); //vector of inside points
  cout << "intersecting points: "<<intersecting_pts.size()<<endl; //would be a good idea to have this as a function -- e.g. disp_points()

  for (int i=0;i<intersecting_pts.size();i++){
    cout<<"Point "<<i+1<<": "<<intersecting_pts[i].v[0]<<"\t"<<intersecting_pts[i].v[1]<<endl;
  }

  //CONSTRUCTION OF NODES AND CONNECTIVITIES OF OVERLAP SHAPE
  vector<Vec3D> Nodes; vector<Int2> Connectivities;

  top.shape_construct(Connectivities,Nodes,intersecting_pts,pts_inside);
  cout<<"Nodes after overlap shape construction"<<endl;
  for (int i=0;i<Nodes.size();i++){
    cout<<"Node "<<i+1<<": "<<Nodes[i].v[0]<<"\t"<<Nodes[i].v[1]<<endl;
  }

  cout<<"Connectivities after overlap shape construction"<<endl;
  for (int i=0;i<Connectivities.size();i++){
    cout<<"Element "<<i+1<<": "<<Connectivities[i].v[0]<<"\t"<<Connectivities[i].v[1]<<endl;
  }

  //FILLING OF GRID NODES FOR 2 ARBITRARY SHAPES
  top.shape_fill(Elements_rect1,Nodes_rect1,xcoords,ycoords,color,shape1_color);
  top.shape_fill(Elements_rect2,Nodes_rect2,xcoords,ycoords,color,shape2_color);

  //WRITING CONSTRUCTED SHAPE TO FILE IN .TOP FORMAT - (would be a good idea to have this as a function)
  if (Connectivities.size() != 0 & Nodes.size() != 0){ //case if shapes are intersecting
    cout<<"Overlap detected"<<endl;
    ofstream myfile;
    myfile.open("overlap_shape.txt");
    if (myfile.is_open()){

      myfile<<"Nodes MySurfaceNodes"<<endl;
      for (int i=0;i<Nodes.size();i++){
        myfile<<i+1<<"  "<<Nodes[i].v[0]<<"  "<<Nodes[i].v[1]<<"  "<<Nodes[i].v[2]<<endl;
      }


      myfile<<"Elements MySurface using MySurfaceNodes"<<endl;
      for (int i=0;i<Connectivities.size();i++){
        myfile<<i+1<<"  "<<1<<"  "<<Connectivities[i].v[0]+1<<"  "<<Connectivities[i].v[1]+1<<endl;
      }
      myfile.close();

    }
 
    top.shape_fill(Connectivities,Nodes,xcoords,ycoords,color,overlap_color); //fills the grid nodes that are in overlap shape
  }

  //---------------------------GRID FLOOD FILL-------------------------------------------


  tool.flood_fill(start_nodex,start_nodey,imax,jmax,i0,j0,color);

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


