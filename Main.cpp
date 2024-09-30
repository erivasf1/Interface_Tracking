// Main file of Interface Tracker - Erick Rivas
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

  //! Setting up SpaceVariables that store Color and Coordinates
  SpaceVariable3D Color(comm, &(dms.ghosted1_1dof));
  SpaceVariable3D& coordinates(spo.GetMeshCoordinates());

  double*** color = Color.GetDataPointer();
  Vec3D*** coords = (Vec3D***)coordinates.GetDataPointer();

  // This part isn't needed, it is only an example to set an arbitrary value to Color (for visualization)
  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  coordinates.GetCornerIndices(&i0, &j0, &k0, &imax, &jmax, &kmax);

  //cout<<"node(0,0):\n"<<"x-coord: "<<coords[0][0][0][0]<<endl
  //<<"y-coord: "<<coords[0][0][0][1]<<endl;  



  //Color.RestoreDataPointerAndInsert();
  //coordinates.RestoreDataPointerToLocalVector();


  //Color.StoreMeshCoordinates(coordinates);
  //Color.WriteToVTRFile("Grid.vtr", "MySol");


  /*//printf("Shape number is: %d\n",iod.shapelist.shapes.dataMap.size());
  printf("Shape1 is: %s\n",iod.shapelist.shapes.dataMap[0]->filename);
  const char* Shape1 = iod.shapelist.shapes.dataMap[0]->filename;
  printf("Shape1 is: %s\n",Shape1);
*/

  //! Storing nodes and connectivities of embedded surfaces (if any)
  Tools tool; 
  int shapes_number = iod.shapelist.shapes.dataMap.size();
  if (shapes_number < 1){
    print("No Embedded Surface Input!\n");
    exit_mpi();
  }
  else {
    vector<vector<Vec3D>> Shape_Nodes(shapes_number);
    vector<vector<Int2>> Shape_Elements(shapes_number);
    for (int n = 0; n<shapes_number; n++){

      tool.ReadMeshFileInTopFormat(iod.shapelist.shapes.dataMap[n]->filename,Shape_Nodes[n],Shape_Elements[n]);
      
    }
  
//-----------------------------SOLID-FLUID INTERFACE OF SHAPES-----------------------------

  vector<Vec3D> intersecting_nodes; //stores the nodes that pertain to the solid-fluid interface
  vector<Int2> intersecting_edges; //stores the edges/connectivities of the nodes that pertain to the solid-fluid interface
  

  // Calculating the nodes and connectivites of the interface and marking the nodes to indicate the solid-fluid interface of both solids
  int start_nodex = 0; //starting at node(i=0,j=0)
  int start_nodey = 0;

  // Creating list containing color values of all shapes involved
  vector<int> interface_colors(shapes_number);
  for (int c = 0; c<shapes_number; c++) interface_colors[c] = c+2;
  

  int shape1_color = 2; int shape2_color = 3; int overlap_color = 4; //color values of each sub-domain/topology

  // filling all solid-fluid interface
  for (int n = 0; n<shapes_number; n++) {
  
  tool.intersect_fill(start_nodex,start_nodey,imax,jmax,i0,j0,color,interface_colors[n],Shape_Nodes[n],Shape_Elements[n],xcoords,ycoords,intersecting_nodes,intersecting_edges); //intersect fill of shape 1

  }

 /*//----------------------------TOPOLOGICAL CHANGE TEST----------------------------------------
 // TODO: Change to only using flood fill 
  Topology top(Shape_Elements[0],Shape_Nodes[0],Shape_Elements[1],Shape_Nodes[1]); 

  // Calculating inside points
  vector<Vec3D> pts_inside;
  pts_inside = top.inside_points(); //vector of inside points
  cout << "inside points: "<<pts_inside.size()<<endl; //TODO: would be a good idea to have this as a function -- e.g. disp_points()
  for (int i=0;i<pts_inside.size();i++){
    cout<<"Point "<<i+1<<": "<<pts_inside[i].v[0]<<"\t"<<pts_inside[i].v[1]<<endl;
  }

  // Calculating intersecting points
  vector<Vec3D> intersecting_pts;
  intersecting_pts = top.intersecting_points(); //vector of inside points
  cout << "intersecting points: "<<intersecting_pts.size()<<endl; //TODO: would be a good idea to have this as a function -- e.g. disp_points()

  for (int i=0;i<intersecting_pts.size();i++){
    cout<<"Point "<<i+1<<": "<<intersecting_pts[i].v[0]<<"\t"<<intersecting_pts[i].v[1]<<endl;
  }

  // Constructing nodes and connectivities of overlap shape
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

  // Filling of grid nodes for 2 arbitrary shapes
  top.shape_fill(Shape1_Elements,Shape1_Nodes,xcoords,ycoords,color,shape1_color,imax,jmax,i0,j0);
  top.shape_fill(Shape2_Elements,Shape2_Nodes,xcoords,ycoords,color,shape2_color,imax,jmax,i0,j0);
  
  double x=21.4313; double y = 43.0572;
  cout<<"\nMANUAL COLOR TEST"<<endl;
  cout<<"xcoord of test node: "<<x<<endl;
  cout<<"ycoord of test node: "<<y<<endl;
  cout<<"Is inside function point check"<<endl;
  Vec3D test_pt{x,y,0};
  //color[0][
  if (top.is_inside(test_pt,Shape2_Elements,Shape2_Nodes)==true) cout<<"Point is inside!"<<endl; 
  */

  /*// Writing constructed shape (if detected) to file in .top format used for visualization - TODO: (would be a good idea to have this as a function)
  if (Connectivities.size() != 0 & Nodes.size() != 0){ //case if shapes are intersecting
    cout<<"Overlap detected"<<endl;
    ofstream myfile;
    myfile.open("overlap_shape.top");
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
 
    top.shape_fill(Connectivities,Nodes,xcoords,ycoords,color,overlap_color,imax,jmax,i0,j0); //fills the grid nodes that are in overlap shape
  }*/

  //---------------------------GRID FLOOD FILL-------------------------------------------
  
  // filling the "untouched" nodes of the grid to mark the fluid domain
  tool.flood_fill(start_nodex,start_nodey,imax,jmax,i0,j0,color);

  Color.RestoreDataPointerAndInsert();
  Color.WriteToVTRFile("Color.vtr", "Color");

  coordinates.RestoreDataPointerToLocalVector();


// Restoring Memory to heap
  Color.Destroy();
  //V.Destroy();
  spo.Destroy();
  dms.DestroyAllDataManagers();
  PetscFinalize();
  MPI_Finalize();
}
  return 0;

}


