// File to Move User-defined Shape - Erick Rivas
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

//------------------------PARAMETERS-----------------------
  Tools tool;
  vector<Vec3D> Nodes;
  vector<Int2> Elements;
  double dx = 10.0; //desired change for nodes in x-dir.
  double dy = 15.0; //desired change for nodes in y-dir.
  double* dx_ptr = &dx; double* dy_ptr = &dy;

//-------------------READING ORIGINAL FILE--------------------
  tool.ReadMeshFileInTopFormat("small_sub_slice.top",Nodes,Elements);

//-----------------NEW POSITIONS FOR NODES--------------------
  for (auto i=0;i<Nodes.size();i++){
    Nodes[i].v[0] += *dx_ptr; //x-dir. update
    Nodes[i].v[1] += *dy_ptr; //y-dir. update
  }

//---------------WRITING TO NEW FILE--------------------------
  ofstream myfile;
  myfile.open("Moved_Shape.top");
  if (myfile.is_open()){

    myfile<<"Nodes MySurfaceNodes"<<endl;
    for (int i=0;i<Nodes.size();i++){
      myfile<<i+1<<"  "<<Nodes[i].v[0]<<"  "<<Nodes[i].v[1]<<"  "<<Nodes[i].v[2]<<endl;
    }


    myfile<<"Elements MySurface using MySurfaceNodes"<<endl;
    for (int i=0;i<Elements.size();i++){
      myfile<<i+1<<"  "<<1<<"  "<<Elements[i].v[0]+1<<"  "<<Elements[i].v[1]+1<<endl;
    }
    myfile.close();


  }
  return 0;
}
