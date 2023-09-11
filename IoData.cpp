/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Utils.h>
#include <IoData.h>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cmath>
#include <unistd.h>
#include <bits/stdc++.h> //INT_MAX
//#include <dlfcn.h>
using namespace std;

double avogadro_number = 6.02214076e23;

RootClassAssigner *nullAssigner = new RootClassAssigner;

//------------------------------------------------------------------------------

MeshResolution1DPointData::MeshResolution1DPointData()
{
  coord = 0.0;
  h     = 0.0;
}

//------------------------------------------------------------------------------

Assigner *MeshResolution1DPointData::getAssigner()
{

  ClassAssigner *ca = new ClassAssigner("normal", 2, nullAssigner);

  new ClassDouble<MeshResolution1DPointData>
    (ca, "Coordinate", this, &MeshResolution1DPointData::coord);
  new ClassDouble<MeshResolution1DPointData>
    (ca, "CellWidth", this, &MeshResolution1DPointData::h);

  return ca;
}

//------------------------------------------------------------------------------

MeshData::MeshData()
{
  type = THREEDIMENSIONAL;
  x0 = 0.0;
  xmax = 1.0;
  y0 = 0.0;
  ymax = 1.0;
  z0 = 0.0;
  zmax = 1.0;
  Nx = -1;
  Ny = -1;
  Nz = -1;

  bc_x0   = NONE;
  bc_xmax = NONE;
  bc_y0   = NONE;
  bc_ymax = NONE;
  bc_z0   = NONE;
  bc_zmax = NONE;
}

//------------------------------------------------------------------------------

void MeshData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 20, father);

  new ClassToken<MeshData>(ca, "Type", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::type), 3,
                               "ThreeDimensional", 0, "Spherical", 1, "Cylindrical", 2);
  new ClassDouble<MeshData>(ca, "X0", this, &MeshData::x0);
  new ClassDouble<MeshData>(ca, "Xmax", this, &MeshData::xmax);
  new ClassDouble<MeshData>(ca, "Y0", this, &MeshData::y0);
  new ClassDouble<MeshData>(ca, "Ymax", this, &MeshData::ymax);
  new ClassDouble<MeshData>(ca, "Z0", this, &MeshData::z0);
  new ClassDouble<MeshData>(ca, "Zmax", this, &MeshData::zmax);
  new ClassInt<MeshData>(ca, "NumberOfCellsX", this, &MeshData::Nx);
  new ClassInt<MeshData>(ca, "NumberOfCellsY", this, &MeshData::Ny);
  new ClassInt<MeshData>(ca, "NumberOfCellsZ", this, &MeshData::Nz);

  xpoints_map.setup("ControlPointX", ca);
  ypoints_map.setup("ControlPointY", ca);
  zpoints_map.setup("ControlPointZ", ca);

  // Inside the code: Farfield0 = Farfield = Inlet, Farfield1 = Outlet
  new ClassToken<MeshData>(ca, "BoundaryConditionX0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_x0), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionXmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_xmax), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionY0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_y0), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionYmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_ymax), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionZ0", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_z0), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
  new ClassToken<MeshData>(ca, "BoundaryConditionZmax", this,
                               reinterpret_cast<int MeshData::*>(&MeshData::bc_zmax), 12,
                               "None", 0, 
                               "Inlet", 1, "Outlet", 2, //option 1
                               "Farfield0", 1, "Farfield1", 2, //option 2,
                               "Farfield", 1,//option 3
                               "Wall", 3, //slip wall
                               "SlipWall", 3, //slip wall,
                               "StickWall", 4, //no-slip wall
                               "NoSlipWall", 4, "Symmetry", 5,
                               "Overset", 6);
 } 

//------------------------------------------------------------------------------

void MeshData::check()
{
  if(type == SPHERICAL) {
    if(Ny != 1 || Nz != 1) {
      print_error("*** Error: For a domain w/ spherical symmetry, the mesh should "
                  "have 1 cell in y- and z-dirs (%d, %d).\n", Ny, Nz); 
      exit_mpi();
    }
    if(x0<0.0) {
      print_error("*** Error: For a domain w/ spherical symmetry, x0 must be nonnegative (%e).\n",
                  x0);
      exit_mpi();
    }
  }
  else if(type == CYLINDRICAL) {
    if(Nz != 1) {
      print_error("*** Error: For a domain w/ cylindrical symmetry, the mesh should "
                  "have 1 cell in z-dir (%d).\n", Nz); 
      exit_mpi();
    }
    if(y0<0.0) {
      print_error("*** Error: For a domain w/ cylindrical symmetry, y0 must be nonnegative (%e).\n",
                  y0);
      exit_mpi();
    }
  }
}

//------------------------------------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  //Should NOT call functions in Utils (e.g., print(), exit_mpi()) because the
  //M2C communicator may have not been properly set up.
  readCmdLine(argc, argv);
  readCmdFile();
}

//------------------------------------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stdout,"\033[0;31m*** Error: Input file not provided!\n\033[0m");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//------------------------------------------------------------------------------

void IoData::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stdout,"\033[0;31m*** Error: could not open \'%s\'\n\033[0m", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stdout,"\033[0;31m*** Error: command file contained parsing errors.\n\033[0m");
    exit(error);
  }
  fclose(cmdFilePtr);
}

//------------------------------------------------------------------------------
// This function is supposed to be called after creating M2C communicator. So, 
// functions in Utils can be used.
void IoData::finalize()
{
  //Check spatial domain (for spherical and cylindrical)
  mesh.check();

}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{

  mesh.setup("Mesh");

}

//------------------------------------------------------------------------------
