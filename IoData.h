/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <Vector2D.h>
#include <Vector3D.h>
#include <Utils.h>

using std::map;

/*********************************************************************
 * class IoData reads and processes the input data provided by the user
 *********************************************************************
*/
//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:
  map<int, DataType *> dataMap;

  void setup(const char *name, ClassAssigner *p) {
    SysMapObj<DataType> *smo = new SysMapObj<DataType>(name, &dataMap);
    if (p) p->addSmb(name, smo);
    else addSysSymbol(name, smo);
  }

  ~ObjectMap() {
    for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      delete it->second;
  }
};

//------------------------------------------------------------------------------

struct MeshResolution1DPointData {

  double coord;
  double h; 

  MeshResolution1DPointData();
  ~MeshResolution1DPointData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct MeshData {

  enum Type {THREEDIMENSIONAL = 0, SPHERICAL = 1, CYLINDRICAL = 2} type;
  double x0, xmax, y0, ymax, z0, zmax;
  int Nx, Ny, Nz;

  // mesh resolution info
  ObjectMap<MeshResolution1DPointData>  xpoints_map;
  ObjectMap<MeshResolution1DPointData>  ypoints_map;
  ObjectMap<MeshResolution1DPointData>  zpoints_map;

  enum BcType {NONE = 0, INLET = 1, OUTLET = 2, SLIPWALL = 3, STICKWALL = 4, SYMMETRY = 5, 
               OVERSET = 6, SIZE = 7};
  BcType bc_x0, bc_xmax, bc_y0, bc_ymax, bc_z0, bc_zmax;

  MeshData();
  ~MeshData() {} 

  void setup(const char *, ClassAssigner * = 0);

  void check(); //!< check input parameters (for spherical & cylindrical domains)
};

//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  MeshData mesh;

public:

  IoData() {}
  IoData(int, char**);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();
  void finalize();

};
#endif
