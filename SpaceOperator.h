/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _SPACEOPERATOR_H_
#define _SPACEOPERATOR_H_
#include<SpaceVariable.h>
#include<GhostPoint.h>
#include<IoData.h>

/*******************************************
 * class SpaceOperator drives computations
 * that require domain/mesh information
 ******************************************/
class SpaceOperator
{
  MPI_Comm&                 comm;
  DataManagers3D&           dm_all;
  IoData&                   iod;

  //! Mesh info
  SpaceVariable3D coordinates;
  SpaceVariable3D delta_xyz;
  SpaceVariable3D volume; //!< volume of node-centered control volumes
  
  SpaceVariable3D Tag;

  vector<GhostPoint> ghost_nodes_inner; //!< ghost nodes inside the physical domain (shared with other subd)
  vector<GhostPoint> ghost_nodes_outer; //!< ghost nodes outside the physical domain

  int i0, j0, k0, imax, jmax, kmax; //!< corners of the real subdomain
  int ii0, jj0, kk0, iimax, jjmax, kkmax; //!< corners of the ghosted subdomain
  int NX, NY, NZ; //!< global size

public:
  SpaceOperator(MPI_Comm &comm_, DataManagers3D &dm_all_, IoData &iod_, 
                vector<double> &x, vector<double> &y, vector<double> &z,
                vector<double> &dx, vector<double> &dy, vector<double> &dz);
  ~SpaceOperator();

  //! Reset the coords of ghost layer nodes  (a NULL pointer means that value does not need to be reset)
  void ResetGhostLayer(double* xminus, double* xplus, double* yminus,  double* yplus,
                       double* zminus, double* zplus, double* dxminus, double* dxplus, double* dyminus,
                       double* dyplus, double* dzminus, double* dzplus);

  SpaceVariable3D& GetMeshCoordinates() {return coordinates;}
  SpaceVariable3D& GetMeshDeltaXYZ()    {return delta_xyz;}
  SpaceVariable3D& GetMeshCellVolumes() {return volume;}

  vector<GhostPoint>* GetPointerToInnerGhostNodes() {return &ghost_nodes_inner;}
  vector<GhostPoint>* GetPointerToOuterGhostNodes() {return &ghost_nodes_outer;}

  void Destroy();


private:

  void SetupMesh(vector<double> &x, vector<double> &y, vector<double> &z,
                 vector<double> &dx, vector<double> &dy, vector<double> &dz);
  void PopulateGhostBoundaryCoordinates();

  void CreateGhostNodeLists(bool screenout);

};


#endif
