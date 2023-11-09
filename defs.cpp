// Function definitions program
#include "defs.h"
 
//----Mesh class function definitions--------------------

Mesh::Mesh(int x,int y) { // constructor def. - creates the vector of x_nodes and y_nodes
  x_nodes = x;y_nodes = y;
}

void Mesh::print_mesh() { //print_mesh def.
  cout << "x_nodes: \n";
  for (unsigned i=0;i<x_list.size();++i) {
    cout << x_list[i] << "\n";
  }

  cout << "y_nodes: \n";
  for (unsigned i=0;i<y_list.size();++i) {
    cout << y_list[i] << "\n";
  }
}

//----Point class function definitions----------------------
Point::Point(int x,int y) //constructor
  :x_coord{x},y_coord{y} {
}

void Point::print_point() {
  cout<<"["<<x_coord<<","<<y_coord<<"]";
}


//----Tools class function definitions--------------------- 
void Tools::point_mesh_valid(Mesh m,Point pt) { //is_valid def. - checks for validity of point
/*
  bool validity = (m.x_nodes>=pt.x_coord && m.y_nodes>=pt.y_coord)? true:false;
  if (validity==false) {
    cout <<"Point is not in bounds: ";
    p.print_point();
    //error("Invalid Point");
  }
*/
}
void Tools::extract_coords(vector<int>& nlist,vector<int>& xlist,vector<int>& ylist,vector<int>& zlist){

  string file = "embedded_surface1.top";

  ifstream ist{file}; // puts all contents of file into ist stream

  char x;
  int i=1;
  while (ist>>x){ // starts reading chars of ist one by one into x
    switch(x){
    case'E':
    {
      ist.putback(x);
      string word;
      ist>>word;
      if (word=="Elements")return;
    }
    case'0':case'1':case'2':case'3':case'4':case'5':case'6':case'7':case'8':case'9':
    case'-':
    {
      switch(i){
      case 1:{
      ist.putback(x); //puts the char x back into the ist stream for int read prep
      int num;ist>>num; //int values will be read into num now
      nlist.push_back(num); //appends the read num into the num array
      i++;break;
      }
      case 2:{
      ist.putback(x); //puts the char x back into the ist stream for int read prep
      int num;ist>>num; //int values will be read into num now
      xlist.push_back(num); //appends the read num into the num array
      i++;break;}
      case 3:{
      ist.putback(x); //outs the char x back into the ist stream for int read prep
      int num;ist>>num; //int values will be read into num now
      ylist.push_back(num); //appends the read num into the num array
      i++;break;}
      case 4:{
      ist.putback(x); //outs the char x back into the ist stream for int read prep
      int num;ist>>num; //int values will be read into num now
      zlist.push_back(num); //appends the read num into the num array
      i=1;break;}
      }

    }
    }
  }
} 
void Tools::ReadMeshFileInTopFormat(const char *filename, vector<Vec3D> &Xs, vector<Int2> &Es){

  // read data from the surface input file.
  FILE *topFile;
  topFile = fopen(filename, "r");
  if(topFile == NULL) {
    print_error("*** Error: Cannot open embedded surface mesh file (%s).\n", filename);
    exit_mpi();
  }
 
  int MAXLINE = 500;
  char line[MAXLINE], key1[MAXLINE], key2[MAXLINE]; //, copyForType[MAXLINE];

  int num0 = 0;
  int num1 = 0;
  double x1, x2, x3;
  int node1, node2;
  int type_reading = 0; //1 means reading node set; 2 means reading element set
  std::deque<std::pair<int, Vec3D>> nodeList;
  std::deque<std::array<int, 3>> elemList; // element ID + three node IDs
  int maxNode = 0, maxElem = 0;
  bool found_nodes = false;
  bool found_elems = false;


  // --------------------
  // Read the file
  // --------------------
  while(fgets(line, MAXLINE, topFile) != 0) {

    sscanf(line, "%s", key1);
    string key1_string(key1);

    if(key1[0] == '#') {
      //Do nothing. This is user's comment
    }
    else if(same_strings_insensitive(key1_string,"Nodes")){
      if(found_nodes) {//already found nodes... This is a conflict
        print_error("*** Error: Found multiple sets of nodes (keyword 'Nodes') in %s.\n", filename);
        exit_mpi();
      }
      sscanf(line, "%*s %s", key2);
      type_reading = 1;
      found_nodes = true;
    }
    else if(same_strings_insensitive(key1_string, "Elements")) {

      if(found_elems) {//already found elements... This is a conflict
        print_error("*** Error: Found multiple sets of elements (keyword 'Elements') in %s.\n", filename);
        exit_mpi();
      }
      type_reading = 2;
      found_elems = true;

    }
    else if(type_reading == 1) { //reading a node (following "Nodes Blabla")
      int count = sscanf(line, "%d %lf %lf %lf", &num1, &x1, &x2, &x3);
      if(count != 4) {
        print_error("*** Error: Cannot interpret line %s (in %s). Expecting a node.\n", line, filename);
        exit_mpi();
      }
      if(num1 < 1) {
        print_error("*** Error: detected a node with index %d in embedded surface file %s.\n", num1, filename);
        exit_mpi();
      }
      if(num1 > maxNode)
        maxNode = num1;

      nodeList.push_back({num1, {x1, x2, x3}});
    }
    else if(type_reading == 2) { // we are reading an element --- HAS TO BE A LINE SEGMENT!
      int count = sscanf(line, "%d %d %d %d", &num0, &num1, &node1, &node2);
      if(count != 4) {
        print_error("*** Error: Cannot interpret line %s (in %s). Expecting a line segment.\n", line, filename);
        exit_mpi();
      }
      if(num0 < 1) {
        print_error("*** Error: detected an element with index %d in embedded surface file %s.\n", num0, filename);
        exit_mpi();
      }
      if(num0 > maxElem)
        maxElem = num0;

      elemList.push_back({num0, node1, node2});
    }
    else { // found something I cannot understand...
      print_error("*** Error: Unable to interpret line %s (in %s).\n", line, filename);
      exit_mpi();
    }

  }

  fclose(topFile);

  if(!found_nodes) {
    print_error("*** Error: Unable to find node set in %s.\n", filename);
    exit_mpi();
  }
  if(!found_elems) {
    print_error("*** Error: Unable to find element set in %s.\n", filename);
    exit_mpi();
  }

  // ----------------------------
  // Now, check and store nodes
  // ----------------------------
  int nNodes = nodeList.size();
  map<int,int> old2new;
  Xs.resize(nNodes);
  int id(-1);
  if(nNodes != maxNode) { // need to renumber nodes, i.e. create "old2new"
    print_warning("Warning: The node indices of an embedded surface may have a gap: "
                  "max index = %d, number of nodes = %d. Renumbering nodes. (%s)\n",
                  maxNode, nNodes, filename);
//    assert(nNodes < maxNode);

    int current_id = 0; 
    vector<bool> nodecheck(maxNode+1, false);
    for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
      id = it1->first;
      if(nodecheck[id]) {
        print_error("*** Error: Found duplicate node (id: %d) in embedded surface file %s.\n", id, filename);
        exit(-1);
      }
      nodecheck[id] = true;
      Xs[current_id] = it1->second; 
      old2new[id] = current_id;
      current_id++;
    }
    assert(current_id==(int)Xs.size());
  } 
  else { //in good shape
    vector<bool> nodecheck(nNodes, false);
    for(auto it1 = nodeList.begin(); it1 != nodeList.end(); it1++) {
      id = it1->first - 1; 
      if(nodecheck[id]) {
        print_error("*** Error: Found duplicate node (id: %d) in embedded surface file %s.\n", id+1, filename);
        exit(-1);
      }
      nodecheck[id] = true;
      Xs[it1->first - 1] = it1->second;
    }
  }


  // ------------------------------
  // check nodes used by elements
  // ------------------------------
  for(auto it = elemList.begin(); it != elemList.end(); it++) {

    id = (*it)[0];
    node1 = (*it)[1];
    node2 = (*it)[2];
      
    if(old2new.empty()) {//node set is original order

      if(node1<=0 || node1 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node1, id, filename);
        exit_mpi();
      }

      if(node2<=0 || node2 > nNodes) {
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node2, id, filename);
        exit_mpi();
      }

    }
    else {// nodes are renumbered

      auto p1 = old2new.find(node1);
      if(p1 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node1, id, filename);
        exit_mpi();
      }

      auto p2 = old2new.find(node2);
      if(p2 == old2new.end()) { 
        print_error("*** Error: Detected unknown node number (%d) in element %d (%s).\n", node2, id, filename);
        exit_mpi();
      }

    }
  }


  // ----------------------------
  // check and store elements
  // ----------------------------
  int nElems = elemList.size();
  Es.resize(nElems);
  if(nElems != maxElem) { // need to renumber elements.
    print_warning("Warning: The element indices of an embedded surface may have a gap: "
                  "max index = %d, number of elements = %d. Renumbering elements. (%s)\n",
                  maxElem, nElems, filename);
//    assert(nElems < maxElem);
    
    int current_id = 0; 
    vector<bool> elemcheck(maxElem+1, false);
    for(auto it = elemList.begin(); it != elemList.end(); it++) {
      id = (*it)[0];
      if(elemcheck[id]) {
        print_error("*** Error: Found duplicate element (id: %d) in embedded surface file %s.\n", id, filename);
        exit_mpi();
      }
      elemcheck[id] = true;

      node1 = (*it)[1];
      node2 = (*it)[2];
      
      if(old2new.empty()) //node set is original order
        Es[current_id] = Int2(node1-1, node2-1);
      else {// nodes are renumbered
        auto p1 = old2new.find(node1);
        auto p2 = old2new.find(node2);
        Es[current_id] = Int2(p1->second, p2->second);
      }      
      current_id++;
    }
  } 
  else { //element numbers in good shape

    vector<bool> elemcheck(nElems, false);
    for(auto it = elemList.begin(); it != elemList.end(); it++) {
      id = (*it)[0] - 1;
      if(elemcheck[id]) {
        print_error("*** Error: Found duplicate element (id: %d) in embedded surface file %s.\n", id, filename);
        exit_mpi();
      }
      elemcheck[id] = true;

      node1 = (*it)[1];
      node2 = (*it)[2];
      
      if(old2new.empty()) //node set is original order
        Es[id] = Int2(node1-1, node2-1);
      else {// nodes are renumbered
        auto p1 = old2new.find(node1);
        auto p2 = old2new.find(node2);
        Es[id] = Int2(p1->second, p2->second);
      }
    }
  }

}

double Tools::min_list(vector<double> &list){
  if (list.size() == 0){ 
    cout<<"Empty List"<<endl;return 0;}
  double minimum=list[0];
  for (unsigned int i=0;i<list.size();i++){
    if (list[i]<minimum) minimum=list[i];
  }
  return minimum; 
}

double Tools::point_distance(Vec3D &p1,double &p2x,double &p2y){
  // Distance formula - d = sqrt[(x2-x1)^2 + (y2-y1)^2]
  double dist = sqrt(pow(p1.v[0]-p2x,2)+pow(p1.v[1]-p2y,2));
  return dist;  
}

void Tools::closest_point(vector<Vec3D> &d,vector<Vec3D> &es,vector<double> &Gx,vector<double> &Gy){ 
  for (unsigned int n=0;n<es.size();n++){
    int idx_min=0; //1st Grid x coord
    int idy_min=0; //1st Grid y coord
    double dist_current=Tools::point_distance(es[n],Gx[0],Gy[0]); //distance with 1st Grid point
    for (unsigned int i=0;i<Gx.size();i++){ //xcoord loop
      for (unsigned int j=0;j<Gy.size();j++){ //ycoord loop
        double dist_new = Tools::point_distance(es[n],Gx[i],Gy[j]); // calculates the distance
        if (dist_new<dist_current){ 
          dist_current=dist_new; //updates dist_current with the new smaller dist_new value 
          idx_min=i; //updates id_min with the id of the grid point with the smallest distance
          idy_min=j;
        }
      }
    }
    Vec3D closest_grid_point{Gx[idx_min],Gy[idy_min],0};
    d.push_back(closest_grid_point); //adds the grid point closest to the embedded surface point
    
  }  
}

void Tools::max_min_coords(double &xmax,double &ymax,double &xmin,double &ymin,vector<Vec3D> &es){
  xmax = es[0].v[0];
  xmin = es[0].v[0];
  ymax = es[0].v[1];
  ymin = es[0].v[1];
  for (unsigned int i=0;i<es.size();i++){
    if (es[i].v[0]>xmax) xmax=es[i].v[0]; //updates xmax with the larger value if the condition is met
    if (es[i].v[1]>ymax) ymax=es[i].v[1]; //updates ymax with the larger value if the condition is met
    if (es[i].v[0]<xmin) xmin=es[i].v[0]; //updates xmin with smaller value
    if (es[i].v[1]<ymin) ymin=es[i].v[1]; //updates ymin with smaller value
  }
}


void Tools::grid_nodes_solid_domain(double &xmax,double &ymax,double &xmin,double &ymin,double &delta_x,double &delta_y,vector<double> &xcoord,vector<double> &ycoord,vector<Vec3D> &es){
  Vec3D node; // grid node that will be filled in
  int num_xcoords = (xmax-xmin)/delta_x; //# of xcoords in solid surface
  int num_ycoords = (ymax-ymin)/delta_y; //# of ycoords in solid surface
  for (int i=0;i>num_xcoords;i++){
    
  }  
}

 
/*void Tools::flood_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color){ //algorithm obtained from https://en.wikipedia.org/wiki/Flood_fill
  //  if (double*** color is empty or DNE) error(); //checks if triple pointer color exists
  //if (color[0][j][i]!=0 || i>imax || j>jmax || i<imin || j<jmin) return;  //condition if node is already colored or out of grid domain
  if (color[0][j][i]!=0 || i>imax || j>jmax || i<imin || j<jmin) return;  //condition if node is already colored or out of grid domain
  color[0][j][i] = 1; //setting current node to a color
  flood_fill(i,j-1,imax,jmax,imin,jmin,color); //moving south by 1 node
  flood_fill(i,j+1,imax,jmax,imin,jmin,color); //moving north by 1 node
  flood_fill(i-1,j,imax,jmax,imin,jmin,color); //moving west by 1 node
  flood_fill(i+1,j,imax,jmax,imin,jmin,color); //moving east by 1 node

  return;
}*/
void Tools::flood_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color){ //algorithm obtained from https://en.wikipedia.org/wiki/Flood_fill
  //  if (double*** color is empty or DNE) error(); //checks if triple pointer color exists
  //if (blocked==true || i>imax || j>jmax || i<imin || j<jmin) return;  //condition if node is already colored or out of grid domain
  if (color[0][j][i]!=0|| i>imax || j>jmax || i<imin || j<jmin) return;  //condition if node is already colored or out of grid domain

//*************************Flood Fill (Fluid Domain)**************************************//
  color[0][j][i] = 1; //only filling in node w/ 1 if it has no intersections with embedded surface
  flood_fill(i,j-1,imax,jmax,imin,jmin,color); //moving south by 1 node
  flood_fill(i,j+1,imax,jmax,imin,jmin,color); //moving north by 1 node
  flood_fill(i-1,j,imax,jmax,imin,jmin,color); //moving east by 1 node
  flood_fill(i+1,j,imax,jmax,imin,jmin,color); //moving west by 1 node

  return;
}

void Tools::intersect_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color,vector<Vec3D> &surface_nodes,vector<Int2> &surface_connectivities,vector<double> &xcoords,vector<double> &ycoords){

  for (int j=jmin;j<jmax;j++){
    for (int i=imin;i<imax;i++){
      
  //*************************Intersection check(Fluid-Structure Interface)******************************//
      if (i==imax||i==jmax||i==imin||j==jmin) continue;
      bool intersection_south = intersect(i,j,i,j-1,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment south of node - intersection check
      if (intersection_south == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=2;color[0][j-1][i]=2;
      }
      bool intersection_north = intersect(i,j,i,j+1,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment north of node - intersection check
      if (intersection_north == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=2;color[0][j+1][i]=2;

      }
      bool intersection_west = intersect(i,j,i-1,j,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment west of node - intersection check
      if (intersection_west == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=2;color[0][j][i-1]=2;

      }
      bool intersection_east = intersect(i,j,i+1,j,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment east of node - intersection check
      if (intersection_east == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=2;color[0][j][i+1]=2;

      }
    }
  }

}

bool Tools::intersect(int grid_node1i,int grid_node1j,int grid_node2i,int grid_node2j,int &imax,int &jmax,int &imin,int &jmin,vector<Vec3D> &surface_nodes,vector<Int2> &surface_connectivities,vector<double> &xcoords,vector<double> &ycoords){ //will return eithe true or false if the nodes of the grid intersect with the embedded surface. Use the line to line intersection formula from wikipedia. Once point of intersection is found, use a condition to check if intersection point is on the same line as the grid line and if is in the bounds of the grid points. Reference: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection 
 
  if (grid_node2i>imax || grid_node2j>jmax || grid_node2i<imin || grid_node2j<jmin) return false; //condition if line segment extends further past the grid boundaries
  bool blocked = false; //setting blocked to false by default

  Vec3D grid_node1(xcoords[grid_node1i],ycoords[grid_node1j],0); //coordinates of grid nodes - Vec3D
  Vec3D grid_node2(xcoords[grid_node2i],ycoords[grid_node2j],0);

  if (grid_node1i==16 && grid_node1j==13){
    cout<<"Here is the problem!"<<endl;
    cout<<"For node:16,13"<<" X-coord. = "<<xcoords[grid_node1i]<<"\t"<<"Y_coord. = "<<ycoords[grid_node1j]<<endl;
    cout<<"For node:17,13"<<" X-coord. = "<<xcoords[grid_node2i]<<"\t"<<"Y_coord. = "<<ycoords[grid_node2j]<<endl;
    cout<<"blocked ="<<blocked<<endl;
  }
  //cout<<"gridnode1i: "<<grid_node1i<<"\t"<<"grid_node1j: "<<grid_node1j<<endl;

  for (int i=0;i<surface_connectivities.size();i++){ //comparing specific grid line segment with all line segmnents of embedded surface
    Int2 surface_lineseg = surface_connectivities[i]; //line segment of embedded surface
    int* p1 = &surface_lineseg.v[0]; int* p2 = &surface_lineseg.v[1];
    //cout<<"node1 id from connectivities = "<<*p1<<endl; 
    Vec3D surface_node1 = surface_nodes[*p1];Vec3D surface_node2 = surface_nodes[*p2]; //coordinates of surface nodes - Vec3D
    
    
    double intersect_point = point_of_intersection(grid_node1,grid_node2,surface_node1,surface_node2);
    if (intersect_point == 1){ //intersection detected
      blocked = true;
      return blocked;
    } 
  }
  return blocked;
}

double Tools::point_of_intersection(Vec3D &grid_node1,Vec3D &grid_node2,Vec3D &surface_node1,Vec3D &surface_node2){ //calculates the point of intersection between 2 line segments using the bezier parameters to compute point of intersection
  double x1=grid_node1.v[0];double x2=grid_node2.v[0];double x3=surface_node1.v[0];double x4=surface_node2.v[0]; //x-coordinates
  double y1=grid_node1.v[1];double y2=grid_node2.v[1];double y3=surface_node1.v[1];double y4=surface_node2.v[1]; //y-coordinates

  double t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //1st bezier parameter
  double u = ((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //2nd bezier parameter
 


  if (t>1 || t<0 || u>1 || u<0) return 0; //bezier parameters condition 
  
 /* cout<<"Grid Node1: x-coord:"<<x1<<"\ty-coord: "<<y1<<endl;
  cout<<"Grid Node2: x-coord:"<<x2<<"\ty-coord: "<<y2<<endl;
  cout<<"Surface Node1: x-coord:"<<x3<<"\ty-coord: "<<y3<<endl;
  cout<<"Surface Node2: x-coord:"<<x4<<"\ty-coord: "<<y4<<endl;*/
  double Px = x1+(t*(x2-x1));double Py = y1+(t*(y2-y1));
  cout<<"intersection point = "<<Px<<","<<Py<<endl;
  Vec2D point(Px,Py);
  return 1;
}

void Tools::SpaceVariable3D_print(double*** color,int &imax,int &jmax){
  cout<<"Color of all grid points:"<<endl;
  for (int i=0;i<imax+1;i++){
    for (int j=0;j<jmax+1;j++){
      cout<<"Point: "<<i<<","<<j<<" = "<<color[0][j][i]<<endl;
    }
  }
}  
