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


//------------TOOLS CLASS FUNCTION DEFINITIONS--------------------- 
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
    exit_mpi(); }

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

double Tools::point_distance(Vec3D &pt1,Vec3D &pt2){
  // Distance formula - d = sqrt[(x2-x1)^2 + (y2-y1)^2]
  double dist = sqrt(pow(pt1.v[0]-pt2.v[0],2)+pow(pt1.v[1]-pt2.v[1],2));
  return dist;  
}

Vec3D Tools::closest_point(Vec3D &pt1,vector<Vec3D> &pts){ 
  double dist=100; Vec3D closest_pt;double dist_compare; 
  for (unsigned int n=0;n<pts.size();n++){
    dist_compare = point_distance(pt1,pts[n]);
    if (dist_compare == 0) continue; //distance of zero means the same point
    if (dist_compare < dist){ //condition if a closer point is found
      dist = dist_compare;
      closest_pt = pts[n];
    }
  }
  return closest_pt;
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

int Tools::getIndex(vector<Vec3D> &v, Vec3D K) 
{ 

  for (int i=0;i<v.size();i++){
    double dist = point_distance(K,v[i]);
    if (dist == 0) return i; //distance of 0 means same point
  }
  return -1; //case if K is not in v
}

vector<Vec3D> Tools::remove_Vec3D(vector<Vec3D> &v, Vec3D &N){
  //if (v.size == 0) return -- perform error handling for this case
  vector<Vec3D> trunc; //truncated vector that is returned
  for (int i=0;i<v.size();i++){
    double dist = point_distance(N,v[i]);
    if (dist==0) continue; //distance of 0 means same point
    trunc.push_back(v[i]);
  }
  return trunc;
}

vector<Int2> Tools::remove_duplicates_Int2(vector<Int2> &rep){
  if (rep.size() < 2) return rep; //case if vector is too small to have repeated elements

  vector<Int2> no_rep;
  int* x1,*y1,*x2,*y2;
  int tag;
  for (int i=0;i<rep.size();i++){ //looping through original repeating vector
    tag = 0;
    if (no_rep.size()==0){ //1st element case
      no_rep.push_back(rep[i]);
      continue;
    }
    x1 = &rep[i].v[0]; y1 = &rep[i].v[1];
    for (int j=0;j<no_rep.size();j++){ //looping through non-repeating list
      x2 = &no_rep[j].v[0]; y2 = &no_rep[j].v[1]; 
      if ((*x1==*x2 || *x1==*y2) && (*y1==*x2 || *y1==*y2)){
        tag = 1;break;
      }
    }
    if (tag==1) continue; //repeated element found 
    no_rep.push_back(rep[i]);
  }
  
  return no_rep;
}

void Tools::grid_nodes_solid_domain(double &xmax,double &ymax,double &xmin,double &ymin,double &delta_x,double &delta_y,vector<double> &xcoord,vector<double> &ycoord,vector<Vec3D> &es){
  Vec3D node; // grid node that will be filled in
  int num_xcoords = (xmax-xmin)/delta_x; //# of xcoords in solid surface
  int num_ycoords = (ymax-ymin)/delta_y; //# of ycoords in solid surface
  for (int i=0;i>num_xcoords;i++){
    
  }  
}

 
void Tools::flood_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color){ //algorithm obtained from https://en.wikipedia.org/wiki/Flood_fill
  //  if (double*** color is empty or DNE) error(); //checks if triple pointer color exists
  //if (blocked==true || i>imax || j>jmax || i<imin || j<jmin) return;  //condition if node is already colored or out of grid domain
  if (i>imax || j>jmax || i<imin || j<jmin) return;  //condition if node is already colored or out of grid domain
  if (color[0][j][i]!=0) return;  //condition if node is already colored or out of grid domain

//*************************Flood Fill (Fluid Domain)**************************************//
  color[0][j][i] = 1; //only filling in node w/ 1 if it has no intersections with embedded surface
  flood_fill(i,j-1,imax,jmax,imin,jmin,color); //moving south by 1 node
  flood_fill(i,j+1,imax,jmax,imin,jmin,color); //moving north by 1 node
  flood_fill(i-1,j,imax,jmax,imin,jmin,color); //moving west by 1 node
  flood_fill(i+1,j,imax,jmax,imin,jmin,color); //moving east by 1 node

  return;
}

void Tools::intersect_fill(int i,int j,int &imax,int &jmax,int &imin,int &jmin,double*** color,int &color_val,vector<Vec3D> &surface_nodes,vector<Int2> &surface_connectivities,vector<double> &xcoords,vector<double> &ycoords,vector<Vec3D> &intersecting_nodes,vector<Int2> &intersecting_edges){

  for (int j=jmin;j<jmax;j++){
    for (int i=imin;i<imax;i++){
      
  //*************************Intersection check(Fluid-Structure Interface)******************************//
  //***BLOCKED meaning - line connnecting one node from another is "blocked" an intersecting line from the embedded surface******
      if (i==imax||i==jmax||i==imin||j==jmin) continue;
      bool intersection_south = intersect(i,j,i,j-1,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment south of node - intersection check
      Vec3D P1, P2; //declaration of intersecting nodes
      Int2 edge; //declaration of intersecting edge

      if (intersection_south == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=color_val;color[0][j-1][i]=color_val;

        P1[0] = xcoords[i]; P1[1] = ycoords[j];
        P2[0] = xcoords[i]; P2[1] = ycoords[j-1];
  
        intersecting_nodes.push_back(P1);       
        intersecting_nodes.push_back(P2);       
        
        edge[0] = intersecting_nodes.size()-1;
        edge[1] = intersecting_nodes.size();
        intersecting_edges.push_back(edge);
      }
      bool intersection_north = intersect(i,j,i,j+1,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment north of node - intersection check
      if (intersection_north == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=color_val;color[0][j+1][i]=color_val;

        P1[0] = xcoords[i]; P1[1] = ycoords[j];
        P2[0] = xcoords[i]; P2[1] = ycoords[j+1];
  
        intersecting_nodes.push_back(P1);       
        intersecting_nodes.push_back(P2);       
        
        edge[0] = intersecting_nodes.size()-1;
        edge[1] = intersecting_nodes.size();
        intersecting_edges.push_back(edge);

      }
      bool intersection_west = intersect(i,j,i-1,j,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment west of node - intersection check
      if (intersection_west == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=color_val;color[0][j][i-1]=color_val;

        P1[0] = xcoords[i]; P1[1] = ycoords[j];
        P2[0] = xcoords[i-1]; P2[1] = ycoords[j];
  
        intersecting_nodes.push_back(P1);       
        intersecting_nodes.push_back(P2);       
        
        edge[0] = intersecting_nodes.size()-1;
        edge[1] = intersecting_nodes.size();
        intersecting_edges.push_back(edge);

      }
      bool intersection_east = intersect(i,j,i+1,j,imax,jmax,imin,jmin,surface_nodes,surface_connectivities,xcoords,ycoords); //grid line segment east of node - intersection check
      if (intersection_east == true){
        //cout<<"Line segment colored\n";
        color[0][j][i]=color_val;color[0][j][i+1]=color_val;

        P1[0] = xcoords[i]; P1[1] = ycoords[j];
        P2[0] = xcoords[i+1]; P2[1] = ycoords[j];
  
        intersecting_nodes.push_back(P1);       
        intersecting_nodes.push_back(P2);       
        
        edge[0] = intersecting_nodes.size()-1;
        edge[1] = intersecting_nodes.size();
        intersecting_edges.push_back(edge);

      }
    }
  }

}

bool Tools::intersect(int grid_node1i,int grid_node1j,int grid_node2i,int grid_node2j,int &imax,int &jmax,int &imin,int &jmin,vector<Vec3D> &surface_nodes,vector<Int2> &surface_connectivities,vector<double> &xcoords,vector<double> &ycoords){ //will return eithe true or false if the nodes of the grid intersect with the embedded surface. Use the line to line intersection formula from wikipedia. Once point of intersection is found, use a condition to check if intersection point is on the same line as the grid line and if is in the bounds of the grid points. Reference: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection 
  if (grid_node2i>imax || grid_node2j>jmax || grid_node2i<imin || grid_node2j<jmin) return false; //condition if line segment extends further past the grid boundaries
  bool blocked = false; //setting blocked to false by default

  Vec3D grid_node1(xcoords[grid_node1i],ycoords[grid_node1j],0); //coordinates of grid nodes - Vec3D
  Vec3D grid_node2(xcoords[grid_node2i],ycoords[grid_node2j],0);

  //cout<<"gridnode1i: "<<grid_node1i<<"\t"<<"grid_node1j: "<<grid_node1j<<endl;

  for (long unsigned int i=0;i<surface_connectivities.size();i++){ //comparing specific grid line segment with all line segmnents of embedded surface
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
// website: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
  double x1=grid_node1.v[0];double x2=grid_node2.v[0];double x3=surface_node1.v[0];double x4=surface_node2.v[0]; //x-coordinates
  double y1=grid_node1.v[1];double y2=grid_node2.v[1];double y3=surface_node1.v[1];double y4=surface_node2.v[1]; //y-coordinates


  //case if grid line and surface lines are parallel to each other
  if (x2<x3 && x1<x3 && x2<x4 && x1<x4) return 0; // parallel in x
  if (x2>x3 && x1>x3 && x2>x4 && x1>x4) return 0;

  if (y2<y3 && y1<y3 && y2<y4 && y1<y4) return 0; // parallel in y
  if (y2>y3 && y1>y3 && y2>y4 && y1>y4) return 0;

  double t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //t value (non-dimensional)
  double u = ((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //u value (non-dimensional)
 


  if (t>1 || t<0 || u>1 || u<0) return 0; //line segement condition 
  
 /* cout<<"Grid Node1: x-coord:"<<x1<<"\ty-coord: "<<y1<<endl;
  cout<<"Grid Node2: x-coord:"<<x2<<"\ty-coord: "<<y2<<endl;
  cout<<"Surface Node1: x-coord:"<<x3<<"\ty-coord: "<<y3<<endl;
  cout<<"Surface Node2: x-coord:"<<x4<<"\ty-coord: "<<y4<<endl;*/
  double Px = x1+(t*(x2-x1));double Py = y1+(t*(y2-y1));
  //cout<<"intersection point = "<<Px<<","<<Py<<endl;
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

//---------------TOPOLOGY TRACKER FUNCTION DEFINITIONS-------------

//MAJOR CHANGE UP FOR THIS!!!: Use algorithm provided by geeks for geeks or at least incorporate some
//of it. Hopefully, it should take care of the collinear points with the edges
bool Topology::is_inside(Vec3D &point,vector<Int2> &elements,vector<Vec3D> &nodes){
  Tools tool;

// ray casting of point
  Vec3D point_cast{100,point.v[1],point.v[2]}; //corresponding point to A in ray cast
  //now comparing if point is inside shape
  int count = 0;
  vector<Int2> elem_collision_history;
  vector<Vec3D> node_collision_history;
  double xmin,xmax,ymin,ymax;
  for (int i=0;i<elements.size();i++){
    Vec3D node1 = nodes[elements[i].v[0]];
    Vec3D node2 = nodes[elements[i].v[1]];
    vector<Vec3D> nodes{node1,node2};
    tool.max_min_coords(xmax,ymax,xmin,ymin,nodes);
    
    if (point.v[1] > ymin){ //checking if point is above y min

      if (point.v[1] <= ymax && point.v[0] <= xmax){ //checking if point is below or eq to ymax & left or eq to xmax
        Vec3D pt_of_intersect = pt_of_intersection(point,point_cast,node1,node2);
        if (point.v[0] == node1.v[0] || point.v[0] <= pt_of_intersect.v[0]){ //checks if point is on same line of element or left or eq to the intersection point
          count++; 
        }
      }
    }
  }

/* repeated elements case
  int same_line = 0;
  if (elem_collision_history.size() >= 2){
    elem_collision_history = tool.remove_duplicates_Int2(elem_collision_history);   
    
    //for lines that are drawn on the same line b/c of STL file :)
    double* x_const = &nodes[elem_collision_history[0].v[0]].v[0];
    double* y_const = &nodes[elem_collision_history[0].v[0]].v[1];
    for (int i=1;i<elem_collision_history.size();i++){
      double* x1 = &nodes[elem_collision_history[i].v[0]].v[0];
      double* y1 = &nodes[elem_collision_history[i].v[0]].v[1];
      if (*x1==*x_const || *y1==*y_const) same_line++; 
    }

  }
  //cout<<"elem_collision_history size: "<<elem_collision_history.size()<<endl;
  for (int i=0;i<elem_collision_history.size();i++){
    //cout<<"Collision with elem_collision_history after non-repeating function:\t"<<"["<<elem_collision_history[i].v[0]<<","<<elem_collision_history[i].v[1]<<"]"<<endl;
  } 
  if (same_line != 0) count = elem_collision_history.size() -1;
  else count = elem_collision_history.size(); //updates the count to the "actual" count*/
  //cout<<"Collision number: "<<count<<endl;
  //cout<<"final count: "<<count<<endl;
  if (count % 2 == 0 || count == 0) return false; // returns false if collision number is even or 0
  return true;
}

vector<Vec3D> Topology::inside_points() {
  vector<Vec3D> in_pts; //creating the inside points vector that will be returned
//checking all shape A points in shape B
  for (long unsigned int i=0;i<shape1_nodes.size();i++){
    Vec3D pt_a = shape1_nodes[i]; //shape 1 point

    if (is_inside(pt_a,shape2_elements,shape2_nodes) == true) in_pts.push_back(pt_a); // adds point to inside points list condition
  }

//checking all shape B points in shape A
  for (long unsigned int i=0;i<shape2_nodes.size();i++){
    Vec3D pt_b = shape2_nodes[i]; //shape 1 point

    if (is_inside(pt_b,shape1_elements,shape1_nodes) == true) in_pts.push_back(pt_b); // adds point to inside points list condition
  }

  return in_pts; //returns the list of inside points
}

Vec3D Topology::pt_of_intersection(Vec3D &pt1_A,Vec3D &pt2_A,Vec3D &pt1_B,Vec3D &pt2_B){

  double x1=pt1_A.v[0];double x2=pt2_A.v[0];double x3=pt1_B.v[0];double x4=pt2_B.v[0]; //x-coordinates
  double y1=pt1_A.v[1];double y2=pt2_A.v[1];double y3=pt1_B.v[1];double y4=pt2_B.v[1]; //y-coordinates

  Vec3D no_intersect(0,0,-1); // no intersection value

  //special case if line segments are parallel to each other
  if ((x1==x3 && x2==x4) || (y1==y3 && y2==y4)) return no_intersect; //this type of point will NOT be classified as an intersecting point

  double t = ((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //t value (non-dimensional)
  double u = ((x1-x3)*(y1-y2)-(y1-y3)*(x1-x2))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)); //u value (non-dimensional)
 


  if (t>1 || t<0 || u>1 || u<0){  ; //condition if there is no intersection
    return no_intersect;
  }
  
  double Px = x1+(t*(x2-x1));double Py = y1+(t*(y2-y1));
  Vec3D intersect_point(Px,Py,0);
  return intersect_point;

}

vector<Vec3D> Topology::intersecting_points(){
  vector<Vec3D> intersecting_pts; //vector of intersecting points that gets returned
  for (long unsigned int i=0;i<shape1_elements.size();i++){
    //Vec3D points of shape 1 element
    Vec3D pt1_a = shape1_nodes[shape1_elements[i].v[0]];
    Vec3D pt2_a = shape1_nodes[shape1_elements[i].v[1]];
    for (long unsigned int j=0;j<shape2_elements.size();j++){
      //Vec3D points of shape 2 element
      Vec3D pt1_b = shape2_nodes[shape2_elements[j].v[0]];
      Vec3D pt2_b = shape2_nodes[shape2_elements[j].v[1]];
      //calculates the point of intersection, if any
      Vec3D intersect_pt = pt_of_intersection(pt1_a,pt2_a,pt1_b,pt2_b);
      if (intersect_pt.v[2] == -1){ //no intersection condition
        continue;
      }
      intersecting_pts.push_back(intersect_pt); //appends calculated point of intersection to intersecting points vector
    }  
  }
  return intersecting_pts;
}

void Topology::shape_construct(vector<Int2> &overlap_connectivities,vector<Vec3D> &overlap_nodes,vector<Vec3D> &intersecting_pts,vector<Vec3D> &inside_pts){
  if (intersecting_pts.size() == 0 && inside_pts.size() == 0) { //case if the shapes do not intersect
    cout<<"Shapes do not intersect!"<<endl;
    return;
  }
  //creation of unsorted list of interesecting + inside points - all_unsorted_pts
  for (int i=0;i<intersecting_pts.size();i++){ 
    overlap_nodes.push_back(intersecting_pts[i]);
  }
  for (int i=0;i<inside_pts.size();i++){ 
    overlap_nodes.push_back(inside_pts[i]);
  }
  vector<Vec3D> all_unsorted_pts = overlap_nodes; //vector that contains all unsorted points of overlap shape
  Vec3D pt1,pt_closest;
  Tools tool;
  Int2 element;
  //cout<<"\nbeginning overlap shape construction: \n";
  //construction of overlap nodes vector
  pt1 = all_unsorted_pts[0]; //starts off at 1st point of unsorted list
  while (all_unsorted_pts.size() != 1) {
  
    //cout<<"unsorted list size: "<<all_unsorted_pts.size()<<endl;

    pt_closest = tool.closest_point(pt1,all_unsorted_pts);
    //pt1.print("Point 1");
    //pt_closest.print("Closest point");
 
    //get index of closest point and pt 1
    int index_pt1 = tool.getIndex(overlap_nodes,pt1); //?
    int index_pt_closest = tool.getIndex(overlap_nodes,pt_closest);

    //make element out of them
    element.v[0] = index_pt1; element.v[1] = index_pt_closest;

    //append to the overlap_elements list
    overlap_connectivities.push_back(element);

    all_unsorted_pts = tool.remove_Vec3D(all_unsorted_pts,pt1); //removes the point that was looked at
    pt1 = pt_closest;


  }
  //adds the last element from last point to 1st of overlap_nodes
  element.v[0] = overlap_connectivities.back().v[1]; element.v[1] = overlap_connectivities[0].v[0]; 
  overlap_connectivities.push_back(element);
  //last point equals final point in all_unsorted_list
  //make element of this point and 1st point of 1st element in overlap_elements
  //overlap_connectivities.push_back(element);
  
  return;
}

void Topology::shape_fill(vector<Int2> &overlap_connectivities,vector<Vec3D> &overlap_nodes,vector<double> &xcoords,vector<double> &ycoords,double*** &color, int color_val,int imax,int jmax,int i0,int j0){
  Tools tool; //tool object creation
  //need to define bounding box on overlap shape so grid nodes in bounding box wil be compared
  //cout<<"NOW STARTING SHAPE FILL"<<endl;
  double xmax,xmin,ymax,ymin; //max and min coordinates of overlap shape
  tool.max_min_coords(xmax,ymax,xmin,ymin,overlap_nodes); 
  /*cout<<"OVERLAP NODES: "<<endl;
  for (int i=0;i<overlap_nodes.size();i++){
    cout<<"Node: "<<overlap_nodes[i].v[0]<<" "<<overlap_nodes[i].v[1]<<" "<<overlap_nodes[i].v[2]<<endl;
  }*/

  vector <double> box_xcoords,box_ycoords; //these are the x and y grid coordinates that are in the bounding box
  for (long unsigned int i=0;i<xcoords.size();i++){
    if (xcoords[i] >= xmin && xcoords[i] <= xmax) box_xcoords.push_back(xcoords[i]); //case if grid x coord is in overlap bounding box
    if (ycoords[i] >= ymin && ycoords[i] <= ymax) box_ycoords.push_back(ycoords[i]); //case if grid y coord s in overlp bounding box
  }
 
  //cout<<"box coords number: "<<box_xcoords.size()*box_ycoords.size()<<endl;
  //cout<<"These are the points that are in the bounding box"<<endl;  
  //cout<<"bounding box size: "<<box_xcoords.size()<<"x\t"<<box_ycoords.size()<<endl;
  /*for (int i=0;i<box_xcoords.size();i++){
    cout<<"box point: "<<box_xcoords[i]<<"\t"<<box_ycoords[i]<<endl;
  } */
  //now coloring bounding box points that are inside the overlap shape
  for (int i=0;i<box_xcoords.size();i++){ 
    for (int j=0;j<box_ycoords.size();j++){
      Vec3D grid_point{box_xcoords[i],box_ycoords[j],0};
   //   cout<<"Grid Point: "<<grid_point.v[0]<<"\t"<<grid_point.v[1]<<endl;
  //    cout<<"grid point before the inside check:"<<endl;
    //  cout<<"box point: "<<box_xcoords[i]<<"\t"<<box_ycoords[j]<<endl;
   
/*   for (int i=i0;i<imax;i++){
     for (int j=j0;j<jmax;j++){
       Vec3D grid_point{xcoords[i],ycoords[j],0};
       bool inside = is_inside(grid_point,overlap_connectivities,overlap_nodes); 
       if (inside==true)color[0][j][i]=color_val;
     }
   }
*/
      bool inside = is_inside(grid_point,overlap_connectivities,overlap_nodes); 
      if (inside == true){
        //finds the index of the box coord in the grid coords list to determine node # for color
        auto itx = find(xcoords.begin(),xcoords.end(),box_xcoords[i]);
        auto ity = find(ycoords.begin(),ycoords.end(),box_ycoords[j]);
        int xindex = itx - xcoords.begin();
        int yindex = ity - ycoords.begin();
        //cout<<"xindex: "<<xindex<<endl;
        //cout<<"yindex: "<<yindex<<endl;
        color[0][yindex][xindex] = color_val;
       }
     } 
   } 
 
  return;
}

