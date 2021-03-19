// This program starts from an icosahedron and refines its mesh.
// Each initial face is subdivided in n^2 faces (n integer >2).
// The mesh of the sphere with given radius is written to a mesh file.
// The corresponding red blood cell mesh with the same major radius is written to another mesh file.

#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "p_vec_ten.hh"

#define SQ2(x) ((x) * (x))

using namespace std;

/*************************************************************
*********** MAIN FUNCTION ************************************
*************************************************************/

int main() {

  /// SETUP
  /// =====

for(int num_division = 1; num_division <= 30; ++num_division) {

  /// Set radius and degree of refinement

  double radius = 1; // radius of sphere and red blood cell
//  int num_division = 6; // number of triangles along one edge ('1' corresponds to no refinement)

  cout << "divide each edge into " << num_division << " segments" << endl;

  /// Compute initial and final number of nodes, faces, and edges

  int num_faces_i = 20; // initial number of faces for icosahedron
  int num_faces_f = 20 * SQ2(num_division); // final number of faces
//  int num_nodes_i = (num_faces_i + 4) / 2; // initial number of nodes (12 for icosahedron)
  int num_nodes_f = (num_faces_f + 4) / 2; // final number of nodes
  int num_edges_i = (3 * num_faces_i) / 2; // initial number of edges (30 for icosahedron)
//  int num_edges_f = (3 * num_faces_f) / 2; // final number of edges

  cout << "resulting in " << num_faces_f << " faces" << endl;

  /// Allocate memory for node identification look-up table
  // this is a central look-up table
  // it is used to find the correct faces at the end
  // during the refinement, it is successively filled with the node indexes within the face
  // the first array index is used for the index of the initial face (0-19)
  // the second and third array indexes are used for triangular integer coordinates
  // within the face ('n_1' and 'n_2' later on)
  // the array is initialized to negative values (= not used) in order to find problems faster

  int node_identification[num_faces_i][num_division + 1][num_division + 1];

  for(int f_i = 0; f_i < num_faces_i; ++f_i) {
    for(int n_1 = 0; n_1 < num_division + 1; ++ n_1) {
      for(int n_2 = 0; n_2 < num_division + 1; ++ n_2) {
        node_identification[f_i][n_1][n_2] = -1;
      }
    }
  }

  /// Allocate memory and initialize look-up table 'face => nodes in face'
  // each face has three nodes
  // value '-1' means 'not yet used'
  // define the first 20 known faces of the icosahedron (indexes 0-19)
  // this data will be overwritten later when the final mesh is constructed

  int nodes_in_face[num_faces_f][3];

  for(int f_i = 0; f_i < num_faces_f; ++f_i) {
    nodes_in_face[num_faces_f][0] = -1;
    nodes_in_face[num_faces_f][1] = -1;
    nodes_in_face[num_faces_f][2] = -1;
  }

  nodes_in_face[ 0][0] =  0, nodes_in_face[ 0][1] =  8, nodes_in_face[ 0][2] =  4;
  nodes_in_face[ 1][0] =  0, nodes_in_face[ 1][1] =  5, nodes_in_face[ 1][2] = 10;
  nodes_in_face[ 2][0] =  2, nodes_in_face[ 2][1] =  4, nodes_in_face[ 2][2] =  9;
  nodes_in_face[ 3][0] =  2, nodes_in_face[ 3][1] = 11, nodes_in_face[ 3][2] =  5;
  nodes_in_face[ 4][0] =  1, nodes_in_face[ 4][1] =  6, nodes_in_face[ 4][2] =  8;
  nodes_in_face[ 5][0] =  1, nodes_in_face[ 5][1] = 10, nodes_in_face[ 5][2] =  7;
  nodes_in_face[ 6][0] =  3, nodes_in_face[ 6][1] =  9, nodes_in_face[ 6][2] =  6;
  nodes_in_face[ 7][0] =  3, nodes_in_face[ 7][1] =  7, nodes_in_face[ 7][2] = 11;
  nodes_in_face[ 8][0] =  0, nodes_in_face[ 8][1] = 10, nodes_in_face[ 8][2] =  8;
  nodes_in_face[ 9][0] =  1, nodes_in_face[ 9][1] =  8, nodes_in_face[ 9][2] = 10;
  nodes_in_face[10][0] =  2, nodes_in_face[10][1] =  9, nodes_in_face[10][2] = 11;
  nodes_in_face[11][0] =  3, nodes_in_face[11][1] =  9, nodes_in_face[11][2] = 11;
  nodes_in_face[12][0] =  4, nodes_in_face[12][1] =  2, nodes_in_face[12][2] =  0;
  nodes_in_face[13][0] =  5, nodes_in_face[13][1] =  0, nodes_in_face[13][2] =  2;
  nodes_in_face[14][0] =  6, nodes_in_face[14][1] =  1, nodes_in_face[14][2] =  3;
  nodes_in_face[15][0] =  7, nodes_in_face[15][1] =  3, nodes_in_face[15][2] =  1;
  nodes_in_face[16][0] =  8, nodes_in_face[16][1] =  6, nodes_in_face[16][2] =  4;
  nodes_in_face[17][0] =  9, nodes_in_face[17][1] =  4, nodes_in_face[17][2] =  6;
  nodes_in_face[18][0] = 10, nodes_in_face[18][1] =  5, nodes_in_face[18][2] =  7;
  nodes_in_face[19][0] = 11, nodes_in_face[19][1] =  7, nodes_in_face[19][2] =  5;

  /// Put initial nodes in array 'node_identification'
  // the initial nodes can already be put into the identification array
  // by convention,
  // the first node in the face (0) has coordinates (0, 0)
  // the second node in the face (1) has coordinates (0, num_division)
  // the third node in the face (2) has coordinates (num_division, 0)

  for(int f_i = 0; f_i < num_faces_i; ++f_i) {
    node_identification[f_i][0][0] = nodes_in_face[f_i][0];
    node_identification[f_i][0][num_division] = nodes_in_face[f_i][1];
    node_identification[f_i][num_division][0] = nodes_in_face[f_i][2];
  }

  /// Allocate memory for initial edges and initialize look-up table 'edge => nodes in edge' and 'face => edges in face'
  // in the following, it is easier to switch to an edge-based approach
  // for this purpose, one has to know which nodes are in which edges
  // and which faces belong to which edge and vice versa
  // use three look-up tables for that and initialize them for the icosahedron
  // there are two nodes at the ends of an edge (nodes_in_edges),
  // three edges in a face (edges_in_faces),
  // and two faces at an edge (face_at_edge)

  int nodes_in_edges[num_edges_i][2];
  int edges_in_faces[num_faces_i][3];
  int face_at_edge[num_edges_i][2];

  nodes_in_edges[ 0][0] =  0, nodes_in_edges[ 0][1] =  8;
  nodes_in_edges[ 1][0] =  0, nodes_in_edges[ 1][1] =  4;
  nodes_in_edges[ 2][0] =  4, nodes_in_edges[ 2][1] =  8;
  nodes_in_edges[ 3][0] =  0, nodes_in_edges[ 3][1] =  5;
  nodes_in_edges[ 4][0] =  0, nodes_in_edges[ 4][1] = 10;
  nodes_in_edges[ 5][0] =  5, nodes_in_edges[ 5][1] = 10;
  nodes_in_edges[ 6][0] =  2, nodes_in_edges[ 6][1] =  4;
  nodes_in_edges[ 7][0] =  2, nodes_in_edges[ 7][1] =  9;
  nodes_in_edges[ 8][0] =  4, nodes_in_edges[ 8][1] =  9;
  nodes_in_edges[ 9][0] =  2, nodes_in_edges[ 9][1] = 11;
  nodes_in_edges[10][0] =  2, nodes_in_edges[10][1] =  5;
  nodes_in_edges[11][0] =  5, nodes_in_edges[11][1] = 11;
  nodes_in_edges[12][0] =  1, nodes_in_edges[12][1] =  6;
  nodes_in_edges[13][0] =  1, nodes_in_edges[13][1] =  8;
  nodes_in_edges[14][0] =  6, nodes_in_edges[14][1] =  8;
  nodes_in_edges[15][0] =  1, nodes_in_edges[15][1] = 10;
  nodes_in_edges[16][0] =  1, nodes_in_edges[16][1] =  7;
  nodes_in_edges[17][0] =  7, nodes_in_edges[17][1] = 10;
  nodes_in_edges[18][0] =  3, nodes_in_edges[18][1] =  9;
  nodes_in_edges[19][0] =  3, nodes_in_edges[19][1] =  6;
  nodes_in_edges[20][0] =  6, nodes_in_edges[20][1] =  9;
  nodes_in_edges[21][0] =  3, nodes_in_edges[21][1] =  7;
  nodes_in_edges[22][0] =  3, nodes_in_edges[22][1] = 11;
  nodes_in_edges[23][0] =  7, nodes_in_edges[23][1] = 11;
  nodes_in_edges[24][0] =  8, nodes_in_edges[24][1] = 10;
  nodes_in_edges[25][0] =  9, nodes_in_edges[25][1] = 11;
  nodes_in_edges[26][0] =  0, nodes_in_edges[26][1] =  2;
  nodes_in_edges[27][0] =  1, nodes_in_edges[27][1] =  3;
  nodes_in_edges[28][0] =  4, nodes_in_edges[28][1] =  6;
  nodes_in_edges[29][0] =  5, nodes_in_edges[29][1] =  7;

  edges_in_faces[ 0][0] =  2, edges_in_faces[ 0][1] =  1, edges_in_faces[ 0][2] =  0;
  edges_in_faces[ 1][0] =  5, edges_in_faces[ 1][1] =  4, edges_in_faces[ 1][2] =  3;
  edges_in_faces[ 2][0] =  8, edges_in_faces[ 2][1] =  7, edges_in_faces[ 2][2] =  6;
  edges_in_faces[ 3][0] = 11, edges_in_faces[ 3][1] = 10, edges_in_faces[ 3][2] =  9;
  edges_in_faces[ 4][0] = 14, edges_in_faces[ 4][1] = 13, edges_in_faces[ 4][2] = 12;
  edges_in_faces[ 5][0] = 17, edges_in_faces[ 5][1] = 16, edges_in_faces[ 5][2] = 15;
  edges_in_faces[ 6][0] = 20, edges_in_faces[ 6][1] = 19, edges_in_faces[ 6][2] = 18;
  edges_in_faces[ 7][0] = 23, edges_in_faces[ 7][1] = 22, edges_in_faces[ 7][2] = 21;
  edges_in_faces[ 8][0] = 24, edges_in_faces[ 8][1] =  0, edges_in_faces[ 8][2] =  4;
  edges_in_faces[ 9][0] = 24, edges_in_faces[ 9][1] = 15, edges_in_faces[ 9][2] = 13;
  edges_in_faces[10][0] = 25, edges_in_faces[10][1] =  9, edges_in_faces[10][2] =  7;
  edges_in_faces[11][0] = 25, edges_in_faces[11][1] = 22, edges_in_faces[11][2] = 18;
  edges_in_faces[12][0] = 26, edges_in_faces[12][1] =  1, edges_in_faces[12][2] =  6;
  edges_in_faces[13][0] = 26, edges_in_faces[13][1] = 10, edges_in_faces[13][2] =  3;
  edges_in_faces[14][0] = 27, edges_in_faces[14][1] = 19, edges_in_faces[14][2] = 12;
  edges_in_faces[15][0] = 27, edges_in_faces[15][1] = 16, edges_in_faces[15][2] = 21;
  edges_in_faces[16][0] = 28, edges_in_faces[16][1] =  2, edges_in_faces[16][2] = 14;
  edges_in_faces[17][0] = 28, edges_in_faces[17][1] = 20, edges_in_faces[17][2] =  8;
  edges_in_faces[18][0] = 29, edges_in_faces[18][1] = 17, edges_in_faces[18][2] =  5;
  edges_in_faces[19][0] = 29, edges_in_faces[19][1] = 11, edges_in_faces[19][2] = 23;

  face_at_edge[ 0][0] =  0, face_at_edge[ 0][1] =  8;
  face_at_edge[ 1][0] =  0, face_at_edge[ 1][1] = 12;
  face_at_edge[ 2][0] =  0, face_at_edge[ 2][1] = 16;
  face_at_edge[ 3][0] =  1, face_at_edge[ 3][1] = 13;
  face_at_edge[ 4][0] =  1, face_at_edge[ 4][1] =  8;
  face_at_edge[ 5][0] =  1, face_at_edge[ 5][1] = 18;
  face_at_edge[ 6][0] =  2, face_at_edge[ 6][1] = 12;
  face_at_edge[ 7][0] =  2, face_at_edge[ 7][1] = 10;
  face_at_edge[ 8][0] =  2, face_at_edge[ 8][1] = 17;
  face_at_edge[ 9][0] =  3, face_at_edge[ 9][1] = 10;
  face_at_edge[10][0] =  3, face_at_edge[10][1] = 13;
  face_at_edge[11][0] =  3, face_at_edge[11][1] = 19;
  face_at_edge[12][0] =  4, face_at_edge[12][1] = 14;
  face_at_edge[13][0] =  4, face_at_edge[13][1] =  9;
  face_at_edge[14][0] =  4, face_at_edge[14][1] = 16;
  face_at_edge[15][0] =  5, face_at_edge[15][1] =  9;
  face_at_edge[16][0] =  5, face_at_edge[16][1] = 15;
  face_at_edge[17][0] =  5, face_at_edge[17][1] = 18;
  face_at_edge[18][0] =  6, face_at_edge[18][1] = 11;
  face_at_edge[19][0] =  6, face_at_edge[19][1] = 14;
  face_at_edge[20][0] =  6, face_at_edge[20][1] = 17;
  face_at_edge[21][0] =  7, face_at_edge[21][1] = 15;
  face_at_edge[22][0] =  7, face_at_edge[22][1] = 11;
  face_at_edge[23][0] =  7, face_at_edge[23][1] = 19;
  face_at_edge[24][0] =  8, face_at_edge[24][1] =  9;
  face_at_edge[25][0] = 10, face_at_edge[25][1] = 11;
  face_at_edge[26][0] = 12, face_at_edge[26][1] = 13;
  face_at_edge[27][0] = 14, face_at_edge[27][1] = 15;
  face_at_edge[28][0] = 16, face_at_edge[28][1] = 17;
  face_at_edge[29][0] = 18, face_at_edge[29][1] = 19;

  /// Allocate memory for node positions
  // those positions define the spatial shape of the mesh at the end
  // the mesh can later be defined for a sphere, a red blood cell or other shapes

  p_vec<> node_pos[num_nodes_f];

  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
    node_pos[n_i].reset();
  }

  /// Set the positions of the initial nodes
  // define the first 12 known nodes of the icosahedron (indexes 0-11)
  // note: the radius of the icosahedron is not unity at the beginning
  // the radius will be normalized when all nodes are known

  double phi = (1. + sqrt(5.)) / 2.;

  node_pos[ 0].set_to(+phi, +1., 0.);
  node_pos[ 1].set_to(-phi, +1., 0.);
  node_pos[ 2].set_to(+phi, -1., 0.);
  node_pos[ 3].set_to(-phi, -1., 0.);
  node_pos[ 4].set_to(+1., 0., +phi);
  node_pos[ 5].set_to(+1., 0., -phi);
  node_pos[ 6].set_to(-1., 0., +phi);
  node_pos[ 7].set_to(-1., 0., -phi);
  node_pos[ 8].set_to(0., +phi, +1.);
  node_pos[ 9].set_to(0., -phi, +1.);
  node_pos[10].set_to(0., +phi, -1.);
  node_pos[11].set_to(0., -phi, -1.);

  for(int n_i = 0; n_i < 12; ++n_i) {
    node_pos[n_i].norm_to(radius);
  }

  /// MESH GENERATION
  /// ===============

  /// Create new nodes on edges
  // on each of the initial edges, new nodes are created
  // there are '(num_division - 1)' new nodes on each initial edge
  // take subsequent node indexes after 11 for those, i.e., starting at 12
  // run through all 30 edges
  // on each edge, run through all 'num_division - 1' new nodes

  int n_j = 12;

  for(int e_i = 0; e_i < num_edges_i; ++e_i) {
    // identify the two end nodes belonging to the edge

    int node_start = nodes_in_edges[e_i][0];
    int node_end   = nodes_in_edges[e_i][1];

    // identify the two faces belonging to the edge

    int face_0 = face_at_edge[e_i][0];
    int face_1 = face_at_edge[e_i][1];

    for(int n = 0; n < num_division - 1; ++n) {
      // the node position is computed from the positions of the nodes
      // at the end of the edge

      node_pos[n_j] = (node_pos[node_start] * (num_division - 1 - n) + node_pos[node_end] * (n + 1)) / num_division;

      // the more complicated part is to correctly update the identification array
      // here, the correct coordinates (n_1, n_2) of the new edge nodes is computed
      // the coordinates (n_1, n_2) depend on the location of the edge in the face
      // and the location of the end nodes in the edge

      // do this for the first face belonging to the edge

      if(e_i == edges_in_faces[face_0][0]) {
        if(node_start == nodes_in_face[face_0][1]) {
          node_identification[face_0][n + 1][num_division - 1 - n] = n_j;
        }
        else if(node_start == nodes_in_face[face_0][2]) {
          node_identification[face_0][num_division - 1 - n][n + 1] = n_j;
        }
      }
      else if(e_i == edges_in_faces[face_0][1]) {
        if(node_start == nodes_in_face[face_0][0]) {
          node_identification[face_0][n + 1][0] = n_j;
        }
        else if(node_start == nodes_in_face[face_0][2]) {
          node_identification[face_0][num_division - 1 - n][0] = n_j;
        }
      }
      else if(e_i == edges_in_faces[face_0][2]) {
        if(node_start == nodes_in_face[face_0][0]) {
          node_identification[face_0][0][n + 1] = n_j;
        }
        else if(node_start == nodes_in_face[face_0][1]) {
          node_identification[face_0][0][num_division - 1 - n] = n_j;
        }
      }

      // do this for the second face belonging to the edge

      if(e_i == edges_in_faces[face_1][0]) {
        if(node_start == nodes_in_face[face_1][1]) {
          node_identification[face_1][n + 1][num_division - 1 - n] = n_j;
        }
        else if(node_start == nodes_in_face[face_1][2]) {
          node_identification[face_1][num_division - 1 - n][n + 1] = n_j;
        }
      }
      else if(e_i == edges_in_faces[face_1][1]) {
        if(node_start == nodes_in_face[face_1][0]) {
          node_identification[face_1][n + 1][0] = n_j;
        }
        else if(node_start == nodes_in_face[face_1][2]) {
          node_identification[face_1][num_division - 1 - n][0] = n_j;
        }
      }
      else if(e_i == edges_in_faces[face_1][2]) {
        if(node_start == nodes_in_face[face_1][0]) {
          node_identification[face_1][0][n + 1] = n_j;
        }
        else if(node_start == nodes_in_face[face_1][1]) {
          node_identification[face_1][0][num_division - 1 - n] = n_j;
        }
      }

      /// Shift the nodes to the sphere
      // there are different ways to do this
      // 1. simply shift along the radial vector
      // 2. shift in such a way that the arc length is unchanged

      /// Method 1

//      node_pos[n_j].norm_to(radius);

      /// Method 2

      double psi = acos(node_pos[node_start] * node_pos[node_end]); // total arc length
      double theta = psi * (n + 1) / num_division; // arc length between new node and corner node
      double beta = sin(theta) / sin(psi); // coefficients for linear combination...
      double alpha = cos(theta) - beta * cos(psi); // ...of new node position

      node_pos[n_j] = node_pos[node_start] * alpha + node_pos[node_end] * beta; // apply linear combination

      n_j++;
    }
  }

  /// Create new nodes in faces
  // there are '(num_division + 1) * (num_division + 2) / 2 - 3 * num_division' new nodes on each initial face
  // take subsequent node indexes after the last used from the former computations
  // run through all 20 initial faces
  // on each face, run through all '(num_division + 1) * (num_division + 2) / 2 - 3 * num_division' new nodes
  // note: this number is only positive for 'num_division > 2'

  for(int f_i = 0; f_i < num_faces_i; ++f_i) {
    // identify nodes in face

    int node_0 = nodes_in_face[f_i][0];
    int node_1 = nodes_in_face[f_i][1];
    int node_2 = nodes_in_face[f_i][2];

    // define unit vectors for the coordinates (n_1, n_2) in order to correctly compute
    // the coordinates of the new nodes

    p_vec<> unit_vector_1 = (node_pos[node_1] - node_pos[node_0]) / num_division;
    p_vec<> unit_vector_2 = (node_pos[node_2] - node_pos[node_0]) / num_division;

    // run through the 2D face with coordinates (n_1, n_2)

    for(int n_1 = 0; n_1 < num_division - 2; ++n_1) {
      for(int n_2 = 0; n_2 < num_division - 2 - n_1; ++n_2) {
        // distribute the new nodes on the initial face in such a way
        // that they are located in a regular triangular mesh
        // update the node identification table at the same time

        node_pos[n_j] = node_pos[node_0] + unit_vector_1 * (n_2 + 1) + unit_vector_2 * (n_1 + 1);
        node_identification[f_i][n_1 + 1][n_2 + 1] = n_j;

        /// Shift the nodes to the sphere
        // there are two methods to project the nodes on the sphere
        // 1. simply shift along the radial vector
        // 2. shift in such a way that the arc length is as constant as possible

        /// Method 1

//        node_pos[n_j].norm_to(radius);

        /// Method 2

        int n1 = n_1 + 1;
        int n2 = n_2 + 1;

        // define three lines in the initial face which contain the point (n1, n2)
        // each line is parallel to one of the edges of the initial face,
        // and it is defined by the start and end points

        int node_l0_e1 = node_identification[f_i][n1 + n2][0];            // parallel to edge...
        int node_l0_e2 = node_identification[f_i][0][n1 + n2];            // ...opposite of node 0
        int node_l1_e1 = node_identification[f_i][0][n2];                 // parallel to edge...
        int node_l1_e2 = node_identification[f_i][num_division - n2][n2]; // ...opposite of node 1
        int node_l2_e1 = node_identification[f_i][n1][0];                 // parallel to edge...
        int node_l2_e2 = node_identification[f_i][n1][num_division - n1]; // ...opposite of node 2

        // define the normal vectors of the three planes containing the above lines and the origin

        p_vec<> normal_l0 = node_pos[node_l0_e1] % node_pos[node_l0_e2];
        p_vec<> normal_l1 = node_pos[node_l1_e1] % node_pos[node_l1_e2];
        p_vec<> normal_l2 = node_pos[node_l2_e1] % node_pos[node_l2_e2];

        // define three direction vectors which are normal to two of the three above planes, respectively

        p_vec<> direction_0 = normal_l2 % normal_l1;
        p_vec<> direction_1 = normal_l0 % normal_l2;
        p_vec<> direction_2 = normal_l0 % normal_l1;

        direction_0.norm_to(radius);
        direction_1.norm_to(radius);
        direction_2.norm_to(radius);

        // make sure that all direction vectors point outwards

        if(direction_0 * node_pos[n_j] < 0) {
          direction_0 *= -1;
        }
        if(direction_1 * node_pos[n_j] < 0) {
          direction_1 *= -1;
        }
        if(direction_2 * node_pos[n_j] < 0) {
          direction_2 *= -1;
        }

        // the projected position of the new node is the average of the intersections
        // of each of the direction vectors and the surface of the sphere

        node_pos[n_j] = direction_0 + direction_1 + direction_2;
        node_pos[n_j].norm_to(radius);

        n_j++;
      }
    }
  }

  // at this point, all new nodes are correctly created including
  // their physical position on the icosahedron (node_pos),
  // their geometrical position on each initial face (node_identification)

  /// IDENTIFY FACES
  /// ==============

  // as the last major step, the faces have to be constructed from the known nodes
  // the look-up table 'node_identification' is used for this purpose
  // each face consist of three neighboring nodes in anti-cyclic order
  // such that all normal vectors point outwards

  int f_j = 0;

  for(int f_i = 0; f_i < num_faces_i; ++f_i) {

    /// 'up-standing' faces
    // these are the faces with the baseline at the bottom (triangle points 'upwards')

    for(int n_1 = 0; n_1 < num_division; ++n_1) {
      for(int n_2 = 0; n_2 < num_division - n_1; ++n_2) {
        nodes_in_face[f_j][0] = node_identification[f_i][n_1][n_2];
        nodes_in_face[f_j][1] = node_identification[f_i][n_1][n_2 + 1];
        nodes_in_face[f_j][2] = node_identification[f_i][n_1 + 1][n_2];
        f_j++;
      }
    }

    /// 'down-standing' faces
    // these are the faces with the baseline at the top (triangle points 'downwards')

    for(int n_1 = 1; n_1 < num_division + 1; ++n_1) {
      for(int n_2 = 0; n_2 < num_division - n_1; ++n_2) {
        nodes_in_face[f_j][0] = node_identification[f_i][n_1][n_2];
        nodes_in_face[f_j][1] = node_identification[f_i][n_1 - 1][n_2 + 1];
        nodes_in_face[f_j][2] = node_identification[f_i][n_1][n_2 + 1];
        f_j++;
      }
    }
  }

  // now, the geometrical mesh is completely known, but the nodes can still be shifted
  // in physical space in order to match a given shape

  /// CREATE SPHERICAL SHAPE AND WRITE TO MESH FILE
  /// =============================================

  /// Shift node positions to unit sphere

//  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
//    node_pos[n_i].norm_to(radius);
//  }

  /// Write geometry to mesh file

  stringstream filename;
  filename << "sph_ico_" << num_faces_f << ".msh";
  ofstream meshfile(filename.str().c_str());

  meshfile << "$MeshFormat\n";
  meshfile << "2 0 8\n";
  meshfile << "$EndMeshFormat\n";
  meshfile << "$Nodes\n";
  meshfile << num_nodes_f << "\n";
  meshfile << fixed << setprecision(15);

  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
    meshfile << n_i + 1 << " " << node_pos[n_i] << "\n";
  }

  meshfile << "$EndNodes\n";
  meshfile << "$Elements\n";
  meshfile << num_faces_f << "\n";

  for(int f_i = 0; f_i < num_faces_f; ++f_i) {
    meshfile << f_i + 1 << " 2 3 0 1 0 " << nodes_in_face[f_i][0] + 1 << " " << nodes_in_face[f_i][1] + 1 << " " << nodes_in_face[f_i][2] + 1 << "\n";
  }

  meshfile << "$EndElements\n";
  meshfile.close();

  /// CREATE ELLIPSOIDAL SHAPE AND WRITE TO MESH FILE
  /// ===============================================

  /// Define geometry of an ellipsoid

  double R = radius;
  double r = 5. * R;
  double x, y, z, X, Y, Z; // coordinates
  p_vec<> coord_sph, coord_ell_old, coord_ell_new;

  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
    // start with coordinates for a sphere

    coord_sph = node_pos[n_i];

    // deform the sphere to the initial ellipsoid

    coord_ell_old.x = coord_sph.x;
    coord_ell_old.y = coord_sph.y;

    if(abs(coord_sph.z) > 1e-4) {
      coord_ell_old.z = r / R * sqrt(SQ2(R) - SQ2(coord_sph.x) - SQ2(coord_sph.y));
    }
    else {
      coord_ell_old.z = 0;
    }

    // set correct sign for z-coordinate of initial ellipsoid
    // correct for possible numerical problems from shape function

    if(coord_sph.z < 0) {
      coord_ell_old.z *= -1;
    }

    if(isnan(coord_ell_old.z) || abs(coord_ell_old.z) < 1e-4) {
      coord_ell_old.z = 0;
    }

//     // compute inclination angles of the
//     // 1. spherical node relative to the xy-plane
//     // 2. ellipsoidal node relative to the xy-plane
// 
//     double inclination_circle;
//     double inclination_ellipsoid;
// 
//     if(sqrt(SQ2(coord_sph.x) + SQ2(coord_sph.y)) > 1e-4) {
//       inclination_circle = atan(coord_sph.z / sqrt(SQ2(coord_sph.x) + SQ2(coord_sph.y)));
//       inclination_ellipsoid = atan(coord_ell_old.z / sqrt(SQ2(coord_ell_old.x) + SQ2(coord_ell_old.y)));
//     }
//     else {
//       inclination_circle = 0.5 * M_PI;
//       inclination_ellipsoid = 0.5 * M_PI;
//     }
// 
//     if(sqrt(SQ2(coord_sph.x) + SQ2(coord_sph.y)) > 1e-4) {
//       double weight_s = 0.7;
//       double weight_m = 0.75;
//       double weight_e = 1.0;
// 
//       double c = weight_s;
//       double a = 2 * ((weight_e - c) - 2 * (weight_m - c));
//       double b = weight_e - c - a;
//       double angle = 2 * abs(inclination_circle) / M_PI;
// 
//       inclination_circle *= (a * SQ2(angle) + b * angle + c);
//     }
// 
//     // find new radius for the node position on the ellipsoid
// 
//     if(abs(coord_ell_old.z) > 1e-4 && sqrt(SQ2(coord_sph.x) + SQ2(coord_sph.y)) > 1e-4) {
//       coord_ell_new.z = sqrt(SQ2(r) / (1 + SQ2(r / (R * tan(inclination_circle)))));
// 
//       double length_new = abs(coord_ell_new.z / tan(inclination_circle));
// 
//       coord_ell_new.x = sqrt(SQ2(length_new) / (1 + SQ2(coord_ell_old.y / coord_ell_old.x)));
//       coord_ell_new.y = sqrt(SQ2(length_new) / (1 + SQ2(coord_ell_old.x / coord_ell_old.y)));
// 
//       if(coord_sph.x < 0) {
//         coord_ell_new.x *= -1;
//       }
// 
//       if(coord_sph.y < 0) {
//         coord_ell_new.y *= -1;
//       }
// 
//       if(coord_sph.z < 0) {
//         coord_ell_new.z *= -1;
//       }
//     }
//     else {
      coord_ell_new = coord_ell_old;
//     }

    node_pos[n_i] = coord_ell_new;
  }

  /// Write geometry to mesh file

  stringstream filename_ell;
  filename_ell << "ell_ico_" << num_faces_f << ".msh";
  ofstream meshfile_ell(filename_ell.str().c_str());

  meshfile_ell << "$MeshFormat\n";
  meshfile_ell << "2 0 8\n";
  meshfile_ell << "$EndMeshFormat\n";
  meshfile_ell << "$Nodes\n";
  meshfile_ell << num_nodes_f << "\n";
  meshfile_ell << fixed << setprecision(15);

  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
    meshfile_ell << n_i + 1 << " " << node_pos[n_i] << "\n";
  }

  meshfile_ell << "$EndNodes\n";
  meshfile_ell << "$Elements\n";
  meshfile_ell << num_faces_f << "\n";

  for(int f_i = 0; f_i < num_faces_f; ++f_i) {
    meshfile_ell << f_i + 1 << " 2 3 0 1 0 " << nodes_in_face[f_i][0] + 1 << " " << nodes_in_face[f_i][1] + 1 << " " << nodes_in_face[f_i][2] + 1 << "\n";
  }

  meshfile_ell << "$EndElements\n";
  meshfile_ell.close();

  /// CREATE RED BLOOD CELL SHAPE AND WRITE TO MESH FILE
  /// ==================================================

  /// Define geometry of a red blood cell

  double C_0 = 0.207; // those are shape parameters
  double C_1 = 2.0;
  double C_2 = -1.123;
//  double x, y, z, Z; // coordinates

  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
    x = node_pos[n_i].x;
    y = node_pos[n_i].y;
    z = node_pos[n_i].z;

    // correct the z-coordinate as function of the x- and y-coordinates

    Z = (0.5 * sqrt(1. - (SQ2(x) + SQ2(y))) * (C_0 + C_1 * (SQ2(x) + SQ2(y)) + C_2 * SQ2(SQ2(x) + SQ2(y))));

    // set correct sign for z-coordinate

    if(z < 0) {
      Z *= -1;
    }

    // correct for possible numerical problems from shape function

    if(isnan(Z) || abs(Z) < 1e-4) {
      Z = 0;
    }

    // overwrite the coordinate by the computed one

    node_pos[n_i].z = Z;
  }

  /// Write geometry to mesh file

  stringstream filename_rbc;
  filename_rbc << "rbc_ico_" << num_faces_f << ".msh";
  ofstream meshfile_rbc(filename_rbc.str().c_str());

  meshfile_rbc << "$MeshFormat\n";
  meshfile_rbc << "2 0 8\n";
  meshfile_rbc << "$EndMeshFormat\n";
  meshfile_rbc << "$Nodes\n";
  meshfile_rbc << num_nodes_f << "\n";
  meshfile_rbc << fixed << setprecision(15);

  for(int n_i = 0; n_i < num_nodes_f; ++n_i) {
    meshfile_rbc << n_i + 1 << " " << node_pos[n_i] << "\n";
  }

  meshfile_rbc << "$EndNodes\n";
  meshfile_rbc << "$Elements\n";
  meshfile_rbc << num_faces_f << "\n";

  for(int f_i = 0; f_i < num_faces_f; ++f_i) {
    meshfile_rbc << f_i + 1 << " 2 3 0 1 0 " << nodes_in_face[f_i][0] + 1 << " " << nodes_in_face[f_i][1] + 1 << " " << nodes_in_face[f_i][2] + 1 << "\n";
  }

  meshfile_rbc << "$EndElements\n";
  meshfile_rbc.close();

  /// Report successful completion of program

  cout << "done" << endl;

}

  return 0;
}
