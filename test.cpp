/*
Next steps:
I've added the functions to find the midpoint between two points as Lat/Long.
Stuff to do:
-@ Ensure rounding errors don't happen at vertex mapping conversion.
-@ ability to define Tessalation level.
-@ Output final tessalation to Wavefront (.OBJ) 3D model file.
-@ Add height generation algorithm
*/

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <fstream>

using namespace std;

double pi = 3.14159265358979311599796346854;
int Tessalation_Level = 1;
const double radius = 20.0;
const double height = radius;
bool triOrQuad = 0; // if 0 the output will be quads, if 1 the output will be triangles

// Defines the latitude, longitude, and height of each vertex.
struct struct_VertexArray {
//  int v_Index;
  double v_Lat, v_Long, v_Height;
};

struct struct_FaceArray { // used to define the variables used in the various face arrays; FaceArray_current, FaceArray_initial, FaceArray_new
    int v1, v2, v3, v4;
};

struct struct_EdgeArray { // used to track what edges have been subdivided. v_Start & v_End are the start & end points by index. v_Mid is the midpoint.
// the order of the v_Start and v_End from the vertex index is W-E or N-S.
  int v_Start, v_End, v_Mid;
};

struct struct_ScratchArray { // don't need
  // the scratch array is used to temporarily store the vertex indexes for each face during calculation.
  int v1s, v2s, v3s, v4s, v5s, v6s, v7s, v8s, v9s;
};

struct struct_VectorReturn { // probably don't need
   int v1f, v2f, v3f, v4f; 
};

/*
v1s is the North point of the parent face.
v2s is the East point of the parent face.
v3s is the South point of the parent face.
v4s is the West point of the parent face.
v5s is the midpoint between v1s & v2s
v6s is the midpoint between v2s & v3s
v7s is the midpoint between v3s & v4s
v8s is the midpoint between v4s & v1s
v9s is the midpoint between v2s & v4s
*/
// ******************************* Functions ********************

// Generates the 12 vertices of a regular icosahedron, aligned so that one vertex is at the north pole of the sphere and two vertices are aligned opposite each other on the z-axis
vector<struct_VertexArray> generate_initial_icosahedron_vertices() { // function named "generate_icosahedron_vertices" using vector<Vertex> instead of void.
  vector<struct_VertexArray> VertexArray; // initializes vartype:vector using struct_VertexArray.

// vertices above and below the equator
double ico_cAngle = acos(sqrt(5.0)/5.0); // central angle of an icosahedron in radians, approx 63.4 degrees.
double x1 = (pi / 2.0) - ico_cAngle; // typical latitude of non-polar vertices of an icosahedron. Will be negative for Southern Vertices. approx 26.56 degrees.

/* 
North_Long1 = 0;
North_Long2 = (2*pi)/5; // central angle between points of a regular pentagon: 72 deg or (2pi)/5 radians
North_Long3 = (4*pi)/5;
North_Long4 = (6*pi)/5;
North_Long5 = (8*pi)/5;
South_Long1 = pi/5;
South_Long2 = (3*pi)/5;
South_Long3 = (5*pi)/5;
South_Long4 = (7*pi)/5;
South_Long5 = (9*pi)/5;
90 degrees = pi/2
*/

// Northern vertices
    VertexArray.push_back({ pi/2.0, 0.0, height });        // North pole      (0)
    VertexArray.push_back({ x1, 0.0, height });            // North point 1   (1)
    VertexArray.push_back({ x1, (2.0*pi)/5.0, height });   // North point 2   (2)
    VertexArray.push_back({ x1, (4.0*pi)/5.0, height });   // North point 3   (3)
    VertexArray.push_back({ x1, (6.0*pi)/5.0, height });   // North point 4   (4)
    VertexArray.push_back({ x1, (8.0*pi)/5.0, height });   // North point 5   (5)
// Southern vertices
    VertexArray.push_back({ -x1, pi/5.0, height });        // South point 1.5 (6)
    VertexArray.push_back({ -x1, (3.0*pi)/5.0, height });  // South point 2.5 (7)
    VertexArray.push_back({ -x1, (5.0*pi)/5.0, height });  // South point 3.5 (8)
    VertexArray.push_back({ -x1, (7.0*pi)/5.0, height });  // South point 4.5 (9)
    VertexArray.push_back({ -x1, (9.0*pi)/5.0, height });  // South point 5.5 (10)
    VertexArray.push_back({ -pi/2.0, 0.0, height });       // South pole      (11)

 return VertexArray;   
}

vector<struct_FaceArray> associate_initial_faces() { // This is a function.
    vector<struct_FaceArray> FaceArray_current;
    // Faces are associated by band in the order of West-to-East then North-to-South direction
    // before moving to the next band to the East.
    // variables: v1, v2, v3, v4; Rhombus:North/East/West/South
    FaceArray_current.push_back({ 0, 2,  6,  1 }); // index numbers were marked down by 1
    FaceArray_current.push_back({ 1, 6,  11, 10}); // since the vector starts at 0
    FaceArray_current.push_back({ 0, 3,  7,  2 });
    FaceArray_current.push_back({ 2, 7,  11, 6 });
    FaceArray_current.push_back({ 0, 4,  8,  3 });
    FaceArray_current.push_back({ 3, 8,  11, 7 });
    FaceArray_current.push_back({ 0, 5,  9,  4 });
    FaceArray_current.push_back({ 4, 9,  11, 8 });
    FaceArray_current.push_back({ 0, 1,  10, 5 });
    FaceArray_current.push_back({ 5, 10, 11, 9 });

 return FaceArray_current;
 }

// Function to find the midpoint between two points : output is vertex<double> (Lat, Long)
// input and output in radians
vector<double> midpointCalc(double Lat_1, double Long_1, double Lat_2, double Long_2) { // Function
  double Bx = cos(Lat_2)*cos(Long_2-Long_1);
  double By = cos(Lat_2)*sin(Long_2-Long_1);
  // vector<double> midpointoutput = { atan2(sin(Lat_1)+sin(Lat_2), pow(sqrt(cos(Lat_1)+Bx),2)+pow(By,2)), atan2(By,cos(Lat_1+Bx)) };
  //return midpointoutput;   
  double z = sin(Lat_1) + sin(Lat_2);
  double x = cos(Lat_1) + Bx;
  //double y = By;

  double Lat_mid = atan2(z, sqrt(x * x + By * By));
  double Long_mid = Long_1 + atan2(By, x);

  /* Check for rounding errors with the limitations of "double" values and irrational numbers. Trig functions 
  that should result in an answer of "0" might not be true due to numbers being truncated during calculation 
  because of the limited nature of digital computing with irrational numbers. The minimum value expected is 
  within the 1e-8 range given a chord length of a 10th of a meter with Dezengarth's diameter. This checks if the 
  value is less than 1e-10, which is 100 times smaller than the smallest expected value that isn't zero. 
  If the absolute value of the output is less than 1e-10, then this corrects the value to be zero.
  The minimum value is calculated to be the <CHORD LENGTH / RADIUS OF PLANET>.
  Dezengarth's diameter is 11020.155106 km (11020155.106 m) and the chord lenth is the size of the voxel (0.1 m).
  */
  if (abs(Lat_mid) < 1e-10)  // checks for rounding error if Lat_mid "should" be zero
  {
    Lat_mid = 0.0;
  }
  if (abs(Long_mid) < 1e-10) // checks for rounding error if Long_mid "should" be zero
  {
    Long_mid = 0.0;
  }

  return {Lat_mid, Long_mid};
}
 
// function to check if an edge has been subdivided and return the vertex index of the midpoint if it has
// OR return a value of -1 if it hasn't.
 int checkEdgeDivide(vector<struct_EdgeArray> &edgeList, int tempOne, int tempTwo) // 3 input arguments (name of array to check (EdgeArray), target start, target end)
 {
    // std::find_if searches for the first element in edgeList that satisfies the lambda condition
    int targetOne, targetTwo;
    if (tempOne < tempTwo) {
      targetOne = tempOne;
      targetTwo = tempTwo;
    } else {
      targetOne = tempTwo;
      targetTwo = tempOne;
    }
   
    auto it = find_if(edgeList.begin(), edgeList.end(),
        [targetOne, targetTwo](const struct_EdgeArray& edge) {
            // Return true if the current edge has the exact same start and end as the target
            return edge.v_Start == targetOne && edge.v_End == targetTwo;
        });

    // Check if an edge was found (i.e., iterator is not at the end)
    if (it != edgeList.end()) {
        return it->v_Mid;  // Return the midpoint value of the matching edge
    } else {
        return -1;         // Return -1 to indicate edge was not found
    }
}

// Function 



// **************************************************************************************

int main() {
  
// Generate the inital 12 vertices and original 10 faces.  
  vector<struct_VertexArray> VertexArray = generate_initial_icosahedron_vertices(); // Generate initial vertices
  vector<struct_FaceArray> FaceArray_current = associate_initial_faces(); // Associate initial faces
  
  
  /*/ **************************    initialization testing
  int vcount = 1;
  int fcount = 1;
  // Print the index, Lat, Long, and Height of each vertex
  for (struct_VertexArray vertex_loop : VertexArray) {
    cout << "v" << vcount << " " << vertex_loop.v_Index << " " << vertex_loop.v_Lat << " " << vertex_loop.v_Long << " " << vertex_loop.v_Height << endl;
    vcount++;  };
        
  for (struct_FaceArray face_loop : FaceArray_current) { // int f_Index, v1, v2, v3, v4;
      cout << "f" << fcount << " " << face_loop.f_Index << " " << face_loop.v1 << " " << face_loop.v2 << " " << face_loop.v3 << " " << face_loop.v4 << endl;
      fcount++;  };
  
 // midpoint test.
  cout << "Lat_1:0.0, Long_1:1.0, Lat_2:1.0, Long_2:2.0" << " : " << midpoint_Lat(0.0, 1.0, 1.0, 2.0) << endl;
  cout << "Lat_1:0.0, Long_1:1.0, Lat_2:1.0, Long_2:2.0" << " : " << midpoint_Long(0.0, 1.0, 1.0, 2.0) << endl;
  // End of initialization testing */ 

  // ******************** Start of Tessalation *********************
  
  // Step 1: Determine number of faces to subdivide and apply that to a count number
  int fcountMax = FaceArray_current.size();
  vector <struct_EdgeArray> EdgeArray;
  vector <struct_FaceArray> FaceArray_new;
    
  // Step 2: Face Subdivide Loop
for ( int fcount = 0; fcount < fcountMax; fcount++)
{
      int v_I1 = FaceArray_current.at(fcount).v1; // index of vertex 1 on parent face (North)
      int v_I2 = FaceArray_current.at(fcount).v2; // index of vertex 2 on parent face (East)
      int v_I3 = FaceArray_current.at(fcount).v3; // index of vertex 3 on parent face (South)
      int v_I4 = FaceArray_current.at(fcount).v4; // index of vertex 4 on parent face (West) 
      int v_I5 = -1;
      int edgecheck = checkEdgeDivide(EdgeArray, v_I1, v_I2);
      if (edgecheck == -1)
      {
          vector<double> midpoint_temp = midpointCalc(VertexArray.at(v_I1).v_Lat, VertexArray.at(v_I1).v_Long, VertexArray.at(v_I2).v_Lat, VertexArray.at(v_I2).v_Long);
          VertexArray.push_back({ midpoint_temp.at(0), midpoint_temp.at(1), height }); // Add vertices to VertexArray
          midpoint_temp.clear();
          v_I5 = VertexArray.size()-1;
          if (v_I1 < v_I2)
          {
          EdgeArray.push_back({ v_I1, v_I2, v_I5 }); // add calculated edge to edgeArray
          }
          else 
          {
            EdgeArray.push_back({ v_I2, v_I1, v_I5 });
          }
      }
      else
      {
        v_I5 = checkEdgeDivide(EdgeArray, v_I1, v_I2);
      }
    
      int v_I6 = -1;
      if (checkEdgeDivide(EdgeArray, v_I2,v_I3)==-1)
      {
          vector<double> midpoint_temp = midpointCalc(VertexArray.at(v_I2).v_Lat, VertexArray.at(v_I2).v_Long, VertexArray.at(v_I3).v_Lat, VertexArray.at(v_I3).v_Long);
          VertexArray.push_back({ midpoint_temp.at(0), midpoint_temp.at(1), height });
          midpoint_temp.clear();
          v_I6 = VertexArray.size()-1;
          EdgeArray.push_back({ v_I2, v_I3, v_I6 }); // add calculated edge to edgeArray
      }
      else
      {
        v_I6 = checkEdgeDivide(EdgeArray, v_I2,v_I3);
      }
      int v_I7 = -1;
      if (checkEdgeDivide(EdgeArray, v_I3,v_I4)==-1)
      {
          vector<double> midpoint_temp = midpointCalc(VertexArray.at(v_I3).v_Lat, VertexArray.at(v_I3).v_Long, VertexArray.at(v_I4).v_Lat, VertexArray.at(v_I4).v_Long);
          VertexArray.push_back({ midpoint_temp.at(0), midpoint_temp.at(1), height });
          midpoint_temp.clear();
          v_I7 = VertexArray.size()-1;
          EdgeArray.push_back({ v_I3, v_I4, v_I7 }); // add calculated edge to edgeArray
      }
      else
      {
        v_I7 = checkEdgeDivide(EdgeArray, v_I3,v_I4);
      }
      int v_I8 = -1;
      if (checkEdgeDivide(EdgeArray, v_I4,v_I1)==-1)
      {
          vector<double> midpoint_temp = midpointCalc(VertexArray.at(v_I4).v_Lat, VertexArray.at(v_I4).v_Long, VertexArray.at(v_I1).v_Lat, VertexArray.at(v_I1).v_Long);
          VertexArray.push_back({ midpoint_temp.at(0), midpoint_temp.at(1), height });
          midpoint_temp.clear();
          v_I8 = VertexArray.size()-1;
          EdgeArray.push_back({ v_I4, v_I1, v_I8 }); // add calculated edge to edgeArray
      }
      else
      {
        v_I8 = checkEdgeDivide(EdgeArray, v_I4,v_I1);
      }
// int v_I9 = -1;
          vector<double> midpoint_temp = midpointCalc(VertexArray.at(v_I2).v_Lat, VertexArray.at(v_I2).v_Long, VertexArray.at(v_I4).v_Lat, VertexArray.at(v_I4).v_Long);
          VertexArray.push_back({ midpoint_temp.at(0), midpoint_temp.at(1), height });
          midpoint_temp.clear();
          int v_I9 = VertexArray.size()-1;
          EdgeArray.push_back({ v_I2, v_I4, v_I9});
// add new faces to FaceArray_new
FaceArray_new.push_back ({v_I1, v_I5, v_I9, v_I8}); // Face 1
FaceArray_new.push_back ({v_I5, v_I2, v_I6, v_I9}); // Face 2
FaceArray_new.push_back ({v_I9, v_I6, v_I3, v_I7}); // Face 3
FaceArray_new.push_back ({v_I8, v_I9, v_I7, v_I4}); // Face 4

}

FaceArray_current.clear();
FaceArray_current = FaceArray_new;
FaceArray_new.clear();
EdgeArray.clear();
  
/*
// test vertex output
  int vcount = 0;
  for (struct_VertexArray vertex_loop : VertexArray) {
    cout << "v" << vcount << " " << vertex_loop.v_Lat << " " << vertex_loop.v_Long << " " << vertex_loop.v_Height;
    cout << " lat: " << vertex_loop.v_Lat * (180/pi) << " long: " << vertex_loop.v_Long * (180/pi);
    cout << endl;
    vcount++;  };
// test face output
  int fcount = 0;
  cout << endl << "FaceArray_current:" << endl;
  for (struct_FaceArray face_loop : FaceArray_current) {
    cout << "f" << fcount << " " << face_loop.v1 << " " << face_loop.v2 << " " << face_loop.v3 << " " << face_loop.v4 << endl;
    fcount++; };
// convert vertices to 3D Cartesian coordinates
// rounding error needs to be corrected for numbers approaching zero
  vcount = 0;
  cout << endl;
  for (struct_VertexArray vertex_loop: VertexArray) {
    cout << "v" << vcount << " X: " << vertex_loop.v_Height * cos(vertex_loop.v_Lat) * cos(vertex_loop.v_Long) << " Y: " << vertex_loop.v_Height * cos(vertex_loop.v_Lat) * sin(vertex_loop.v_Long) << " Z: " << vertex_loop.v_Height * sin(vertex_loop.v_Lat) << endl;
    vcount++;
  }
 
 cout << endl << "End of Program" << endl << endl;
 */  

// outputting the strings to a text file

/*
Wavefront .OBJ File Format - Quick Reference
============================================

# Comment line (ignored by parsers)
# This is a comment

v x y z [w]
    Vertex position. 'w' is optional and defaults to 1.0.

vt u [v] [w]
    Texture coordinate. 'v' is optional (defaults to 0.0). 'w' rarely used.

vn x y z
    Vertex normal vector.

vp u [v] [w]
    Parameter space vertex (used for curves and surfaces).

f v1 v2 v3 ...
    Face using vertex indices only.

f v1/vt1 v2/vt2 v3/vt3 ...
    Face using vertex + texture indices.

f v1//vn1 v2//vn2 v3//vn3 ...
    Face using vertex + normal indices.

f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...
    Face using vertex + texture + normal indices.

o object_name
    Object name — groups subsequent elements into a logical object.

g group_name
    Group name — useful for organizing parts of a model.

usemtl material_name
    Use the specified material for subsequent faces.


mtllib filename.mtl
    Reference an external material library.
s 1
s off
    Enable (s 1, s 2, etc.) or disable (s off) smoothing groups.

Additional Notes:
-----------------
- Indices are 1-based (first vertex is '1', not '0').
- Negative indices refer to relative positions from the current line.
- Faces can be triangles, quads, or n-gons (3+ vertices).
- Lines can be continued with a backslash `\` at the end.
- vertices are indexed by the order they appear, starting at 1
*/

ofstream outFile("test.obj"); // Create and open a file named output.txt

    if (outFile.is_open()) {
        outFile << "# icosahedron test\n";
        outFile << "# This is your first file output.\n";
        int vcount = 0;
        cout << endl;

        // add vertices to output file.
        for (struct_VertexArray vertex_loop: VertexArray) {
          outFile << "v " << vertex_loop.v_Height * cos(vertex_loop.v_Lat) * cos(vertex_loop.v_Long) << " " << vertex_loop.v_Height * cos(vertex_loop.v_Lat) * sin(vertex_loop.v_Long) << " " << vertex_loop.v_Height * sin(vertex_loop.v_Lat) << endl;
          vcount++;
        }

        // add faces to output file
        int fcount = 0;
        if (triOrQuad = 1 ) {
        fcount = 0;
        outFile << endl << "# Faces - Triangles" << endl;
        for (struct_FaceArray face_loop : FaceArray_current) {
          // 1,4,2
          // 4,3,2
        }
        }
        else {
        fcount = 0;
        outFile << endl << "# Faces - Quads" << endl;
        for (struct_FaceArray face_loop : FaceArray_current) {
          outFile << "f " << face_loop.v1+1 << " " << face_loop.v4+1 << " " << face_loop.v3+1 << " " << face_loop.v2+1 << endl;
          fcount++; };
        }
      

        outFile.close(); // Always close the file when done
        cout << "File written successfully.\n";
    } else {
        cerr << "Unable to open file for writing.\n";
    }


  return 0;
}
