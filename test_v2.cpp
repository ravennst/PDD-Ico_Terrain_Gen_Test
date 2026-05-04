/*
 * test_v2.cpp - Icosahedron Tessellation with Structured Grid Addressing
 *
 * Rewrite of test.cpp incorporating the Array.txt methodology.
 *
 * KEY ARCHITECTURAL CHANGE:
 * The original code used a flat VertexArray with insertion-order indices,
 * requiring an EdgeArray + O(n) linear scan (checkEdgeDivide) to detect
 * already-computed shared-edge midpoints.
 *
 * This version assigns every vertex a DETERMINISTIC index based on its
 * structural position in the icosahedron tessellation:
 *
 *   Indices 0..11               : 12 primary icosahedron corner vertices
 *   Indices 12..12+20*(t-1)-1   : edge-interior vertices, grouped by edge
 *   Indices 12+20*(t-1)..end    : face-interior vertices, grouped by face
 *
 * For tessellation level t, given face f, row r, col c (0-indexed):
 *   - Corners  (r==0||r==t) && (c==0||c==t) -> primary vertex via FACE_CORNERS
 *   - Edges    exactly one of r==0, r==t, c==0, c==t   -> FACE_EDGES lookup
 *   - Interior otherwise                                -> arithmetic formula
 *
 * Result: midpoint lookup is O(1) arithmetic. No EdgeArray. No search.
 * Each vertex's planet() call is guarded by a computed flag so shared
 * edge vertices are never evaluated twice.
 *
 * VERTEX LAYOUT within a face at level t (example t=3):
 *
 *     col: 0    1    2    3
 * row 0:   NW   .    .    NE    <- north edge (r=0)
 * row 1:   .    i    i    .
 * row 2:   .    i    i    .
 * row 3:   SW   .    .    SE    <- south edge (r=t)
 *          ^                ^
 *        west             east
 *        edge             edge
 *
 * North corner = face corner N, East = E, South = S, West = W
 * (matching the N/E/S/W convention from associate_initial_faces)
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cassert>

extern "C" {
  #include "libs/Planet/planet.h"
}

// ---------------------------------------------------------------------------
// planet.h globals
// ---------------------------------------------------------------------------
double rseed   = 0.21;
double M       = 0.021;
planet_vertex tetra[4];

using namespace std;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
const double PI        = 3.141592653589793;
const double RADIUS    = 1.0;
const double HEIGHT_MOD = 1.0;

int  Tessellation_Level = 7;   // target tessellation depth
bool triOrQuad          = true; // true = triangles, false = quads

// Calc_Level controls planet() recursion depth.
// Higher values = more terrain detail but slower per-call.
// The +15 offset from the original is preserved here; tune as needed.
int Calc_Level; // set in main after Tessellation_Level is known

// ---------------------------------------------------------------------------
// Structures
// ---------------------------------------------------------------------------
struct Vertex {
    double lat, lon, height;
    bool   computed = false; // has planet() been called for this slot?
};

struct Face {
    int v[4]; // indices: v[0]=N, v[1]=E, v[2]=S, v[3]=W
};

// ---------------------------------------------------------------------------
// Coordinate helpers
// ---------------------------------------------------------------------------
struct XYZ { double x, y, z; };

XYZ ll_to_xyz(double lat, double lon) {
    return { cos(lat)*sin(lon), sin(lat), cos(lat)*cos(lon) };
}

// Great-circle midpoint, with rounding-error correction near zero.
// Returns {lat, lon} in radians.
pair<double,double> midpoint(double lat1, double lon1,
                              double lat2, double lon2)
{
    double dLon = lon2 - lon1;
    double Bx   = cos(lat2) * cos(dLon);
    double By   = cos(lat2) * sin(dLon);
    double z    = sin(lat1) + sin(lat2);
    double x    = cos(lat1) + Bx;

    double lat_m = atan2(z, sqrt(x*x + By*By));
    double lon_m = lon1 + atan2(By, x);

    // Rounding-error correction: values that should be exactly zero
    // can have floating-point residue. Threshold is well below the
    // smallest expected non-zero value for Dezengarth's resolution.
    if (fabs(lat_m) < 1e-10) lat_m = 0.0;
    if (fabs(lon_m) < 1e-10) lon_m = 0.0;
    if (lon_m >  PI) lon_m -= 2.0*PI;
    if (lon_m < -PI) lon_m += 2.0*PI;

    return {lat_m, lon_m};
}

// ---------------------------------------------------------------------------
// Topology tables
// ---------------------------------------------------------------------------

// 12 primary icosahedron vertices as {lat, lon} in radians.
// Indices match the original code's VertexArray order.
//   0  = North pole
//   1..5  = Northern ring
//   6..10 = Southern ring
//   11 = South pole
const double ICO_C = acos(sqrt(5.0)/5.0); // central angle ~63.43 deg
const double X1    = PI/2.0 - ICO_C;      // ring latitude ~26.57 deg

struct LatLon { double lat, lon; };
const LatLon PRIMARY[12] = {
    {  PI/2.0,          0.0           }, // 0  North pole
    {  X1,              0.0           }, // 1
    {  X1,    (2.0*PI)/5.0            }, // 2
    {  X1,    (4.0*PI)/5.0            }, // 3
    {  X1,    (6.0*PI)/5.0            }, // 4
    {  X1,    (8.0*PI)/5.0            }, // 5
    { -X1,     PI/5.0                 }, // 6
    { -X1,    (3.0*PI)/5.0            }, // 7
    { -X1,    (5.0*PI)/5.0            }, // 8
    { -X1,    (7.0*PI)/5.0            }, // 9
    { -X1,    (9.0*PI)/5.0            }, // 10
    { -PI/2.0,          0.0           }, // 11 South pole
};

// 10 faces: each row is {N, E, S, W} primary vertex indices.
// Matches associate_initial_faces() from the original code exactly.
const int FACE_CORNERS[10][4] = {
    //  N   E   S   W
    {   0,  2,  6,  1 }, // face 0
    {   1,  6, 11, 10 }, // face 1
    {   0,  3,  7,  2 }, // face 2
    {   2,  7, 11,  6 }, // face 3
    {   0,  4,  8,  3 }, // face 4
    {   3,  8, 11,  7 }, // face 5
    {   0,  5,  9,  4 }, // face 6
    {   4,  9, 11,  8 }, // face 7
    {   0,  1, 10,  5 }, // face 8
    {   5, 10, 11,  9 }, // face 9
};

// 20 canonical edges as {lo_vertex, hi_vertex} with lo < hi.
// Edge index is the position in this array.
const int EDGES[20][2] = {
    {  0,  1 }, // e0
    {  0,  2 }, // e1
    {  0,  3 }, // e2
    {  0,  4 }, // e3
    {  0,  5 }, // e4
    {  1,  6 }, // e5
    {  1, 10 }, // e6
    {  2,  6 }, // e7
    {  2,  7 }, // e8
    {  3,  7 }, // e9
    {  3,  8 }, // e10
    {  4,  8 }, // e11
    {  4,  9 }, // e12
    {  5,  9 }, // e13
    {  5, 10 }, // e14
    {  6, 11 }, // e15
    {  7, 11 }, // e16
    {  8, 11 }, // e17
    {  9, 11 }, // e18
    { 10, 11 }, // e19
};

// Per-face edge table.
// For each face (0..9) and each of its 4 sides (NE, ES, SW, WN),
// stores {edge_index, forward}.
// 'forward' = true  means the canonical edge direction (lo→hi) runs
//                   in the same direction as the face traversal.
// 'forward' = false means it runs reversed.
//
// Face side order: 0=NE (north→east), 1=ES (east→south),
//                  2=SW (south→west), 3=WN (west→north)
//
// Traversal direction within the face grid:
//   NE side: row 0,    col 0→t   (west to east along north row)
//   ES side: col t,    row 0→t   (north to south along east col)
//   SW side: row t,    col t→0   (east to west along south row) [reversed]
//   WN side: col 0,    row t→0   (south to north along west col)[reversed]

struct EdgeRef { int edge_id; bool forward; };

const EdgeRef FACE_EDGES[10][4] = {
    // face 0: N=0, E=2, S=6, W=1
    //   NE: 0→2 = e1 fwd,  ES: 2→6 = e7 fwd,  SW: 6→1 = e5 rev,  WN: 1→0 = e0 rev
    { {1,true}, {7,true}, {5,false}, {0,false} },
    // face 1: N=1, E=6, S=11, W=10
    //   NE: 1→6=e5 fwd, ES: 6→11=e15 fwd, SW: 11→10=e19 rev, WN: 10→1=e6 rev
    { {5,true}, {15,true}, {19,false}, {6,false} },
    // face 2: N=0, E=3, S=7, W=2
    //   NE: 0→3=e2 fwd, ES: 3→7=e9 fwd, SW: 7→2=e8 rev, WN: 2→0=e1 rev
    { {2,true}, {9,true}, {8,false}, {1,false} },
    // face 3: N=2, E=7, S=11, W=6
    //   NE: 2→7=e8 fwd, ES: 7→11=e16 fwd, SW: 11→6=e15 rev, WN: 6→2=e7 rev
    { {8,true}, {16,true}, {15,false}, {7,false} },
    // face 4: N=0, E=4, S=8, W=3
    //   NE: 0→4=e3 fwd, ES: 4→8=e11 fwd, SW: 8→3=e10 rev, WN: 3→0=e2 rev
    { {3,true}, {11,true}, {10,false}, {2,false} },
    // face 5: N=3, E=8, S=11, W=7
    //   NE: 3→8=e10 fwd, ES: 8→11=e17 fwd, SW: 11→7=e16 rev, WN: 7→3=e9 rev
    { {10,true}, {17,true}, {16,false}, {9,false} },
    // face 6: N=0, E=5, S=9, W=4
    //   NE: 0→5=e4 fwd, ES: 5→9=e13 fwd, SW: 9→4=e12 rev, WN: 4→0=e3 rev
    { {4,true}, {13,true}, {12,false}, {3,false} },
    // face 7: N=4, E=9, S=11, W=8
    //   NE: 4→9=e12 fwd, ES: 9→11=e18 fwd, SW: 11→8=e17 rev, WN: 8→4=e11 rev
    { {12,true}, {18,true}, {17,false}, {11,false} },
    // face 8: N=0, E=1, S=10, W=5
    //   NE: 0→1=e0 fwd, ES: 1→10=e6 fwd, SW: 10→5=e14 rev, WN: 5→0=e4 rev
    { {0,true}, {6,true}, {14,false}, {4,false} },
    // face 9: N=5, E=10, S=11, W=9
    //   NE: 5→10=e14 fwd, ES: 10→11=e19 fwd, SW: 11→9=e18 rev, WN: 9→5=e13 rev
    { {14,true}, {19,true}, {18,false}, {13,false} },
};

// ---------------------------------------------------------------------------
// Index computation
// ---------------------------------------------------------------------------
// Total vertex count at tessellation level t: 10*t*t + 10 (+ 12 primaries, but
// the primaries are 0..11 and the formula 10t²+2 counts them — here we keep
// the primaries separate for clarity).
//
// Layout:
//   [0..11]                         : 12 primary corners
//   [12 .. 12 + 20*(t-1) - 1]       : edge-interior vertices, 20 edges × (t-1)
//   [12 + 20*(t-1) .. total-1]      : face-interior vertices, 10 faces × (t-1)²

int edge_interior_base(int t) { return 12; }
int face_interior_base(int t) { return 12 + 20*(t-1); }
int total_vertices(int t)     { return 12 + 20*(t-1) + 10*(t-1)*(t-1); }

// Index of the p-th interior point along edge e (0-indexed, p in 0..t-2).
// 'forward' controls which end of the edge p=0 starts from.
// forward=true:  p=0 is near the lo-index corner, p=t-2 is near hi-index corner.
// forward=false: reversed.
int edge_vertex_index(int e, int p, int t, bool forward) {
    int pos = forward ? p : (t-2-p);
    return edge_interior_base(t) + e*(t-1) + pos;
}

// Index of the interior vertex at (row r, col c) within face f, at level t.
// r and c here are the interior row/col, 0-indexed within the interior grid,
// so r in [0..t-2] and c in [0..t-2].
int face_vertex_index(int f, int r_int, int c_int, int t) {
    return face_interior_base(t) + f*(t-1)*(t-1) + r_int*(t-1) + c_int;
}

// Main lookup: given face f, grid row r, grid col c, tessellation level t,
// return the global vertex index. r,c in [0..t].
int vertex_index(int f, int r, int c, int t) {
    bool on_north = (r == 0);
    bool on_south = (r == t);
    bool on_west  = (c == 0);
    bool on_east  = (c == t);

    // --- Corners ---
    if (on_north && on_west)  return FACE_CORNERS[f][0]; // N
    if (on_north && on_east)  return FACE_CORNERS[f][1]; // E
    if (on_south && on_east)  return FACE_CORNERS[f][2]; // S (note: SE = S corner)
    if (on_south && on_west)  return FACE_CORNERS[f][3]; // W (note: SW = W corner)

    // --- Edges ---
    // NE side: row=0, col=1..t-1. Position along edge: col-1 (0..t-2).
    if (on_north) {
        EdgeRef er = FACE_EDGES[f][0]; // NE edge
        return edge_vertex_index(er.edge_id, c-1, t, er.forward);
    }
    // ES side: col=t, row=1..t-1. Position along edge: row-1 (0..t-2).
    if (on_east) {
        EdgeRef er = FACE_EDGES[f][1]; // ES edge
        return edge_vertex_index(er.edge_id, r-1, t, er.forward);
    }
    // SW side: row=t, col=t-1..1. Position along edge: (t-1-col)-1 = t-1-col.
    // col runs t→0 on this side; p=0 is at col=t-1, p=t-2 is at col=1.
    if (on_south) {
        EdgeRef er = FACE_EDGES[f][2]; // SW edge
        int p = (t-1-c) - 1; // col t-1 → p=0, col 1 → p=t-2
        return edge_vertex_index(er.edge_id, p, t, er.forward);
    }
    // WN side: col=0, row=t-1..1. Position along edge: (t-1-row)-1 = t-1-row-... 
    // row runs t→0 on this side; p=0 is at row=t-1, p=t-2 is at row=1.
    if (on_west) {
        EdgeRef er = FACE_EDGES[f][3]; // WN edge
        int p = (t-1-r) - 1; // row t-1 → p=0, row 1 → p=t-2
        return edge_vertex_index(er.edge_id, p, t, er.forward);
    }

    // --- Interior ---
    // Interior row/col: subtract 1 from each (strip the border ring).
    return face_vertex_index(f, r-1, c-1, t);
}

// ---------------------------------------------------------------------------
// Vertex computation
// ---------------------------------------------------------------------------
// Compute the lat/lon of vertex at face f, grid (r,c) at level t,
// by interpolating between the face's four corners using great-circle midpoints.
//
// Strategy: bilinear-like interpolation on the sphere.
// At level 1 the four corners are already known primaries.
// At level t we need the grid point at normalized position (r/t, c/t).
//
// We walk down: at each level we halve the search box.
// But for a flat index function, it's cleaner to interpolate directly:
//   1. Interpolate along the north edge (N to E): fraction c/t → point P_north
//   2. Interpolate along the south edge (W to S): fraction c/t → point P_south
//      Wait — south edge goes S→W in the current convention, so:
//      south edge runs from SE corner to SW corner, fraction c/t from east.
//      Actually simpler: use fraction along each horizontal row.
//   3. Interpolate between P_north and P_south: fraction r/t → final point.
//
// This is the spherical equivalent of bilinear interpolation.
// It matches what midpointCalc was doing iteratively in the original code.
//
// For a vertex at (r, c) in a face with corners N(row0,col0), E(row0,colt),
// S(rowt,colt), W(rowt,col0):
//   top    = great-circle lerp(N, E, c/t)
//   bottom = great-circle lerp(W, S, c/t)
//   result = great-circle lerp(top, bottom, r/t)

pair<double,double> gc_lerp(double lat1, double lon1,
                             double lat2, double lon2,
                             double frac)
{
    // frac=0 → point 1, frac=1 → point 2
    // Implemented as repeated midpoint halving would be expensive;
    // instead use the direct spherical interpolation (slerp on great circle).
    // For small arcs this is numerically equivalent to the midpoint formula.
    double dLon = lon2 - lon1;
    double Bx   = cos(lat2) * cos(dLon);
    double By   = cos(lat2) * sin(dLon);

    // Bearing from point 1 to point 2
    double bearing = atan2(By, cos(lat1) + Bx);
    // Great-circle distance
    double dist = atan2(sqrt((cos(lat1)*sin(dLon))*(cos(lat1)*sin(dLon)) +
                             (cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dLon))*
                             (cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dLon))),
                        sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dLon));

    double d = dist * frac;
    double lat_r = asin(sin(lat1)*cos(d) + cos(lat1)*sin(d)*cos(bearing));
    double lon_r = lon1 + atan2(sin(bearing)*sin(d)*cos(lat1),
                                cos(d) - sin(lat1)*sin(lat_r));

    if (fabs(lat_r) < 1e-10) lat_r = 0.0;
    if (fabs(lon_r) < 1e-10) lon_r = 0.0;
    if (lon_r >  PI) lon_r -= 2.0*PI;
    if (lon_r < -PI) lon_r += 2.0*PI;

    return {lat_r, lon_r};
}

pair<double,double> grid_latlon(int f, int r, int c, int t) {
    // Corners of this face
    double lat_N = PRIMARY[FACE_CORNERS[f][0]].lat;
    double lon_N = PRIMARY[FACE_CORNERS[f][0]].lon;
    double lat_E = PRIMARY[FACE_CORNERS[f][1]].lat;
    double lon_E = PRIMARY[FACE_CORNERS[f][1]].lon;
    double lat_S = PRIMARY[FACE_CORNERS[f][2]].lat;
    double lon_S = PRIMARY[FACE_CORNERS[f][2]].lon;
    double lat_W = PRIMARY[FACE_CORNERS[f][3]].lat;
    double lon_W = PRIMARY[FACE_CORNERS[f][3]].lon;

    double fc = (double)c / (double)t; // east-west fraction (0=west, 1=east)
    double fr = (double)r / (double)t; // north-south fraction (0=north, 1=south)

    // Top row: lerp N→E by fc
    auto [lat_top, lon_top] = gc_lerp(lat_N, lon_N, lat_E, lon_E, fc);
    // Bottom row: lerp W→S by fc
    auto [lat_bot, lon_bot] = gc_lerp(lat_W, lon_W, lat_S, lon_S, fc);
    // Column: lerp top→bottom by fr
    return gc_lerp(lat_top, lon_top, lat_bot, lon_bot, fr);
}

// ---------------------------------------------------------------------------
// Ensure a vertex slot is populated (call planet() if not yet computed).
// ---------------------------------------------------------------------------
void ensure_vertex(vector<Vertex>& VA, int idx,
                   double lat, double lon)
{
    if (VA[idx].computed) return;

    XYZ c = ll_to_xyz(lat, lon);
    planet_out result = planet(tetra[0], tetra[1], tetra[2], tetra[3],
                               c.x, c.y, c.z, Calc_Level);
    VA[idx].lat      = lat;
    VA[idx].lon      = lon;
    VA[idx].height   = result.h * HEIGHT_MOD * RADIUS;
    VA[idx].computed = true;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main()
{
    initialize_vertices();
    Calc_Level = Tessellation_Level + 15;

    int t = Tessellation_Level;
    int total = total_vertices(t);

    cout << "Tessellation level : " << t << endl;
    cout << "Calc_Level         : " << Calc_Level << endl;
    cout << "Total vertex slots : " << total << endl << endl;

    // Allocate vertex array with all slots pre-sized, uncomputed.
    vector<Vertex> VA(total);

    // --- Phase 1: populate the 12 primary corner vertices ---
    for (int i = 0; i < 12; i++) {
        ensure_vertex(VA, i, PRIMARY[i].lat, PRIMARY[i].lon);
    }
    cout << "Primary vertices computed." << endl;

    // --- Phase 2: populate all remaining vertices face by face ---
    // For each face, iterate over its (t+1)×(t+1) grid.
    // Corners and edge vertices will be skipped if already computed
    // (they're shared between faces; the first face to visit them fills them).

    long total_cells = 10LL * (t+1) * (t+1);
    long cell_count  = 0;
    int  report_interval = max(1, (int)(total_cells / 200));

    for (int f = 0; f < 10; f++) {
        for (int r = 0; r <= t; r++) {
            for (int c = 0; c <= t; c++) {
                cell_count++;
                if (cell_count % report_interval == 0) {
                    double pct = 100.0 * cell_count / total_cells;
                    cout << "\r  Populating vertices: face " << f
                         << " | " << fixed << setprecision(1) << pct << "%   " << flush;
                }

                int idx = vertex_index(f, r, c, t);
                if (VA[idx].computed) continue; // shared vertex already done

                auto [lat, lon] = grid_latlon(f, r, c, t);
                ensure_vertex(VA, idx, lat, lon);
            }
        }
    }
    cout << "\r  Vertices populated.                              " << endl;

    // Sanity check: all slots should be computed.
    int uncomputed = 0;
    for (int i = 0; i < total; i++) {
        if (!VA[i].computed) uncomputed++;
    }
    if (uncomputed > 0) {
        cerr << "WARNING: " << uncomputed << " vertex slots were not computed!" << endl;
    }
    cout << "Total vertices computed: " << total - uncomputed << " / " << total << endl << endl;

    // --- Phase 3: build face list ---
    // At tessellation level t, each original face subdivides into t² quads.
    // We enumerate them by iterating over the (t×t) grid of cells per face.
    // Each cell (r,c) → corners at grid positions (r,c),(r,c+1),(r+1,c+1),(r+1,c).
    // In N/E/S/W convention: N=(r,c), E=(r,c+1), S=(r+1,c+1), W=(r+1,c)

    vector<Face> Faces;
    Faces.reserve(10 * t * t);

    for (int f = 0; f < 10; f++) {
        for (int r = 0; r < t; r++) {
            for (int c = 0; c < t; c++) {
                Face face;
                face.v[0] = vertex_index(f, r,   c,   t); // N
                face.v[1] = vertex_index(f, r,   c+1, t); // E
                face.v[2] = vertex_index(f, r+1, c+1, t); // S
                face.v[3] = vertex_index(f, r+1, c,   t); // W
                Faces.push_back(face);
            }
        }
    }

    cout << "Faces built: " << Faces.size() << endl;
    cout << "Expected   : " << 10*t*t << endl << endl;

    // --- Phase 4: OBJ output ---
    ostringstream OFN;
    OFN << "T" << t << "_Tri" << triOrQuad << "_Output.OBJ";
    string OutputFileName = OFN.str();
    ofstream outFile(OutputFileName);

    if (!outFile.is_open()) {
        cerr << "Unable to open output file." << endl;
        return 1;
    }

    outFile << fixed << setprecision(12);
    outFile << "# Icosahedron terrain test - structured grid addressing\n";
    outFile << "# Tessellation level: " << t << "\n";
    outFile << "# Total vertices: " << total << "\n\n";

    // Vertices: convert lat/lon/height to Cartesian XYZ.
    // The sphere has radius RADIUS; height is added radially.
    for (int i = 0; i < total; i++) {
        double r   = RADIUS + VA[i].height;
        double ux  = cos(VA[i].lat) * cos(VA[i].lon);
        double uy  = cos(VA[i].lat) * sin(VA[i].lon);
        double uz  = sin(VA[i].lat);

        double x = r * ux;
        double y = r * uy;
        double z = r * uz;

        if (fabs(x) < 1e-10) x = 0.0;
        if (fabs(y) < 1e-10) y = 0.0;
        if (fabs(z) < 1e-10) z = 0.0;

        // OBJ uses right-handed Y-up; negate Z to match original code's convention.
        outFile << "v " << x << " " << y << " " << -z << "\n";
    }

    outFile << "\n";

    // Faces: OBJ indices are 1-based.
    if (triOrQuad) {
        outFile << "# Faces - Triangles\n";
        for (const Face& face : Faces) {
            // Split quad N/E/S/W into two triangles matching original winding.
            outFile << "f " << face.v[0]+1 << " " << face.v[3]+1 << " " << face.v[1]+1 << "\n";
            outFile << "f " << face.v[3]+1 << " " << face.v[2]+1 << " " << face.v[1]+1 << "\n";
        }
    } else {
        outFile << "# Faces - Quads\n";
        for (const Face& face : Faces) {
            outFile << "f " << face.v[0]+1 << " " << face.v[3]+1
                    << " " << face.v[2]+1 << " " << face.v[1]+1 << "\n";
        }
    }

    outFile.close();
    cout << OutputFileName << " written successfully." << endl;

    return 0;
}
