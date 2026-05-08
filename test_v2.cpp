/*
 * test_v2.cpp - Icosahedron Tessellation: Array.txt Structured Addressing
 *
 * Implements the vertex array layout defined in Documentation/Array.txt.
 *
 * ==========================================================================
 * ADDRESSING SCHEME (from Array.txt)
 * ==========================================================================
 *
 * The vertex array is a single flat array of size 10*t*t + 10.
 * All primary vertices sit at indices that are multiples of t:
 *
 *   [0]   = North pole
 *   [1t]  = Northern ring vertex 1   (Long = 0)
 *   [2t]  = Southern ring vertex 1   (Long = pi/5)
 *   [3t]  = South pole  (canonical)
 *   [4t]  = Northern ring vertex 2   (Long = 2pi/5)
 *   [5t]  = Southern ring vertex 2   (Long = 3pi/5)
 *   [6t]  = null (South pole duplicate -> points to [3t])
 *   [7t]  = Northern ring vertex 3   (Long = 4pi/5)
 *   [8t]  = Southern ring vertex 3   (Long = pi)
 *   [9t]  = null (South pole duplicate)
 *   [10t] = Northern ring vertex 4   (Long = 6pi/5)
 *   [11t] = Southern ring vertex 4   (Long = 7pi/5)
 *   [12t] = null (South pole duplicate)
 *   [13t] = Northern ring vertex 5   (Long = 8pi/5)
 *   [14t] = Southern ring vertex 5   (Long = 9pi/5)
 *   [15t] = null (South pole duplicate)
 *
 * Edge-interior vertices fill the gaps between primary anchors:
 *   e1:  slots [1,      1t-1]   (between [0]  and [1t])
 *   e2:  slots [1t+1,   2t-1]   (between [1t] and [2t])
 *   ...
 *   e15: slots [14t+1,  15t-1]
 *   e16: slots [15t+1,  16t-1]  (geometric: 1t -> 14t)
 *   e17: slots [16t+1,  17t-1]  (geometric: 4t -> 2t)
 *   e18: slots [17t+1,  18t-1]  (geometric: 7t -> 5t)
 *   e19: slots [18t+1,  19t-1]  (geometric: 10t -> 8t)
 *   e20: slots [19t+1,  20t-1]  (geometric: 13t -> 11t)
 *
 * Null slots:
 *   6t, 9t, 12t, 15t  = south pole duplicates (mirror 3t)
 *   16t, 17t, 18t, 19t = pure padding nulls (no geometric meaning;
 *                        keep t-offset arithmetic consistent for e16-e20)
 *   20t = first face interior vertex of f0 (real data, not null)
 *
 * Face-interior vertices start at index 20t:
 *   f0:  slots [20t,               20t + (t-1)^2 - 1]
 *   fn:  slots [20t + n*(t-1)^2,   20t + (n+1)*(t-1)^2 - 1]
 *
 * Within each face, interior points are stored in a (t-1)x(t-1) grid,
 * row-major from the NE side inward. Row ln (0..t-2), col lp (0..t-2).
 * Slot = face_base + ln*(t-1) + lp
 *
 * ==========================================================================
 * HAVERSINE INTERMEDIATE POINT FORMULA (Array.txt)
 * ==========================================================================
 *
 *   delta = angular distance between point 1 and point 2
 *   A = sin((1-f)*delta) / sin(delta)
 *   B = sin(f*delta) / sin(delta)
 *   x = A*cos(lat1)*cos(lon1) + B*cos(lat2)*cos(lon2)
 *   y = A*cos(lat1)*sin(lon1) + B*cos(lat2)*sin(lon2)
 *   z = A*sin(lat1)           + B*sin(lat2)
 *   lat_i = atan2(z, sqrt(x^2+y^2))
 *   lon_i = atan2(y, x)
 *
 *   f = p/t for the p-th point along an edge of t divisions (p = 1..t-1)
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

extern "C" {
  #include "libs/Planet/planet.h"
}

double rseed    = 0.21;
double M        = 0.021;
planet_vertex tetra[4];

using namespace std;

const double PI         = 3.141592653589793;
const double RADIUS     = 1.0;
const double HEIGHT_MOD = 1.0;

int  Tessellation_Level = 300;
bool triOrQuad          = true; // if true, output triangles; if false, output quads (one per face grid cell)
int  Calc_Level; // internal calc level for planet() calls; set to Tessellation_Level+15 for extra detail

// ---------------------------------------------------------------------------
// Vertex and face structs
// ---------------------------------------------------------------------------
struct Vertex {
    double lat = 0.0, lon = 0.0, height = 0.0;
    bool computed = false;
};

struct Face { int v[4]; }; // N, E, S, W

// ---------------------------------------------------------------------------
// Primary vertex definitions
// primary_anchor[n]: the vertex whose array slot is n*t.
// Anchors 6,9,12,15 are south-pole duplicates (is_null=true).
// ---------------------------------------------------------------------------
const double ICO_C = acos(sqrt(5.0)/5.0);
const double X1    = PI/2.0 - ICO_C;

struct LatLon { double lat, lon; bool is_null; };

const LatLon PRIMARY[16] = {
    {  PI/2.0,              0.0, false }, //  0: North pole
    {  X1,                  0.0, false }, //  1: N-ring 1   Long=0
    { -X1,             PI/5.0,   false }, //  2: S-ring 1   Long=pi/5
    { -PI/2.0,              0.0, false }, //  3: South pole (canonical)
    {  X1,       (2.0*PI)/5.0,   false }, //  4: N-ring 2
    { -X1,    (3.0*PI)/5.0,      false }, //  5: S-ring 2
    {  0.0,                 0.0, true  }, //  6: null -> [3t]
    {  X1,       (4.0*PI)/5.0,   false }, //  7: N-ring 3
    { -X1,    (5.0*PI)/5.0,      false }, //  8: S-ring 3
    {  0.0,                 0.0, true  }, //  9: null -> [3t]
    {  X1,       (6.0*PI)/5.0,   false }, // 10: N-ring 4
    { -X1,    (7.0*PI)/5.0,      false }, // 11: S-ring 4
    {  0.0,                 0.0, true  }, // 12: null -> [3t]
    {  X1,       (8.0*PI)/5.0,   false }, // 13: N-ring 5
    { -X1,    (9.0*PI)/5.0,      false }, // 14: S-ring 5
    {  0.0,                 0.0, true  }, // 15: null -> [3t]
};

int resolve_anchor(int n) {
    if (n==6||n==9||n==12||n==15) return 3;
    return n;
}

// ---------------------------------------------------------------------------
// 20 edges defined in Array.txt.
//
// sa        = geometric start anchor (haversine interpolation)
// ea        = geometric end anchor
// Slot ranges are 20 sequential t-sized windows: edge e (0-indexed)
// occupies slots e*t+1 .. e*t+(t-1). The geometric start/end vertices
// are purely for haversine interpolation and are independent of slot math.
//
// Null padding slots: 6t,9t,12t,15t mirror south pole (3t).
//                     16t,17t,18t,19t are pure padding nulls.
//                     20t is the first face interior vertex (real data).
// ---------------------------------------------------------------------------
struct EdgeDef { int sa, ea; };
const EdgeDef EDGES[20] = {
    {0,1}, // e1:  slots 1..1t-1,       interp 0   -> 1t
    {1,2}, // e2:  slots 1t+1..2t-1,    interp 1t  -> 2t
    {2,3}, // e3:  slots 2t+1..3t-1,    interp 2t  -> 3t
    {0,4}, // e4:  slots 3t+1..4t-1,    interp 0   -> 4t
    {4,5}, // e5:  slots 4t+1..5t-1,    interp 4t  -> 5t
    {5,3}, // e6:  slots 5t+1..6t-1,    interp 5t  -> 3t
    {0,7}, // e7:  slots 6t+1..7t-1,    interp 0   -> 7t
    {7,8}, // e8:  slots 7t+1..8t-1,    interp 7t  -> 8t
    {8,3}, // e9:  slots 8t+1..9t-1,    interp 8t  -> 3t
    {0,10}, // e10: slots 9t+1..10t-1,   interp 0   -> 10t
    {10,11}, // e11: slots 10t+1..11t-1,  interp 10t -> 11t
    {11,3}, // e12: slots 11t+1..12t-1,  interp 11t -> 3t
    {0,13}, // e13: slots 12t+1..13t-1,  interp 0   -> 13t
    {13,14}, // e14: slots 13t+1..14t-1,  interp 13t -> 14t
    {14,3}, // e15: slots 14t+1..15t-1,  interp 14t -> 3t
    {1,14}, // e16: slots 15t+1..16t-1,  interp 1t  -> 14t
    {4,2}, // e17: slots 16t+1..17t-1,  interp 4t  -> 2t
    {7,5}, // e18: slots 17t+1..18t-1,  interp 7t  -> 5t
    {10,8}, // e19: slots 18t+1..19t-1,  interp 10t -> 8t
    {13,11}, // e20: slots 19t+1..20t-1,  interp 13t -> 11t
};
// ---------------------------------------------------------------------------
// Face topology: for each face, N/E/S/W corner anchors and NE/ES/SW/WN edges.
// edge_fwd[i]=true means the face traverses that edge in the same direction
// as start_anchor->end_anchor (i.e., p=0 is nearest the N or W corner).
// ---------------------------------------------------------------------------
struct FaceDef {
    int  corner[4];   // N, E, S, W anchor indices
    int  edge[4];     // NE, ES, SW, WN edge indices (0-based)
    bool fwd[4];      // traversal direction
};

// Derivation from Array.txt diagram:
// Face f0: N=0, E=4t, S=2t, W=1t  -> anchors 0,4,2,1
//   NE=e4(0->4t) fwd, ES=e17(4t->2t) fwd,
//   SW=e2(1t->2t) rev (face goes S->W = 2t->1t = end->start of e2),
//   WN=e1(0->1t) rev (face goes W->N = 1t->0 = end->start of e1)
// Face f1: N=1t, E=2t, S=3t, W=14t -> anchors 1,2,3,14
//   NE=e2 fwd, ES=e3 fwd, SW=e15(14t->3t) rev, WN=e16(1t->14t) rev
// Face f2: N=0, E=7t, S=5t, W=4t -> anchors 0,7,5,4
//   NE=e7 fwd, ES=e18(7t->5t) fwd, SW=e5(4t->5t) rev, WN=e4(0->4t) rev
// Face f3: N=4t, E=5t, S=3t, W=2t -> anchors 4,5,3,2
//   NE=e5 fwd, ES=e6(5t->3t) fwd, SW=e3(2t->3t) rev, WN=e17(4t->2t) rev
// Face f4: N=0, E=10t, S=8t, W=7t -> anchors 0,10,8,7
//   NE=e10 fwd, ES=e19(10t->8t) fwd, SW=e8(7t->8t) rev, WN=e7(0->7t) rev
// Face f5: N=7t, E=8t, S=3t, W=5t -> anchors 7,8,3,5
//   NE=e8 fwd, ES=e9(8t->3t) fwd, SW=e6(5t->3t) rev, WN=e18(7t->5t) rev
// Face f6: N=0, E=13t, S=11t, W=10t -> anchors 0,13,11,10
//   NE=e13 fwd, ES=e20(13t->11t) fwd, SW=e11(10t->11t) rev, WN=e10(0->10t) rev
// Face f7: N=10t, E=11t, S=3t, W=8t -> anchors 10,11,3,8
//   NE=e11 fwd, ES=e12(11t->3t) fwd, SW=e9(8t->3t) rev, WN=e19(10t->8t) rev
// Face f8: N=0, E=1t, S=14t, W=13t -> anchors 0,1,14,13
//   NE=e1 fwd, ES=e16(1t->14t) fwd, SW=e14(13t->14t) rev, WN=e13(0->13t) rev
// Face f9: N=13t, E=14t, S=3t, W=11t -> anchors 13,14,3,11
//   NE=e14 fwd, ES=e15(14t->3t) fwd, SW=e12(11t->3t) rev, WN=e20(13t->11t) rev

const FaceDef FACES[10] = {
    {{0, 4, 2, 1},  { 3,16, 1, 0},  {true, true,false,false}}, // f0
    {{1, 2, 3,14},  { 1, 2,14,15},  {true, true,false,false}}, // f1
    {{0, 7, 5, 4},  { 6,17, 4, 3},  {true, true,false,false}}, // f2
    {{4, 5, 3, 2},  { 4, 5, 2,16},  {true, true,false,false}}, // f3
    {{0,10, 8, 7},  { 9,18, 7, 6},  {true, true,false,false}}, // f4
    {{7, 8, 3, 5},  { 7, 8, 5,17},  {true, true,false,false}}, // f5
    {{0,13,11,10},  {12,19,10, 9},  {true, true,false,false}}, // f6
    {{10,11, 3, 8}, {10,11, 8,18},  {true, true,false,false}}, // f7
    {{0, 1,14,13},  { 0,15,13,12},  {true, true,false,false}}, // f8
    {{13,14, 3,11}, {13,14,11,19},  {true, true,false,false}}, // f9
};

// ---------------------------------------------------------------------------
// Index helpers
// ---------------------------------------------------------------------------
int total_vertices(int t) { return 10*t*t + 10; }
int face_base(int f, int t) { return 20*t + f*(t-1)*(t-1); }

// Slot for the p-th interior point of edge e.
// e is 0-indexed. p is 1-indexed (1..t-1), matching f=p/t from Array.txt.
// Slot ranges are 20 sequential t-sized windows: edge e occupies e*t+1..e*t+(t-1).
int edge_slot(int e, int p, int t) {
    return e * t + p;
}

// ---------------------------------------------------------------------------
// Haversine intermediate point (Array.txt formula)
// ---------------------------------------------------------------------------
pair<double,double> hav_interp(double lat1, double lon1,
                                double lat2, double lon2,
                                double f)
{
    double sin_dlat = sin((lat2-lat1)/2.0);
    double sin_dlon = sin((lon2-lon1)/2.0);
    double a = sin_dlat*sin_dlat + cos(lat1)*cos(lat2)*sin_dlon*sin_dlon;
    double delta = 2.0*atan2(sqrt(a), sqrt(1.0-a));

    if (delta < 1e-12) return {lat1, lon1};

    double A = sin((1.0-f)*delta) / sin(delta);
    double B = sin(f*delta)       / sin(delta);

    double x = A*cos(lat1)*cos(lon1) + B*cos(lat2)*cos(lon2);
    double y = A*cos(lat1)*sin(lon1) + B*cos(lat2)*sin(lon2);
    double z = A*sin(lat1)           + B*sin(lat2);

    double lat_i = atan2(z, sqrt(x*x+y*y));
    double lon_i = atan2(y, x);

    if (fabs(lat_i)<1e-10) lat_i=0.0;
    if (fabs(lon_i)<1e-10) lon_i=0.0;
    if (lon_i> PI) lon_i-=2.0*PI;
    if (lon_i<-PI) lon_i+=2.0*PI;

    return {lat_i, lon_i};
}

// ---------------------------------------------------------------------------
// Compute and store a vertex (skip if already done)
// ---------------------------------------------------------------------------
void compute_vertex(vector<Vertex>& VA, int slot, double lat, double lon) {
    if (VA[slot].computed) return;
    double cx = cos(lat)*sin(lon);
    double cy = sin(lat);
    double cz = cos(lat)*cos(lon);
    planet_out r = planet(tetra[0],tetra[1],tetra[2],tetra[3], cx,cy,cz, Calc_Level);
    VA[slot] = {lat, lon, r.h*HEIGHT_MOD*RADIUS, true};
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main()
{
    initialize_vertices();
    int t  = Tessellation_Level;
    Calc_Level = min(t + 15, 30);

    int total = total_vertices(t);
    cout << "Tessellation level : " << t          << "\n";
    cout << "Calc_Level         : " << Calc_Level  << "\n";
    cout << "Total vertex slots : " << total       << "\n\n";

    vector<Vertex> VA(total);

    // ------------------------------------------------------------------
    // Phase 1: Primary vertices at slots n*t
    // ------------------------------------------------------------------
    for (int n = 0; n < 16; n++) {
        if (PRIMARY[n].is_null) continue;
        compute_vertex(VA, n*t, PRIMARY[n].lat, PRIMARY[n].lon);
    }
    // Mirror south-pole data into null duplicate slots
    for (int n : {6,9,12,15}) VA[n*t] = VA[3*t];
    // Mark padding null slots (16t-19t) as computed with zero data
    // so the computed flag check never fires on them accidentally.
    for (int n : {16,17,18,19}) VA[n*t].computed = true;

    cout << "Primary vertices done.\n";

    // ------------------------------------------------------------------
    // Phase 2: Edge-interior vertices
    // Slots: e*t+p where p=1..t-1, fraction f=p/t
    // Geometry: interpolate between sa and ea (ea resolved through nulls).
    // ------------------------------------------------------------------
    for (int e = 0; e < 20; e++) {
        int sa = EDGES[e].sa;
        int ea = resolve_anchor(EDGES[e].ea);
        double lat1=PRIMARY[sa].lat, lon1=PRIMARY[sa].lon;
        double lat2=PRIMARY[ea].lat, lon2=PRIMARY[ea].lon;
        for (int p = 1; p < t; p++) {
            double f = (double)p/(double)t;
            auto [lat,lon] = hav_interp(lat1,lon1,lat2,lon2,f);
            compute_vertex(VA, edge_slot(e,p,t), lat, lon);
        }
    }
    cout << "Edge vertices done.\n";

    // ------------------------------------------------------------------
    // Phase 3: Face-interior vertices
    //
    // For face f, interior line ln (0..t-2) runs from the ln-th interior
    // point of the WN edge to the ln-th interior point of the ES edge.
    // (WN and ES are the left and right sides as you look at the face with
    // N at top.)
    //
    // Along that line, interior points lp (0..t-2) are at fraction
    // f=(lp+1)/t from the WN endpoint toward the ES endpoint.
    // ------------------------------------------------------------------
    long total_fc = 10L*(t-1)*(t-1);
    long fc = 0;
    int  rpt = max(1L, total_fc/200);

    for (int f = 0; f < 10; f++) {
        const FaceDef& fd = FACES[f];
        int wn_e = fd.edge[3];  bool wn_fwd = fd.fwd[3];
        int es_e = fd.edge[1];  bool es_fwd = fd.fwd[1];

        int wn_sa = EDGES[wn_e].sa, wn_ea = resolve_anchor(EDGES[wn_e].ea);
        int es_sa = EDGES[es_e].sa, es_ea = resolve_anchor(EDGES[es_e].ea);

        double wn_lat1=PRIMARY[wn_sa].lat, wn_lon1=PRIMARY[wn_sa].lon;
        double wn_lat2=PRIMARY[wn_ea].lat, wn_lon2=PRIMARY[wn_ea].lon;
        double es_lat1=PRIMARY[es_sa].lat, es_lon1=PRIMARY[es_sa].lon;
        double es_lat2=PRIMARY[es_ea].lat, es_lon2=PRIMARY[es_ea].lon;

        int fb = face_base(f,t);

        for (int ln = 0; ln < t-1; ln++) {
            // Fraction along WN/ES edges for this line's boundary endpoints.
            // ln=0 is the line nearest the N corner (p=0 on the edge = f=1/t).
            double f_edge = (double)(ln+1)/(double)t;
            double wn_f = wn_fwd ? (1.0-f_edge) : f_edge;
            double es_f = es_fwd ? f_edge : (1.0-f_edge);

            auto [wn_lat,wn_lon] = hav_interp(wn_lat1,wn_lon1,wn_lat2,wn_lon2, wn_f);
            auto [es_lat,es_lon] = hav_interp(es_lat1,es_lon1,es_lat2,es_lon2, es_f);

            for (int lp = 0; lp < t-1; lp++) {
                fc++;
                if (fc%rpt==0)
                    cout << "\r  Face interiors: f" << f << " ln" << ln
                         << " | " << fixed << setprecision(1)
                         << (100.0*fc/total_fc) << "%    " << flush;

                double f_line = (double)(lp+1)/(double)t;
                auto [lat,lon] = hav_interp(wn_lat,wn_lon,es_lat,es_lon, f_line);
                compute_vertex(VA, fb + ln*(t-1) + lp, lat, lon);
            }
        }
    }
    cout << "\r  Face interiors done.                          \n";

    int uncomputed=0;
    for (int i=0;i<total;i++) if (!VA[i].computed) uncomputed++;
    if (uncomputed) cerr << "WARNING: " << uncomputed << " uncomputed slots!\n";
    cout << "Vertices: " << total-uncomputed << " / " << total << "\n\n";

    // ------------------------------------------------------------------
    // Phase 4: Build face quad list
    //
    // grid_slot(f, row, col) returns the VA index for the grid point at
    // (row, col) within face f, where row,col in [0..t].
    // ------------------------------------------------------------------
    auto grid_slot = [&](int f, int row, int col) -> int {
        const FaceDef& fd = FACES[f];
        bool nn=(row==0), ss=(row==t), ww=(col==0), ee=(col==t);

        // Corners
        if (nn&&ww) return fd.corner[0]*t; // N
        if (nn&&ee) return fd.corner[1]*t; // E
        if (ss&&ee) return resolve_anchor(fd.corner[2])*t; // S
        if (ss&&ww) return resolve_anchor(fd.corner[3])*t; // W

        // NE edge (row=0, col 1..t-1): p steps from N corner
        if (nn) {
            int e=fd.edge[0]; bool fwd=fd.fwd[0];
            int p = fwd ? col : (t-col);
            return edge_slot(e,p,t);
        }
        // ES edge (col=t, row 1..t-1)
        if (ee) {
            int e=fd.edge[1]; bool fwd=fd.fwd[1];
            int p = fwd ? row : (t-row);
            return edge_slot(e,p,t);
        }
        // SW edge (row=t, col t-1..1)
        if (ss) {
            int e=fd.edge[2]; bool fwd=fd.fwd[2];
            // col=t-1 is one step from the S/E corner; col=1 is one step from W
            // fwd=false (typical): edge goes end->start = S->W direction in face
            // p counts from edge's start_anchor:
            //   fwd=false means face's col=1 is near start_anchor, col=t-1 near end_anchor
            int p = fwd ? (t-col) : col;
            return edge_slot(e,p,t);
        }
        // WN edge (col=0, row t-1..1)
        if (ww) {
            int e=fd.edge[3]; bool fwd=fd.fwd[3];
            int p = fwd ? (t-row) : row;
            return edge_slot(e,p,t);
        }

        // Interior
        return face_base(f,t) + (row-1)*(t-1) + (col-1);
    };

    vector<Face> Faces;
    Faces.reserve(10*t*t);
    for (int f=0;f<10;f++)
        for (int r=0;r<t;r++)
            for (int c=0;c<t;c++)
                Faces.push_back({ {grid_slot(f,r,c), grid_slot(f,r,c+1),
                                   grid_slot(f,r+1,c+1), grid_slot(f,r+1,c)} });

    cout << "Faces built: " << Faces.size() << " (expected " << 10*t*t << ")\n\n";

    // ------------------------------------------------------------------
    // Phase 5: OBJ output
    // ------------------------------------------------------------------
    ostringstream OFN;
    OFN << "T" << t << "_Tri" << triOrQuad << "_Output.OBJ";
    ofstream out(OFN.str());
    if (!out.is_open()) { cerr << "Cannot open output file.\n"; return 1; }

    out << fixed << setprecision(12);
    out << "# Icosahedron terrain - Array.txt structured addressing\n";
    out << "# Tessellation level: " << t << "\n";
    out << "# Total vertices: " << total << "\n\n";

    for (int i=0;i<total;i++) {
        double r = RADIUS + VA[i].height;
        double x = r*cos(VA[i].lat)*cos(VA[i].lon);
        double y = r*cos(VA[i].lat)*sin(VA[i].lon);
        double z = r*sin(VA[i].lat);
        if (fabs(x)<1e-10) x=0.0;
        if (fabs(y)<1e-10) y=0.0;
        if (fabs(z)<1e-10) z=0.0;
        out << "v " << x << " " << -z << " " << -y << "\n";
    }
    out << "\n";

    if (triOrQuad) {
        out << "# Faces - Triangles\n";
        for (auto& face : Faces) {
            out << "f " << face.v[0]+1 << " " << face.v[3]+1 << " " << face.v[1]+1 << "\n";
            out << "f " << face.v[3]+1 << " " << face.v[2]+1 << " " << face.v[1]+1 << "\n";
        }
    } else {
        out << "# Faces - Quads\n";
        for (auto& face : Faces) {
            out << "f " << face.v[0]+1 << " " << face.v[3]+1
                << " " << face.v[2]+1 << " " << face.v[1]+1 << "\n";
        }
    }

    out.close();
    cout << OFN.str() << " written successfully.\n";
    return 0;
}
