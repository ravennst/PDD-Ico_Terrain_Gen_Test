/* planet_height.h - self-contained version of Torben AE. Mogensen's terrain height algorithm Simplified for header-only use in voxel engines and spherical planet terrain queries. Input: latitude and longitude (in radians), and a seed Output: height value in range [-1, +1] */

#ifndef PLANET_HEIGHT_H 
#define PLANET_HEIGHT_H

#include <cmath>

namespace planetgen {

static constexpr double PI = 3.141592653589793;

struct vertex { double h; // altitude double s; // seed double x, y, z; // position };

// Internal tetrahedron static vertex tetra[4];

// Configurable planet parameters static double M = -0.02;     // baseline height static double dd1 = 0.45;    // weight for altitude difference static double POWA = 1.0;    // power for altitude difference static double dd2 = 0.035;   // weight for distance static double POW = 0.47;    // power for distance static int Depth = 10;       // subdivision depth

static double rand2(double p, double q) { double r = (p + PI) * (q + PI); return 2.0 * (r - (int)r) - 1.0; }

static double dist2(vertex a, vertex b) { double dx = a.x - b.x; double dy = a.y - b.y; double dz = a.z - b.z; return dx * dx + dy * dy + dz * dz; }

static vertex ssa, ssb, ssc, ssd;

static double planet(vertex a, vertex b, vertex c, vertex d, double x, double y, double z, int level) { if (level == 0) { return 0.25 * (a.h + b.h + c.h + d.h); }

// Always split longest edge ab double lab = dist2(a, b); double lac = dist2(a, c); double lad = dist2(a, d); double lbc = dist2(b, c); double lbd = dist2(b, d); double lcd = dist2(c, d); double maxl = lab; if (lac > maxl) return planet(a, c, b, d, x, y, z, level); if (lad > maxl) return planet(a, d, b, c, x, y, z, level); if (lbc > maxl) return planet(b, c, a, d, x, y, z, level); if (lbd > maxl) return planet(b, d, a, c, x, y, z, level); if (lcd > maxl) return planet(c, d, a, b, x, y, z, level);

vertex e; e.s = rand2(a.s, b.s); double es1 = rand2(e.s, e.s); double es2 = 0.5 + 0.1 * rand2(es1, es1); double es3 = 1.0 - es2; if (a.s < b.s) { e.x = es2 * a.x + es3 * b.x; e.y = es2 * a.y + es3 * b.y; e.z = es2 * a.z + es3 * b.z; } else { e.x = es3 * a.x + es2 * b.x; e.y = es3 * a.y + es2 * b.y; e.z = es3 * a.z + es2 * b.z; }

double l = std::sqrt(e.x * e.x + e.y * e.y + e.z * e.z); if (l == 0.0) l = 1.0; e.h = 0.5 * (a.h + b.h) + e.s * dd1 * std::pow(std::abs(a.h - b.h), POWA) + es1 * dd2 * std::pow(std::sqrt(lab), POW);

// Choose sub-tetrahedron double ax = a.x - e.x, ay = a.y - e.y, az = a.z - e.z; double cx = c.x - e.x, cy = c.y - e.y, cz = c.z - e.z; double dx = d.x - e.x, dy = d.y - e.y, dz = d.z - e.z; double px = x - e.x, py = y - e.y, pz = z - e.z; double det = ax*(cydz - czdy) + ay*(czdx - cxdz) + az*(cxdy - cydx); double detp = px*(cydz - czdy) + py*(czdx - cxdz) + pz*(cxdy - cydx); if (det * detp > 0.0) return planet(c, d, a, e, x, y, z, level - 1); else return planet(c, d, b, e, x, y, z, level - 1); }

static double planet1(double x, double y, double z) { return planet(tetra[0], tetra[1], tetra[2], tetra[3], x, y, z, Depth); }

static void init_tetrahedron(double seed) { double r1 = seed; r1 = rand2(r1, r1); double r2 = rand2(r1, r1); double r3 = rand2(r1, r2); double r4 = rand2(r2, r3);

double R = std::sqrt(3.0); tetra[0] = {M, r1, -R - 0.20, -R - 0.22, -R - 0.23}; tetra[1] = {M, r2, -R - 0.19, +R + 0.18, +R + 0.17}; tetra[2] = {M, r3, +R + 0.21, -R - 0.24, +R + 0.15}; tetra[3] = {M, r4, +R + 0.24, +R + 0.22, -R - 0.25}; }

static double get_planet_height(double lat, double lon, double seed) { init_tetrahedron(seed); double x = std::cos(lat) * std::cos(lon); double y = std::sin(lat); double z = std::cos(lat) * std::sin(lon); return planet1(x, y, z); }

} // namespace planetgen

#endif // PLANET_HEIGHT_H

