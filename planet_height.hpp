#ifndef PLANET_HEIGHT_H
#define PLANET_HEIGHT_H

#include <cmath>
#include <algorithm>

namespace planetgen {

#define PI 3.14159265358979323846

// Adjustable terrain generation parameters
static const double M = -0.02;     // Base altitude (slightly below sea level)
static const double dd1 = 0.45;    // Weight for altitude difference
static const double POWA = 1.0;    // Exponent for altitude difference
static const double dd2 = 0.035;   // Weight for distance
static const double POW = 0.47;    // Exponent for distance function

// Planetary terrain resolution depth:
// Depth = ceil(log2(2 * Radius / DesiredChordLength)) * 3
// For Radius = 5510.077553 km and DesiredChordLength = 0.001 km (1 meter)
// Depth ≈ 25
static const int Depth = 25;

static const double SQRT3 = std::sqrt(3.0);

// Biome table (from planet.c)
static const char biomes[45][46] = {
  "IIITTTTTGGGGGGGGDDDDDDDDDDDDDDDDDDDDDDDDDDDDD",
  "IIITTTTTGGGGGGGGDDDDGGDSDDSDDDDDDDDDDDDDDDDDD",
  "IITTTTTTTTTBGGGGGGGGGGGSSSSSSDDDDDDDDDDDDDDDD",
  "IITTTTTTTTBBBBBBGGGGGGGSSSSSSSSSWWWWWWWDDDDDD",
  "IITTTTTTTTBBBBBBGGGGGGGSSSSSSSSSSWWWWWWWWWWDD",
  "IIITTTTTTTBBBBBBFGGGGGGSSSSSSSSSSSWWWWWWWWWWW",
  "IIIITTTTTTBBBBBBFFGGGGGSSSSSSSSSSSWWWWWWWWWWW",
  "IIIIITTTTTBBBBBBFFFFGGGSSSSSSSSSSSWWWWWWWWWWW",
  "IIIIITTTTTBBBBBBBFFFFGGGSSSSSSSSSSSWWWWWWWWWW",
  "IIIIIITTTTBBBBBBBFFFFFFGGGSSSSSSSSWWWWWWWWWWW",
  "IIIIIIITTTBBBBBBBFFFFFFFFGGGSSSSSSWWWWWWWWWWW",
  "IIIIIIIITTBBBBBBBFFFFFFFFFFGGSSSSSWWWWWWWWWWW",
  "IIIIIIIIITBBBBBBBFFFFFFFFFFFFFSSSSWWWWWWWWWWW",
  "IIIIIIIIIITBBBBBBFFFFFFFFFFFFFFFSSEEEWWWWWWWW",
  "IIIIIIIIIITBBBBBBFFFFFFFFFFFFFFFFFFEEEEEEWWWW",
  "IIIIIIIIIIIBBBBBBFFFFFFFFFFFFFFFFFFEEEEEEEEWW",
  "IIIIIIIIIIIBBBBBBRFFFFFFFFFFFFFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIBBBBBBRFFFFFFFFFFFFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIIBBBBBRRRFFFFFFFFFFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIIIIBBBRRRRRFFFFFFFFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIIIIIIBRRRRRRRFFFFFFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIIIIIIRRRRRRRRRRFFFFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIIIIIIIRRRRRRRRRRRRFFFFFEEEEEEEEEE",
  "IIIIIIIIIIIIIIIIIIIRRRRRRRRRRRRRFRREEEEEEEEEE",
  "IIIIIIIIIIIIIIIIIIIIIRRRRRRRRRRRRRRRREEEEEEEE",
  "IIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRRRRRRROOEEEEEE",
  "IIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRRROOOOOEEEE",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRRROOOOOOEEE",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRRROOOOOOOEE",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRRROOOOOOOEE",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRRRROOOOOOOOE",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIRRRRROOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIRROOOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIROOOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIROOOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOO",
  "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIOOOOOOO"
};

// Vertex structure representing a point in 3D space with height and seed
struct Vertex {
  double h; // height
  double s; // seed
  double x, y, z;
};

// Tetrahedron vertices (initialized per request)
static Vertex tetra[4];

static double rand2(double p, double q) {
  double r = (p + PI) * (q + PI);
  return 2.0 * (r - static_cast<int>(r)) - 1.0;
}

static double dist2(const Vertex& a, const Vertex& b) {
  double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

// Recursive planet subdivision function
static double planet(Vertex a, Vertex b, Vertex c, Vertex d, double x, double y, double z, int level) {
  if (level <= 0) return 0.25 * (a.h + b.h + c.h + d.h);

  double lab = dist2(a, b);
  if (lab < dist2(a, c)) return planet(a, c, b, d, x, y, z, level);
  if (lab < dist2(a, d)) return planet(a, d, b, c, x, y, z, level);
  if (lab < dist2(b, c)) return planet(b, c, a, d, x, y, z, level);
  if (lab < dist2(b, d)) return planet(b, d, a, c, x, y, z, level);
  if (lab < dist2(c, d)) return planet(c, d, a, b, x, y, z, level);

  Vertex e;
  e.s = rand2(a.s, b.s);
  double es1 = rand2(e.s, e.s);
  double es2 = 0.5 + 0.1 * rand2(es1, es1);
  double es3 = 1.0 - es2;
  e.x = es2 * a.x + es3 * b.x;
  e.y = es2 * a.y + es3 * b.y;
  e.z = es2 * a.z + es3 * b.z;
  e.h = 0.5 * (a.h + b.h) + e.s * dd1 * pow(fabs(a.h - b.h), POWA) + es1 * dd2 * pow(sqrt(lab), POW);

  double det = ((a.x - e.x) * ((c.y - e.y) * (d.z - e.z) - (c.z - e.z) * (d.y - e.y))
              - (a.y - e.y) * ((c.x - e.x) * (d.z - e.z) - (c.z - e.z) * (d.x - e.x))
              + (a.z - e.z) * ((c.x - e.x) * (d.y - e.y) - (c.y - e.y) * (d.x - e.x)));
  double test = ((x - e.x) * ((c.y - e.y) * (d.z - e.z) - (c.z - e.z) * (d.y - e.y))
               - (y - e.y) * ((c.x - e.x) * (d.z - e.z) - (c.z - e.z) * (d.x - e.x))
               + (z - e.z) * ((c.x - e.x) * (d.y - e.y) - (c.y - e.y) * (d.x - e.x)));

  if (det * test > 0)
    return planet(c, d, a, e, x, y, z, level - 1);
  else
    return planet(c, d, b, e, x, y, z, level - 1);
}

// Calls recursive subdivision using preset tetrahedron and depth
static double planet1(double x, double y, double z) {
  return planet(tetra[0], tetra[1], tetra[2], tetra[3], x, y, z, Depth);
}

// Sets up tetrahedron from baseSeed and converts lat/lon to XYZ for sampling
inline double get_planet_height(double latRad, double lonRad, double baseSeed) {
  double r1 = rand2(baseSeed, baseSeed);
  double r2 = rand2(r1, baseSeed);
  double r3 = rand2(r2, r1);
  double r4 = rand2(r3, r2);

  tetra[0] = {M, r1, -SQRT3 - 0.20, -SQRT3 - 0.22, -SQRT3 - 0.23};
  tetra[1] = {M, r2, -SQRT3 - 0.19,  SQRT3 + 0.18,  SQRT3 + 0.17};
  tetra[2] = {M, r3,  SQRT3 + 0.21, -SQRT3 - 0.24,  SQRT3 + 0.15};
  tetra[3] = {M, r4,  SQRT3 + 0.24,  SQRT3 + 0.22, -SQRT3 - 0.25};

  double x = cos(latRad) * cos(lonRad);
  double y = sin(latRad);
  double z = cos(latRad) * sin(lonRad);

  return planet1(x, y, z);
}

inline double get_temperature(double latitude, double altitude) {
  double y = sin(latitude);
  double sun = sqrt(std::max(0.0, 1.0 - y * y));
  return (altitude < 0)
    ? sun / 8.0 + altitude * 0.3
    : sun / 8.0 - altitude * 1.2;
}

inline double get_rainfall(double latitude, double temperature) {
  double y = fabs(sin(latitude));
  double y2 = y - 0.5;
  return std::max(0.0, temperature * 0.65 + 0.1 - 0.011 / (y2 * y2 + 0.1));
}

inline char get_biome(double latitude, double altitude) {
  double temp = get_temperature(latitude, altitude);
  double rain = get_rainfall(latitude, temp);
  int ti = std::min(44, std::max(0, static_cast<int>(rain * 300.0 - 9)));
  int ri = std::min(44, std::max(0, static_cast<int>(temp * 300.0 + 10)));
  return biomes[ti][ri];
}

} // namespace planetgen

#endif // PLANET_HEIGHT_H
