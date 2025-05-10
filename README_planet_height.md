
# `planet_height.h` - Procedural Planet Heightfield and Biome Generator

`planet_height.h` is a single-header, C++ terrain generator using recursive tetrahedral subdivision to simulate planetary surface elevation, climate, and biome classification from latitude, longitude, and seed input.

## 📦 Features

- Deterministic heightfield generation using 3D tetrahedron subdivision
- Procedural temperature and rainfall simulation
- Whittaker-style biome classification
- All logic encapsulated in the `planetgen` namespace
- No external dependencies beyond `<cmath>` and `<algorithm>`

## 🔧 Usage

### 1. **Include the Header**

```cpp
#include "planet_height.h"
```

### 2. **Get Terrain Height**

```cpp
double latRad = ...; // latitude in radians
double lonRad = ...; // longitude in radians
double seed   = 0.12345;

double height = planetgen::get_planet_height(latRad, lonRad, seed);
```

### 3. **Get Climate Data**

```cpp
double temp  = planetgen::get_temperature(latRad, height);
double rain  = planetgen::get_rainfall(latRad, temp);
char biomeId = planetgen::get_biome(latRad, height);
```

## 🌍 API Reference

### `double get_planet_height(double latRad, double lonRad, double baseSeed)`
Returns the height (altitude) at the given coordinates on the unit sphere.

| Parameter   | Type    | Description                            |
|-------------|---------|----------------------------------------|
| `latRad`    | double  | Latitude in radians                    |
| `lonRad`    | double  | Longitude in radians                   |
| `baseSeed`  | double  | Seed to ensure repeatable terrain      |
| **Returns** | double  | Terrain elevation at that point        |

### `double get_temperature(double latitude, double altitude)`
Returns the temperature at the given latitude and altitude.

| Parameter    | Type    | Description                       |
|--------------|---------|-----------------------------------|
| `latitude`   | double  | Latitude in radians               |
| `altitude`   | double  | Elevation from `get_planet_height`|
| **Returns**  | double  | Simulated temperature             |

### `double get_rainfall(double latitude, double temperature)`
Estimates rainfall level based on location and temperature.

| Parameter    | Type    | Description              |
|--------------|---------|--------------------------|
| `latitude`   | double  | Latitude in radians      |
| `temperature`| double  | Temperature from previous |
| **Returns**  | double  | Rainfall level (0.0–1.0+) |

### `char get_biome(double latitude, double altitude)`
Returns a biome classification code based on the Whittaker diagram.

| Output | Biome Type                    |
|--------|-------------------------------|
| `T`    | Tundra                        |
| `G`    | Grasslands                    |
| `B`    | Taiga / Boreal Forest         |
| `D`    | Desert                        |
| `S`    | Savanna                       |
| `F`    | Temperate Forest              |
| `R`    | Temperate Rainforest          |
| `W`    | Xeric Shrubland / Dry Forest  |
| `E`    | Tropical Dry Forest           |
| `O`    | Tropical Rainforest           |
| `I`    | Icecap                        |

## 🧪 Example

```cpp
#include <iostream>
#include "planet_height.h"

int main() {
    double lat = 0.5;      // Radians
    double lon = 1.0;      // Radians
    double seed = 42.0;

    double height = planetgen::get_planet_height(lat, lon, seed);
    double temp = planetgen::get_temperature(lat, height);
    double rain = planetgen::get_rainfall(lat, temp);
    char biome = planetgen::get_biome(lat, height);

    std::cout << "Height: " << height << "\n";
    std::cout << "Temp:   " << temp << "\n";
    std::cout << "Rain:   " << rain << "\n";
    std::cout << "Biome:  " << biome << "\n";

    return 0;
}
```

## 🧭 Notes

- Angles must be provided in **radians**, not degrees.
- The system uses recursive midpoint displacement inside a tetrahedron for fractal terrain.
- The output is deterministic per `(lat, lon, seed)`.

## 📜 License

Portions adapted from Torben Ægidius Mogensen’s `planet.c` (public domain).
