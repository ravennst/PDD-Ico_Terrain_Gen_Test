#include <iostream>
#include <iomanip>
#include <cmath>
#include "planet_height.hpp"

int main() {
    // Latitude and Longitude in degrees
    double latitude_deg = 67.89;
    double longitude_deg = -12.345;
    double seed = 0.123;
    double PI2 = 3.14159265358979323846;
    
    // Convert to radians
    double latitude_rad = latitude_deg * PI2 / 180.0;
    double longitude_rad = longitude_deg * PI2 / 180.0;

    // Get height, temperature, rainfall, and biome
    double height = planetgen::get_planet_height(latitude_rad, longitude_rad, seed);
    double temperature = planetgen::get_temperature(latitude_rad, height);
    double rainfall = planetgen::get_rainfall(latitude_rad, temperature);
    char biome = planetgen::get_biome(latitude_rad, height);

    // Print results
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Latitude:  " << latitude_deg << "\n";
    std::cout << "Longitude: " << longitude_deg << "\n";
    std::cout << "Seed:      " << seed << "\n";
    std::cout << "Height:    " << height << "\n";
    std::cout << "Temp:      " << temperature << "\n";
    std::cout << "Rainfall:  " << rainfall << "\n";
    std::cout << "Biome:     " << biome << "\n";

    return 0;
}
