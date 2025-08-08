# Tesselation v2 Method

## Introduction

This document outlines the development and refinement of a tesselation algorithm for subdividing an ico-sphere, specifically focusing on improvements made after encountering issues in the first version. The goal is to provide a clear explanation of the challenges faced and the solutions implemented in the second iteration of the method.

## Background

The project focuses on subdividing an icosphere—an icosahedron whose edges are projected outward to form a sphere—by recursively subdividing its triangular faces. The original tesselation algorithm functioned correctly, but its performance degraded at higher subdivision levels. This was due to the need to check if an edge had already been created by searching through an array of vertex pairs for each face. As each tesselation level subdivided every triangle into four, the edge array grew rapidly, causing the algorithm’s runtime to increase exponentially with each additional level of subdivision.

## Issues with the First Version

- Searching the entirety of the Edge array for vertex pairs to see which edges had already been subdivided.
- This was done during sub-division as to not duplicate efforts to create a vertex at a midpint where one already existed.

## Proposed Solution (v2)

- Eliminate the need for an edge array and calculate the entirety of the tesselated ico-sphere at the desired tesselation level from the onset instead of doing it in stages.
- We need to determine how to determine how to calculate each vertex coordinate (Lat/Long) in a pre-defined order that can scale up.
- We can then catalog faces sequentially in a Face Array without the need for an Edge Array for vertex pairs.

## Important formulas and values

- We will need to utilize the Haversine Midpoint formula for some of this:
  - φ = Latitude
  - λ = Longitude
  - B_x = cos φ_2 \* cos Δλ
  - B_y = cos φ_2 \* sin Δλ
  - φ_m = atan2(sin φ_1 + sin φ_2, sqrt((cos φ_1 + B_x)^2) + B_y^2)
  - λ_m = λ_1 + atan2(B_y, cos (φ_1) + B_x)
  - atan2(x,y) = 2 \* arctan(y / (sqrt(x^2 + y^2)  x)) = 2 \* arctan((sqrt( x^2 + y^2) - x) / y)
- θ_c = acos(sqrt(5)/5)
  - This is the central angle of a primary edge of an icosahedron
- θ_long = (2*π)/5
  - Difference in Longitude between each primary edge at the poles.

## Implementation Details

- how to implement

## Results and Discussion

- Summarize the improvements.
- Discuss any remaining challenges.

## Conclusion

- Recap the main points.
- Suggest next steps or further improvements.
