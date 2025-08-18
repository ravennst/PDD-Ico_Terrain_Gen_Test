<!-- markdownlint-disable MD033 -->

# <span style="background-color: blue;">Tesselation v2 Method</span>

## <span style="background-color: blue;">Introduction</span>

This document outlines the development and refinement of a tesselation algorithm for subdividing an ico-sphere, specifically focusing on improvements made after encountering issues in the first version. The goal is to provide a clear explanation of the challenges faced and the solutions implemented in the second iteration of the method.

## <span style="background-color: blue;">Background</span>

The project focuses on subdividing an icosphere—an icosahedron whose edges are projected outward to form a sphere—by recursively subdividing its triangular faces. The original tesselation algorithm functioned correctly, but its performance degraded at higher subdivision levels. This was due to the need to check if an edge had already been created by searching through an array of vertex pairs for each face. As each tesselation level subdivided every triangle into four, the edge array grew rapidly, causing the algorithm’s runtime to increase exponentially with each additional level of subdivision.

## <span style="background-color: blue;">Issues with the First Version</span>

- Searching the entirety of the Edge array for vertex pairs to see which edges had already been subdivided.
- This was done during sub-division as to not duplicate efforts to create a vertex at a midpint where one already existed.

## <span style="background-color: blue;">Proposed Solution (v2)</span>

- Eliminate the need for an edge array and calculate the entirety of the tesselated ico-sphere at the desired tesselation level from the onset instead of doing it in stages.
- We need to determine how to determine how to calculate each vertex coordinate (Lat/Long) in a pre-defined order that can scale up.
- We can then catalog faces sequentially in a Face Array without the need for an Edge Array for vertex pairs.

## <span style="background-color: blue;">Important formulas and values</span>

- We will need to utilize the Haversine Midpoint formula for some of this:
  - &phi; = Latitude
  - &lambda; = Longitude
  - B<sub>x</sub> = cos &phi;<sub>2</sub> \* cos &Delta;&lambda;
  - B<sub>y</sub> = cos &phi;<sub>2</sub> \* sin &Delta;&lambda;
  - &phi;<sub>mid</sub> = atan2( sin &phi;<sub>1</sub> + sin &phi;<sub>2</sub>, &radic;(( cos &phi;<sub>1</sub> + B<sub>x</sub> )<sup>2</sup> ) + B<sub>y</sub><sup>2</sup> )
  - &lambda;<sub>mid</sub> = &lambda;<sub>1</sub> + atan2( B<sub>y</sub>, cos &phi;<sub>1</sub> + B<sub>x</sub> )
  - atan2( x, y ) = 2 \* arctan( y / (&radic;( x<sup>2</sup> + y<sup>2</sup>)  x)) = 2 \* arctan(( &radic;( x<sup>2</sup> + y<sup>2</sup>) - x) / y )
- &theta;<sub>c</sub> = acos(&radic;(5)/5)
  - This is the central angle of a primary edge of an icosahedron
- &theta;<sub>long</sub> = 2π/5
  - Difference in Longitude between each primary edge at the poles.

## <span style="background-color: blue;">Implementation Details</span>

- The current plan is to start with calculating the vertices of each triangle pair making a parallelogram or rhombus. With an Icosahedron there are 20 faces, making 10 parallelograms to sub-divide. We can subdivide one primary parallelogram (a quad) and translate the vertices of the remaining 9 faces without too much trouble.
- We can actually speed up calculation by mirroring each vertex about the latitudinal midline of the quad.
  - The vertex array size could then be pre-defined as we could very likely be able to calculate how each vertex gets addressed in the array to maintain consistency during face definition.
    - We will need a vertex array and a face array.
    - We could mirror each of the 10 quads simultaneously via translation.
    - <span style="color: orange;">*could I translate an entire line of vertices at a time?*</span>
      - <span style="color: orange;">\- acos(&radic;(5)/5)/X Lat &phi; from previous row.</span>
      - <span style="color: orange;">\+ &pi;/5 Long &lambda; from previous row. </span>(2&pi;/5)/2 = 2&pi;/10 = &pi;/5
        - <span style="color: red;">Need to verify that these equations work universally to translate each point.</span>
- With the inital quad, we determine the level of subdivision, dividing the primary edge by whatever number we desire. The benefit here is that we can make this number arbitrary. The initial version, this number was how many times we divided a primary edge in half. With this version we can divide the primary edge by any number that doesn't necessarily need to be half.
  - acos(&radic;(5)/5)/X; where X = what we're dividing the primary edge by.
  - This value represents the central angle for each edge along the primary edge as Latitude &phi;. Each row of vertices will be a factor of this number in angle.
- We will use the value for each row as the Latitude &phi; for the start point &phi;<sub>1</sub> and end point &phi;<sub>2</sub>.
- The Longitude &lambda; is calculated as 0 rads for &lambda;<sub>1</sub> and (2 \* &pi;)/5 for &lambda;<sub>2</sub>. This will be true until we reach the bottom of the first primary triangle.

- Ignore all of that BS above.
- We divide the edges of the first quad by the tesselation level. We already have the first four points for those quads so they can be plugged in directly.
  - Before, the tesselation level was how many times an edge was divided in half but now its a direct number to divide by.
- The vertices of the first 4 edges are saved to the vertex array in a clockwise order, starting from the North Pole.
  - This is so we know which vertices serve as outside vertices that get shared with the other 9 quads of the icosahedron.
  - When we know the array addresses of each subdivided vertex, we could determine face associated without even calculating the vertices. All calculating the vertices does is give us the coordinates of each vertex.

## <span style="background-color: blue;">Results and Discussion</span>

- Summarize the improvements.
- Discuss any remaining challenges.

## <span style="background-color: blue;">Conclusion</span>

- Recap the main points.
- Suggest next steps or further improvements.
