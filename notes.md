####12.1 3D coordinate System

 - 6 rectangular coordinates
 - x,y,z.  positive x comes towards you, positive y to the right, positive z up
  - this is right handed
  - left handed is positive x to the right, positive y toward you and positive z up
 - 3 planes, xy (z=0), xz (y=0), yz (x=0), with origin at (0,0,0)
 - 8 octants
  - the first is {x > 0, y > 0, z > 0}
 - equation of a circle (centered on the origin) in xy-plane is `x^2 + y^2 = r^2`, with z=0
  - equation in xyz for a unit circle at z=3 would be `x^2 + y^2 = 1, z=3` 
 - distance between P1(X1, Y1, Z1) and P2(X2, Y2, Z2)
  - `sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1^2))`
 - similarly, a sphere has the equation `(x-x0)^2 + (y-y0)^2 + (z-z0^2) = r^2`
  - so centered at (1, -1, 2) with r = 3 is `(x-1)^2 + (y+1)^2 + (z-2)^2 = 3`

####12.2 Vectors

 - vectors have an initial and terminal point, *AB* has initial A and terminal B
 - the component form of a vector *AB* is A - B
 - if two vectos have the same component form, they are congruent AB = CD
 - these component forms are (x, y, z)
 - length of a vector, |V| = (sqrt(x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2) where _2 is on B and _1 is on A
 - vectors have magnitude and direction
    - 0 vector has no specific direction
 - adding vectors is done by simply adding their components
 - to get a unit vector you divide a vector by its magnitude, so unit(V) = vec(V) / |V|
 - subtraction is ordered and done componentwise A - B
 - midpoint of a line segment PQ is P+Q/2 componentwise

####12.3 Dot Product

 - Called inner or scalar product
 - U <dot> V = U1 * V1 + U2 * V2 + ... + UN * VN
 - U <dot> V = |U| * |V| * cos(theta)
    - so to find the angle between two vectors, theta = arccos((U <dot> V) / (|U|*|V|))
 - U <dot> V is 0 iff they are perpendicular
 - Projection
    - idea is to take a vector U and find its component parallel to vector V
    - P_v(U) = |U|cos(theta)
    - and the direction is Dir(P)
    - thus P_v(U) = |U||V|cos(theta)/|V|**2 * V

####12.4 lines and planes
 - 
