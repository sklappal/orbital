# orbital
Reaktor orbital challenge at https://reaktor.com/orbital-challenge/

Solved with Python. The logic is the following

- Read the data file 
- Transform positions from spherical to cartesian coordinates (origin at center of earth)
- Create a graph of the locations (including start and stop positions) based on visibility of satellites
- Two positions see each other if the line between does not intersect the globe
- Find the shortest path from start to stop with Djikstras algorithm
- Print the path


