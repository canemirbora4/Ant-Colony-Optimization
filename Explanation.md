This is a Java program which finds the quickest delivery route of a Migros delivery car. In the map, there is a Migros store, denoted by the orange circle, and several houses, denoted by the gray circles. Each
house places an order from Migros. The delivery car should visit each house only once and return to Migros. Coordinates of Migros and houses are given in an input text file as shown
where the first line denotes the Migros coordinate, and the rest denote house coordinates. During the delivery, each house should be visited at most once. The program should find
the shortest route and plot the final route.

1. Brute-Force Approach: Generate all possible permutations and try all of them. Pick the one which produces the shortest
path. This method is very slow but finds the optimal solution.

3. Ant Colony Optimization Approach: This optimization algorithm may not find the best solution, but it is considerably
faster than the brute-force approach. Usually, ant colony optimization heuristic finds a very close solution and may be
preferred in some situations.

I also created and plotted the density map of the pheromones left by the ants on the roads.
