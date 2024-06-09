/**
 * Migros Delivery using Ant Colony Optimization
 * @author Can Emir Bora
 * @since Date: 11.05.2024
 */

import java.util.ArrayList;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

public class CanEmirBora {

    public static double shortestPathDistance = Double.MAX_VALUE; // data fields that I will use multiple times for saving the shortest path and distance
    public static Integer[] shortestPath = null;
    public static ArrayList<ArrayList<Double>> locations = new ArrayList<>();

    public static void main(String[] args) throws FileNotFoundException {

        int chosenMethod = 1; // 1-Brute-force method , 2-Ant Colony Optimization
        int chosenPlotting = 1; // For ant colony optimization: 1-Pheromone Level Map , 2- The Shortest Route
        String inputFile = "input03.txt";
        // hyper parameters
        int iterationCount = 300;
        int antCountPerIteration = 300;
        double degradationFactor = 0.9;
        double alpha = 1.2;
        double beta = 1.5;
        double initialPheromoneIntensity = 0.1;
        double qValue = 0.0001;

        File file = new File(inputFile); // opening and saving all the coordinates of houses in locations ArrayList
        Scanner scanner = new Scanner(file);
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine();
            String[] parts = line.split(",");
            ArrayList<Double> coordinates = new ArrayList<>();
            coordinates.add(Double.parseDouble(parts[0]));
            coordinates.add(Double.parseDouble(parts[1]));
            locations.add(coordinates);
        }
        scanner.close();

        if (chosenMethod == 1) {
            bruteForce();
        } else if (chosenMethod == 2) {
            antColonyOptimization(iterationCount, antCountPerIteration, degradationFactor, alpha, beta, initialPheromoneIntensity, qValue, chosenPlotting);
        }
    }

     /**
     * Brute-Force Method:Generate all possible permutations and try all of them. Pick the one which produces the shortest path.
     */

    private static void bruteForce() {
        long startTime = System.currentTimeMillis();
        int migrosIdx = 0; // Migros is always the first home in the inputs
        Integer[] initialRoute = new Integer[locations.size() - 1]; // creating a random sorted initial route in order to make permutations on it
        int index = 0;
        for (int i = 0; i < locations.size(); i++) {
            if (i != migrosIdx) {
                initialRoute[index++] = i;
            }
        }
        permute(initialRoute, 0); // making all the permutations and finding the shortest route
        long endTime = System.currentTimeMillis();
        System.out.println("Method: Brute-Force Method"); // console outputs
        String totalDistanceStr = String.format("Total Distance: %.5f ", shortestPathDistance);
        System.out.println("Shortest Distance: " + totalDistanceStr);
        ArrayList<Integer> outputShortestPath = new ArrayList<>(); // creating output ArrayList in order to add Migros to the output route
        outputShortestPath.add(1);
        for (int num : shortestPath) {
            outputShortestPath.add(num + 1);
        }
        outputShortestPath.add(1);
        System.out.println("Shortest Path: " + Arrays.toString(outputShortestPath.toArray()));
        long duration = (endTime - startTime);
        System.out.println("Time it takes to find the shortest path: " + duration / 1000.0 + " seconds."); // printing how long it takes to use this method
        //StdDraw
        StdDraw.enableDoubleBuffering();
        int canvasWidth = 800;
        int canvasHeight = 800;
        StdDraw.setCanvasSize(canvasWidth, canvasHeight);
        StdDraw.setXscale(0, 1);
        StdDraw.setYscale(0, 1);
        double circleRadius = 0.02;
        StdDraw.setPenColor(StdDraw.BLACK);
        StdDraw.setPenRadius(0.002);
        ArrayList<Double> firstLocation = locations.get(migrosIdx);
        double locationX = firstLocation.get(0);
        double locationY = firstLocation.get(1);
        for (int idx : shortestPath) { // drawing the shortest path lines
            ArrayList<Double> nextLocation = locations.get(idx);
            double nextLocationX = nextLocation.get(0);
            double nextLocationY = nextLocation.get(1);
            StdDraw.line(locationX, locationY, nextLocationX, nextLocationY);
            locationX = nextLocationX;
            locationY = nextLocationY;
        }
        StdDraw.line(locationX, locationY, firstLocation.get(0), firstLocation.get(1));
        for (int i = 0; i < locations.size(); i++) { //drawing all homes as a filled circle.
            ArrayList<Double> location = locations.get(i);
            if (i == migrosIdx) {
                StdDraw.setPenColor(StdDraw.PRINCETON_ORANGE);// painting Migros orange.
            } else {
                StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
            }
            StdDraw.filledCircle(location.get(0), location.get(1), circleRadius);
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.text(location.get(0), location.get(1), Integer.toString(i + 1));
        }
        StdDraw.show();
    }

    /**
     * It finds all possible combinations of given initial route and updating the shortest path accordingly.
     * @param arr is initial sorted route that I make all permutations on.
     * @param k is an auxiliary number that tells me to stop recursion when it reaches the length of the initial route.
     */

    private static void permute(Integer[] arr, int k) {
        if (k == arr.length) { // stopping condition
            double distance = calculateRouteDistance(arr.clone()); // calculating the distance of that permutation
            if (distance < shortestPathDistance) { // checking and updating the shortest path and shortest distance values
                shortestPathDistance = distance;
                shortestPath = arr.clone();
            }
        } else {
            for (int i = k; i < arr.length; i++) { // recursive call for permutation
                Integer temp = arr[i];
                arr[i] = arr[k];
                arr[k] = temp;
                permute(arr, k + 1);
                temp = arr[k];
                arr[k] = arr[i];
                arr[i] = temp;
            }
        }
    }

    /**
     * It calculates the distance of given route.
     * @param route is the route that we want to know the length of
     * @return distance of that route
     */

    private static double calculateRouteDistance(Integer[] route) {
        double distance = 0;
        int prevIdx = 0;
        for (int locationIdx : route) { // summing all the roads between each home in the route
            distance += Math.sqrt(Math.pow(locations.get(locationIdx).get(0) - locations.get(prevIdx).get(0), 2) + Math.pow(locations.get(locationIdx).get(1) - locations.get(prevIdx).get(1), 2));
            prevIdx = locationIdx;
        }
        distance += Math.sqrt(Math.pow(locations.getFirst().getFirst() - locations.get(prevIdx).get(0), 2) + Math.pow(locations.getFirst().get(1) - locations.get(prevIdx).get(1), 2));
        return distance;
    }

    /**
     *
     * @param iterationCount is the number of how many rounds of ants I will send.
     * @param antCountPerIteration is the number of ants that I will send in each iteration.
     * @param degradationFactor is multiplier that I will multiply all the pheromone values at the end of each iteration.
     * @param alpha is the value that is the exponent of the pheromone value in the probability formula used when choosing the next node that the ant will go to.
     * @param beta the value that is the exponent of the distance value in the probability formula used when choosing the next node that the ant will go to.
     * @param initialPheromoneIntensity is the initial pheromone level at every road.
     * @param qValue is a factor of pheromone value update process after an ant completes a route.
     * @param chosenPlotting is the value that indicating whether to draw a pheromone level map or the shortest path.
     */

    private static void antColonyOptimization(int iterationCount, int antCountPerIteration, double degradationFactor, double alpha, double beta, double initialPheromoneIntensity, double qValue, int chosenPlotting) {
        long startTime = System.currentTimeMillis();
        int migrosIdx = 0;
        double[][] pheromoneMatrix = new double[locations.size()][locations.size()]; // creating a 2-D array in order to save all pheromone levels of ROAD between 2 homes.
        for (int i = 0; i < locations.size(); i++) {
            for (int j = 0; j < locations.size(); j++) {
                if (i != j) {
                    pheromoneMatrix[i][j] = initialPheromoneIntensity; // initially setting pheromone levels as 0.1 .
                }
            }
        }

        for (int iteration = 0; iteration < iterationCount; iteration++) { // number of iteration times
            for (int ants = 0; ants < antCountPerIteration; ants++) { // number of ant count per iteration times
                Integer[] currentRoute = new Integer[locations.size() - 1]; // saving current route in order to use it later
                ArrayList<Integer> visitedHomesIdx = new ArrayList<>(); // saving visited homes in order to not go again
                int currentIndex = migrosIdx;
                visitedHomesIdx.add(currentIndex);
                for (int i = 0; i < currentRoute.length; i++) { // until coming back to Migros
                    int nextIndex = nextNodeRandomly(currentIndex, visitedHomesIdx, pheromoneMatrix, alpha, beta); // probabilistically choosing the next house
                    currentRoute[i] = nextIndex; // updating values
                    visitedHomesIdx.add(nextIndex);
                    currentIndex = nextIndex;
                }

                double distance = calculateRouteDistance(currentRoute); // calculate the length of the shortest path recorded
                if (distance < shortestPathDistance) { // update the shortest path
                    shortestPathDistance = distance;
                    shortestPath = currentRoute.clone();
                }

                for (int i = 0; i < visitedHomesIdx.size(); i++) { // update the pheromone level of the path traveled
                    int fIdx = visitedHomesIdx.get(i);
                    int sIdx;
                    if (i + 1 == visitedHomesIdx.size()) {
                        sIdx = migrosIdx;
                    } else {
                        sIdx = visitedHomesIdx.get(i + 1);
                    }
                    pheromoneMatrix[fIdx][sIdx] += qValue / distance; //updating
                    pheromoneMatrix[sIdx][fIdx] += qValue / distance;
                }
            }

            for (int i = 0; i < locations.size(); i++) { // after each iteration, lower the pheromone level on all paths
                for (int j = 0; j < locations.size(); j++) {
                    pheromoneMatrix[i][j] = pheromoneMatrix[i][j] * degradationFactor;
                }
            }

        }

        long endTime = System.currentTimeMillis();
        System.out.println("Method: Ant Colony Optimization Method"); // console outputs
        String totalDistanceStr = String.format("Total Distance: %.5f ", shortestPathDistance);
        System.out.println("Shortest Distance: " + totalDistanceStr);
        ArrayList<Integer> outputShortestPath = new ArrayList<>();
        outputShortestPath.add(1);
        for (int num : shortestPath) {
            outputShortestPath.add(num + 1);
        }
        outputShortestPath.add(1);
        System.out.println("Shortest Path: " + Arrays.toString(outputShortestPath.toArray()));
        long duration = (endTime - startTime);
        System.out.println("Time it takes to find the shortest path: " + duration / 1000.0 + " seconds.");// printing how long it takes to use this method
        //StdDraw
        if (chosenPlotting == 1) {// if I want to draw pheromone level map

            StdDraw.enableDoubleBuffering();
            int canvasWidth = 800;
            int canvasHeight = 800;
            StdDraw.setCanvasSize(canvasWidth, canvasHeight);
            StdDraw.setXscale(0, 1);
            StdDraw.setYscale(0, 1);
            StdDraw.setPenColor(StdDraw.BLACK);
            for (int i = 0; i < pheromoneMatrix.length; i++) { // drawing pheromone levels
                for (int j = 0; j < pheromoneMatrix[i].length; j++) {
                    StdDraw.setPenRadius(pheromoneMatrix[i][j] * 0.1); //changing thickness according to pheromone level
                    StdDraw.line(locations.get(i).get(0), locations.get(i).get(1), locations.get(j).get(0), locations.get(j).get(1));
                }
            }
            StdDraw.setPenRadius(0.002);
            double circleRadius = 0.02;
            for (int i = 0; i < locations.size(); i++) { // drawing homes as filled circle
                ArrayList<Double> location = locations.get(i);
                StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
                StdDraw.filledCircle(location.get(0), location.get(1), circleRadius);
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(location.get(0), location.get(1), Integer.toString(i + 1));
            }
            StdDraw.show();

        } else if (chosenPlotting == 2) { // if I want to draw the shortest path

            StdDraw.enableDoubleBuffering();
            int canvasWidth = 800;
            int canvasHeight = 800;
            StdDraw.setCanvasSize(canvasWidth, canvasHeight);
            StdDraw.setXscale(0, 1);
            StdDraw.setYscale(0, 1);
            double circleRadius = 0.02;
            StdDraw.setPenColor(StdDraw.BLACK);
            StdDraw.setPenRadius(0.002);
            ArrayList<Double> firstLocation = locations.get(migrosIdx);
            double locationX = firstLocation.get(0);
            double locationY = firstLocation.get(1);
            for (int idx : shortestPath) { // drawing the shortest path lines
                ArrayList<Double> nextLocation = locations.get(idx);
                double nextLocationX = nextLocation.get(0);
                double nextLocationY = nextLocation.get(1);
                StdDraw.line(locationX, locationY, nextLocationX, nextLocationY);
                locationX = nextLocationX;
                locationY = nextLocationY;
            }
            StdDraw.line(locationX, locationY, firstLocation.get(0), firstLocation.get(1));
            for (int i = 0; i < locations.size(); i++) { // drawing homes as filled circle
                ArrayList<Double> location = locations.get(i);
                if (i == migrosIdx) {
                    StdDraw.setPenColor(StdDraw.PRINCETON_ORANGE);
                } else {
                    StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
                }
                StdDraw.filledCircle(location.get(0), location.get(1), circleRadius);
                StdDraw.setPenColor(StdDraw.BLACK);
                StdDraw.text(location.get(0), location.get(1), Integer.toString(i + 1));
            }
            StdDraw.show();
        }
    }

    /**
     *It is used to probabilistically go to other unvisited nodes from the given index.
     * @param currentIndex keeps the index of the current node
     * @param visitedHomesIdx prevents returning to visited nodes
     * @param pheromoneMatrix keeps a 2-D matrix showing the pheromone levels of the roads. It is used in the formula while calculating probability.
     * @param alpha is the value that is the exponent of the pheromone value in the probability formula.
     * @param beta the value that is the exponent of the distance value in the probability formula.
     * @return It returns probabilistically chosen next home's index
     */

    private static int nextNodeRandomly(int currentIndex, ArrayList<Integer> visitedHomesIdx, double[][] pheromoneMatrix, double alpha, double beta) {
        double total = 0; // total probability
        for (int i = 0; i < locations.size(); i++) { // calculating probability by formula
            if (!visitedHomesIdx.contains(i)) { // checking whether it is visited
                total += Math.pow(pheromoneMatrix[currentIndex][i], alpha) / Math.pow(calculateDistance(locations.get(currentIndex).get(0), locations.get(currentIndex).get(1), locations.get(i).get(0), locations.get(i).get(1)), beta);
            }
        }
        double randomValue = Math.random() * total; // picking random number between 0 - 1.0
        double probability = 0;
        for (int i = 0; i < locations.size(); i++) { // adding all probabilities of roads until it reaches randomValue.
            if (!visitedHomesIdx.contains(i)) {
                probability += Math.pow(pheromoneMatrix[currentIndex][i], alpha) / Math.pow(calculateDistance(locations.get(currentIndex).get(0), locations.get(currentIndex).get(1), locations.get(i).get(0), locations.get(i).get(1)), beta);
                if (probability >= randomValue) { // if it reaches, choose that road
                    return i;
                }
            }
        }
        return -1;
    }

    /**
     * Calculating the distance between two given points coordinates
     * @param fX is the first points X coordinate
     * @param fY is the first points Y coordinate
     * @param lX is the second points X coordinate
     * @param lY is the second points Y coordinate
     * @return It returns the distance between these two points
     */

    private static double calculateDistance(double fX, double fY, double lX, double lY) {
        return Math.pow(Math.pow(fX - lX, 2) + Math.pow(fY - lY, 2), 0.5);
    }

}