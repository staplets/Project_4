/***************************************
Group: 26
Members: Brett Fedack, Nikolay Goncharenko, Shaun Stapleton
CS325 Winter 2015
Project 4
Description: This program approximates the Traveling Salesman Problem (TSP)
    using a greedy algorithm that builds the solution inwards with
    progressively smaller convex-hulls. Input is read from a file, and the
    results are output to another file.
******************************************************************************/

#include <ctime>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include "fileutils.hpp"
#include "geoutils.hpp"
#include "Vector2D.hpp"

std::list<int> tsp(std::vector<Vector2D> &vertices);

int main(int argc, char **argv)
{
    // Exactly one command-line argument must be provided in addition to the
    // program name.
    if (argc != 2) {
        std::cerr << "usage: tsp <Filename>" << std::endl;
        exit(1);
    }

    std::string inFilename;         // Input filename
    std::string outFilename;        // Output filename
    std::vector<Vector2D> vertices; // Indexed vertex coordinates
    std::list<int> solution;        // TSP solution path
    unsigned long long pathLength;  // Length of TSP solution path

    // Store the input file's name from the command line.
    inFilename = std::string(argv[1]);

    // Name the output file as the input file's name with ".tour" appended.
    outFilename = inFilename + ".tour";

    // Read the vertices from the input file.
    vertices = readInput(inFilename);

    // Approximate the TSP solution path.
    solution = tsp(vertices);

    // Measure the TSP solution path.
    pathLength = measurePath(solution, vertices, true);

    // Record the TSP path.
    writeResults(solution, pathLength, outFilename);

    exit(0);
}


/******************************************************************************
Function: tsp

Purpose:  This function approximates a solution to the traveling salesman
          problem by iteratively processing convex-hulls of successively
          smaller sets of vertices. This function runs through the
          convex hull clockwise and counterclockwise to find the optimum
          choice.

Receives: vertices - indexed input sequence

Returns:  The non-optimal solution path is returned.

Pre:      At least three vertices shall be provided.

Post:     None
******************************************************************************/
std::list<int> tsp(std::vector<Vector2D> &vertices)
{
    std::list<int> solution;   // Non-optimal solution path
    std::list<int> clockSolution;   // Clockwise solution path
    std::list<int> counterSolution;   // Counter-clockwise solution path
    std::list<int> hull;       // Convex-hull
    std::vector<int> indices;  // Sequence of vertex indices
    std::vector<int> interior; // Indices not yet part of the solution path
    clock_t start, stop;       // Time measurements for profiling
    int clockwiseLength;      //length of path with clockwise vertices
    int counterClockLength;   //length of path with counterclockwise vertices

    start = clock();

    // Initialize the solution path as a convex-hull enclosing the vertices.
    for (int i = 0; i < vertices.size(); ++i) indices.push_back(i);
    solution = convexHull2D(indices, interior, vertices);

    //clockwise vs counterclockwise testing - verified to work
    clockSolution = solution;
    counterSolution = solution;
    

    // Route the path through vertices of progressively smaller convex-hulls.
    while (interior.size() >= 3) {

        // Calculate a convex-hull with the remaining vertices.
        indices.swap(interior);
        hull = convexHull2D(indices, interior, vertices);

        //sort the vertices in the hull before insertion
        auto middle = std::next(hull.begin(), hull.size() / 2);

        //create two lists to sort
        std::list<int> left(hull.begin(), middle), right(middle, hull.end());

        //sort left
        left.sort();

        //sort right
        right.sort();
        /*//Testing 
        for (auto iter = right.begin(); iter != right.end(); ++iter) {
             std::cout << "right: " << *iter << std::endl;
        }*/

        //merge the two lists back into the hull
        hull = left;
        hull.merge(right);
        
       
        //SS-- Running through the convex hull both clockwise
        //and counter-clockwise to take the optimal length.
        // Insert each vertex of the convex-hull into the counter-clockwise solution path.
        for (std::list<int>::reverse_iterator riter = hull.rbegin(); riter != hull.rend(); ++riter) {
             checkPath(*riter, counterSolution, vertices, true);
        }
        
        //Get counter-clockwise Length
        counterClockLength = measurePath(counterSolution, vertices, true);

        // Insert each vertex of the convex-hull into the clocksolution path.
        for (auto iter = hull.begin(); iter != hull.end(); ++iter) {
             insertIntoPath(*iter, clockSolution, vertices, true);
        }
        //Get counter-clockwise Length
        clockwiseLength = measurePath(clockSolution, vertices, true);
                
        if (counterClockLength <= clockwiseLength)
        { 
             // Insert if counter clockwise is optimal for iteration
             for (std::list<int>::reverse_iterator riter = hull.rbegin(); riter != hull.rend(); ++riter) {
                  insertIntoPath(*riter, solution, vertices, true);
             }
        }
        else{
             // Insert if clockwise is optimal for iteration
             for (auto iter = hull.begin(); iter != hull.end(); ++iter) {
                  insertIntoPath(*iter, solution, vertices, true);
             }
        }


        //Copy over the best path to the other solutions for the next iteration.
        clockSolution = solution;
        counterSolution = solution;
      
    }

    // Insert the one or two remaining vertices into the solution path.
    if (interior.size() > 0) {
        for (auto iter = interior.begin(); iter != interior.end(); ++iter) {
            insertIntoPath(*iter, solution, vertices, true);
        }
    }

    stop = clock();

    // Report the running-time.
    std::cout << "TSP path calculated in "
              << (stop - start) / static_cast<double> (CLOCKS_PER_SEC)
              << " seconds.\n";

    return solution;
}
