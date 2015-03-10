#include <list>
#include <vector>
#include "Vector2D.hpp"

#pragma once

// Calculates a convex-hull from the given points
std::list<int> convexHull2D(std::vector<int> &indices,
                            std::vector<int> &interior,
                            std::vector<Vector2D> &pts );

// Determines the closest point on a line segment from the given point
Vector2D closestPointOnSegment(Vector2D &pt, Vector2D &endPt1, Vector2D &endPt2);

// Inserts a vertex into the specified path
void insertIntoPath(int vertIdx,
                    std::list<int> &path,
                    std::vector<Vector2D> &vertices,
                    bool pathIsCircuit = false      );

// Measures a path's length, rounding each edge length to the nearest integer.
unsigned long long measurePath(std::list<int> &path,
                               std::vector<Vector2D> &vertices,
                               bool pathIsCircuit = false      );

//Determins the closest point on a line segment on a test list for both
//clockwise and counterclockwise ordering.
void checkPath(int vertIdx,
     std::list<int> &testPath,
     std::vector<Vector2D> &vertices,
     bool pathIsCircuit);