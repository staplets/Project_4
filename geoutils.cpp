#include "geoutils.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <list>
#include <utility>
#include <vector>
#include "PolarComparator.hpp"
#include "Vector2D.hpp"

/******************************************************************************
Function: convexHull2D

Purpose:  This function calculates the convex hull of the given points using
          the Graham Scan algorithm:

            http://en.wikipedia.org/wiki/Graham_scan

Receives: indices - points from which to calculate the convex hull

          interior - indices of points inside of the convex hull

          pts - indexed point coordinates

Returns:  The sequence of point indices along the hull is returned.

Pre:      All point indices shall exist within the bounds of the sequence of
          indexed coordinates.

          No two points shall have the same coordinates.

          At least three points shall be provided.

Post:     "indices" and "interior" are modified.
******************************************************************************/
std::list<int> convexHull2D(std::vector<int> &indices,
                            std::vector<int> &interior,
                            std::vector<Vector2D> &pts )
{
    int tail, corner, head; // Adjacent indices along the convex-hull
    PolarComparator comparator; // Comparator for radial sorting
    std::list<int> hull; // Return value

    // Move the index of the lower-left point to the front of the sequence of
    // points from which to calculate the convex-hull.
    for (int i = 1; i < indices.size(); ++i) {

        // Compare the current point against the one at the front.
        if (pts[indices[i]].y < pts[indices[0]].y        // bottom-most
            || (pts[indices[i]].y == pts[indices[0]].y   // left-most
                && pts[indices[i]].x < pts[indices[0]].x)
        ) {
            std::swap(indices[0], indices[i]);
        }
    }

    // Radially sort the points relative to first point.
    comparator = PolarComparator(indices[0], &pts);
    std::sort(indices.begin() + 1, indices.end(), comparator);

    // Compute the convex-hull.
    indices.push_back(indices[0]); // Wrap the sequence of indices.
    interior.clear();              // No interior points have been found yet.
    hull.push_back(indices[0]);    // The first to points are on the hull.
    hull.push_back(indices[1]);
    tail = indices[1];
    corner = indices[2];
    for (int i = 3; i < indices.size(); ++i) {

        // Advance the head.
        head = indices[i];

        // Concave corners are not part of the convex-hull.
        while (((pts[head] - pts[corner]).cross(pts[tail] - pts[corner]) < 0)) {
            interior.push_back(corner);
            corner = hull.back();
            hull.pop_back();
            tail = hull.back();
        }

        // Convex and collinear corners are part of the convex-hull.
        hull.push_back(corner);

        // Advance the tail and corner.
        tail = corner;
        corner = head;
    }

    return hull;
}


/******************************************************************************
Function: closestPointOnSegment

Purpose:  This function determines the closest point on the line segment from
          the given point.

Receives: pt - point from which to find closest point on the line segment

          endPt1 - first endpoint of the line segment

          endPt2 - second endpoint of the line segment

Returns:  The closest point on the line segment is returned as a Vector2D
          object.

Pre:      The two endpoints shall not have the same coordinates.

Post:     None
******************************************************************************/
Vector2D closestPointOnSegment(Vector2D &pt, Vector2D &endPt1, Vector2D &endPt2)
{
    double m; // Slope of the input line segment
    double dx, dy; // Slope components
    double threshold; // Arbitrary permitted error in slope components
    Vector2D intxnPt; // Closest point on the line from the given point

    // Calculate the slope components.
    dx = endPt2.x - endPt1.x;
    dy = endPt2.y - endPt1.y;

    // Vertical and horizontal lines require special handling to avoid
    // divide-by-zero errors.
    threshold = 0.0001;
    if (abs(dx) < threshold) {
        intxnPt.x = endPt1.x;
        intxnPt.y = pt.y;
    } else if (abs(dy) < threshold) {
        intxnPt.x = pt.x;
        intxnPt.y = endPt1.y;
    }

    // The line segment is neither vertical nor horizontal.
    else {
        // Calculate the slope of the input line segment.
        m = dy / dx;

        // Finding the closest point on a line involves determining where two
        // perpendicular lines intersect.

        // The line equations are:
        // m = (intxnPt.y - endPt1.y) / (intxnPt.x - endPt1.x)
        // -1 / m = (intxnPt.y - pt.y) / (intxnPt.x - pt.x)

        // Running the equations through Mathematica results in the following:
        intxnPt.x = (pt.x + m * (pt.y - endPt1.y + endPt1.x * m)) / (1 + m * m);
        intxnPt.y = m * (intxnPt.x - endPt1.x) + endPt1.y;
    }

    // This intersection is the closest point only if it exists within the
    // bounds of the line segment.
    if ((intxnPt - endPt1).lengthSquared() < (endPt2 - endPt1).lengthSquared()
        && (intxnPt - endPt2).lengthSquared() < (endPt1 - endPt2).lengthSquared()
    ) {
        return intxnPt;
    }

    // Otherwise, return the closer of the two endpoints.
    else if ((endPt1 - pt).lengthSquared() < (endPt2 - pt).lengthSquared()) {
        return endPt1;
    } else {
        return endPt2;
    }
}


/******************************************************************************
Function: insertIntoPath

Purpose:  This function inserts a given vertex into the specified path by
          rerouting the nearest edge.

Receives: vertIdx - index of vertex to insert into the path

          path - sequence of connected vertices

          vertices - indexed vertex coordinates

          pathIsCircuit - flag indicating an edge between the path endpoints

Returns:  None

Pre:      All point indices shall exist within the bounds of the sequence of
          indexed coordinates.

          At least two vertices shall exist in the path.

Post:     "path" is modified.
******************************************************************************/
void insertIntoPath(int vertIdx,
                    std::list<int> &path,
                    std::vector<Vector2D> &vertices,
                    bool pathIsCircuit              )
{
    Vector2D closestPt; // Closest point on an path from given vertex
    double distSquared, minDistSquared; // Squares of distances from vertex
    std::list<int>::iterator insert, prev, curr; // Positions in path

    // Traverse the path.
    minDistSquared = DBL_MAX;
    prev = path.begin();
    curr = ++path.begin();
    while (prev != path.end()) {

        // Account for the edge that closes the path, if necessary.
        if (curr == path.end()) {
            if (pathIsCircuit) {
                curr = path.begin();
            } else {
                break;
            }
        }

        // Determine the shortest distance squared to the current edge.
        closestPt =  closestPointOnSegment(
            vertices[vertIdx], vertices[*prev], vertices[*curr]
        );
        distSquared = (closestPt - vertices[vertIdx]).lengthSquared();

        // Update the minimum distance and insertion point accordingly.
        if (distSquared < minDistSquared) {
            minDistSquared = distSquared;
            insert = curr;
        }

        // Advance to the next edge.
        prev++;
        curr++;
    }

    // Insert the vertex, effectively rerouting the closest edge.
    path.insert(insert, vertIdx);
}


/******************************************************************************
Function: measurePath

Purpose:  This function measures the given path's length as the sum of
          individual edge lengths where each edge length is rounded to the
          nearest integer.

Receives: path - sequence of connected vertices

          vertices - indexed vertex coordinates

          pathIsCircuit - flag indicating an edge between the path endpoints

Returns:  The path length is returned.

Pre:      All point indices shall exist within the bounds of the sequence of
          indexed coordinates.

          At least two vertices shall exist in the path.

Post:     None
******************************************************************************/
unsigned long long measurePath(std::list<int> &path,
                               std::vector<Vector2D> &vertices,
                               bool pathIsCircuit              )
{
    std::list<int>::iterator prev, curr; // Positions in path
    double edgeLength; // Rounded to nearest integer
    unsigned long long pathLength; // Return value

    // Traverse the path.
    pathLength = 0;
    prev = path.begin();
    curr = ++path.begin();
    while (prev != path.end()) {

        // Account for the edge that closes the path, if necessary.
        if (curr == path.end()) {
            if (pathIsCircuit) {
                curr = path.begin();
            } else {
                break;
            }
        }

        // Add the edge length to the path length.
        edgeLength = round(
            sqrt(
                (vertices[*curr] - vertices[*prev]).lengthSquared()
            )
        );
        pathLength += edgeLength;

        // Advance to the next edge.
        prev++;
        curr++;
    }

    return pathLength;
}

/******************************************************************************
Function: checkPath

Purpose:  This function inserts a given vertex into the specified path by
rerouting the nearest edge for our list comparing.

Receives: vertIdx - index of vertex to insert into the path

path - sequence of connected vertices

vertices - indexed vertex coordinates

pathIsCircuit - flag indicating an edge between the path endpoints

Returns:  None

Pre:      All point indices shall exist within the bounds of the sequence of
indexed coordinates.

At least two vertices shall exist in the path.

Post:     "path" is modified.
******************************************************************************/
void checkPath(int vertIdx,
     std::list<int> &testPath,
     std::vector<Vector2D> &vertices,
     bool pathIsCircuit)
{
     Vector2D closestPt; // Closest point on an path from given vertex
     double distSquared, minDistSquared; // Squares of distances from vertex
     std::list<int>::iterator insert, prev, curr; // Positions in path

     // Traverse the path.
     minDistSquared = DBL_MAX;
     prev = testPath.begin();
     curr = ++testPath.begin();
     while (prev != testPath.end()) {

          // Account for the edge that closes the path, if necessary.
          if (curr == testPath.end()) {
               if (pathIsCircuit) {
                    curr = testPath.begin();
               }
               else {
                    break;
               }
          }

          // Determine the shortest distance squared to the current edge.
          closestPt = closestPointOnSegment(
               vertices[vertIdx], vertices[*prev], vertices[*curr]
               );
          distSquared = (closestPt - vertices[vertIdx]).lengthSquared();

          // Update the minimum distance and insertion point accordingly.
          if (distSquared < minDistSquared) {
               minDistSquared = distSquared;
               insert = curr;
          }

          // Advance to the next edge.
          prev++;
          curr++;
     }

     // Insert the vertex, effectively rerouting the closest edge.
     testPath.insert(insert, vertIdx);
}
