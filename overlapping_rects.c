#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// A structure to represent a 2D point
typedef struct {
    float x, y;
} Point;

// Function to calculate the dot product of two vectors
float dotProduct(Point v1, Point v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

// Function to subtract two points to get a vector
Point subtract(Point p1, Point p2) {
    return (Point){p1.x - p2.x, p1.y - p2.y};
}

// Function to calculate the cross product of two vectors
float crossProduct(Point v1, Point v2) {
    return v1.x * v2.y - v1.y * v2.x;
}

// Function to get the perpendicular vector (normal to an edge)
Point perpendicular(Point v) {
    return (Point){-v.y, v.x};
}

// Function to project a rectangle onto an axis and find the min/max projection
void projectRectangle(Point *rect, Point axis, float *min, float *max) {
    *min = dotProduct(rect[0], axis);
    *max = *min;
    for (int i = 1; i < 4; i++) {
        float projection = dotProduct(rect[i], axis);
        if (projection < *min) *min = projection;
        if (projection > *max) *max = projection;
    }
}

// Function to check if two projection ranges overlap
bool overlap(float minA, float maxA, float minB, float maxB) {
    return !(maxA < minB || maxB < minA);
}

// Function to check for overlap using the Separating Axis Theorem
bool checkOverlap(Point *rect1, Point *rect2) {
    // Combine edges from both rectangles
    Point edges[4] = {
        subtract(rect1[1], rect1[0]),
        subtract(rect1[3], rect1[0]),
        subtract(rect2[1], rect2[0]),
        subtract(rect2[3], rect2[0])
    };

    // Check each axis
    for (int i = 0; i < 4; i++) {
        // Get the axis perpendicular to the edge
        Point axis = perpendicular(edges[i]);

        // Project both rectangles onto the axis
        float min1, max1, min2, max2;
        projectRectangle(rect1, axis, &min1, &max1);
        projectRectangle(rect2, axis, &min2, &max2);

        // Check if the projections overlap
        if (!overlap(min1, max1, min2, max2)) {
            return false; // Found a separating axis
        }
    }

    return true; // No separating axis found, rectangles overlap
}

// Function to check if a point is inside a rectangle using cross products
bool pointInRectangle(Point point, Point *rect) {
    for (int i = 0; i < 4; i++) {
        Point edge = subtract(rect[(i + 1) % 4], rect[i]);
        Point toPoint = subtract(point, rect[i]);
        if (crossProduct(edge, toPoint) < 0) {
            return false; // Point is outside this edge
        }
    }
    return true; // Point is inside all edges
}

// Function to check if one rectangle is fully contained within another
bool checkContainment(Point *inner, Point *outer) {
    for (int i = 0; i < 4; i++) {
        if (!pointInRectangle(inner[i], outer)) {
            return false; // At least one point of the inner rectangle is outside
        }
    }
    return true; // All points of the inner rectangle are inside
}

// Main function
int main() {
    // Define a large rectangle (container)
    Point largeRect[4] = {
        {0, 0}, {10, 0}, {10, 8}, {0, 8}
    };

    // Define smaller rectangles (placed inside the large rectangle)
    Point smallRect1[4] = {
        {1, 1}, {3, 1}, {3, 2}, {1, 2}
    };
    Point smallRect2[4] = {
        {4, 1}, {6, 1}, {6, 3}, {4, 3}
    };
    Point smallRect3[4] = {
        {7, 2}, {9, 2}, {9, 4}, {7, 4}
    };

    // Check containment of smaller rectangles in the large rectangle
    Point *smallRects[] = {smallRect1, smallRect2, smallRect3};
    int numSmallRects = 3;

    bool allContained = true;
    for (int i = 0; i < numSmallRects; i++) {
        if (!checkContainment(smallRects[i], largeRect)) {
            allContained = false;
            printf("Small rectangle %d is not fully contained in the large rectangle.\n", i + 1);
        }
    }
    if (allContained) {
        printf("All small rectangles are fully contained within the large rectangle.\n");
    }

    // Check for overlap between smaller rectangles
    bool noOverlap = true;
    for (int i = 0; i < numSmallRects; i++) {
        for (int j = i + 1; j < numSmallRects; j++) {
            if (checkOverlap(smallRects[i], smallRects[j])) {
                noOverlap = false;
                printf("Small rectangle %d overlaps with small rectangle %d.\n", i + 1, j + 1);
            }
        }
    }
    if (noOverlap) {
        printf("No smaller rectangles overlap with each other.\n");
    }

    return 0;
}
