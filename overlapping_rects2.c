#include <stdio.h>
#include <stdbool.h>
#include <float.h>

typedef struct {
    float x, y;
} Point;

typedef struct {
    Point topLeft, bottomRight;
} Rectangle;

// Function to find the maximum of two numbers
float max(float a, float b) {
    return (a > b) ? a : b;
}

// Function to find the minimum of two numbers
float min(float a, float b) {
    return (a < b) ? a : b;
}

// Function to check if two rectangles overlap
bool checkOverlapSimple(Rectangle r1, Rectangle r2) {
    return !(r1.bottomRight.x <= r2.topLeft.x || r1.topLeft.x >= r2.bottomRight.x ||
             r1.bottomRight.y <= r2.topLeft.y || r1.topLeft.y >= r2.bottomRight.y);
}

// Function to check if one rectangle is contained within another
bool checkContainmentSimple(Rectangle inner, Rectangle outer) {
    return (inner.topLeft.x >= outer.topLeft.x && inner.bottomRight.x <= outer.bottomRight.x &&
            inner.topLeft.y >= outer.topLeft.y && inner.bottomRight.y <= outer.bottomRight.y);
}

// Function to calculate the area of a rectangle
float rectangleArea(Rectangle rect) {
    float width = rect.bottomRight.x - rect.topLeft.x;
    float height = rect.bottomRight.y - rect.topLeft.y;
    return (width > 0 && height > 0) ? width * height : 0;
}

// Function to calculate the largest rectangle that fits within a given space
Rectangle findLargestRectangle(Rectangle largeRect, Rectangle *smallRects, int numSmallRects) {
    Rectangle bestRectangle = {{0, 0}, {0, 0}};
    float bestArea = 0;

    // Test all potential subregions by extending the edges of the smaller rectangles
    for (int i = 0; i < numSmallRects; i++) {
        Rectangle current = smallRects[i];

        // Define potential candidate rectangles by using the edges of the current small rectangle
        Rectangle candidates[4] = {
            {{largeRect.topLeft.x, largeRect.topLeft.y}, {current.topLeft.x, largeRect.bottomRight.y}},  // Left region
            {{current.bottomRight.x, largeRect.topLeft.y}, {largeRect.bottomRight.x, largeRect.bottomRight.y}},  // Right region
            {{largeRect.topLeft.x, current.bottomRight.y}, {largeRect.bottomRight.x, largeRect.bottomRight.y}},  // Bottom region
            {{largeRect.topLeft.x, largeRect.topLeft.y}, {largeRect.bottomRight.x, current.topLeft.y}}   // Top region
        };

        // Check each candidate rectangle
        for (int j = 0; j < 4; j++) {
            Rectangle candidate = candidates[j];

            // Check if the candidate is inside the large rectangle
            if (!checkContainmentSimple(candidate, largeRect)) {
                continue;
            }

            // Check if the candidate overlaps with any smaller rectangle
            bool overlaps = false;
            for (int k = 0; k < numSmallRects; k++) {
                if (checkOverlapSimple(candidate, smallRects[k])) {
                    overlaps = true;
                    break;
                }
            }

            // If no overlaps, calculate the area and update the best rectangle
            if (!overlaps) {
                float area = rectangleArea(candidate);
                if (area > bestArea) {
                    bestArea = area;
                    bestRectangle = candidate;
                }
            }
        }
    }

    return bestRectangle;
}

int main() {
    // Define a large rectangle (container)
    Rectangle largeRect = {{0, 0}, {10, 8}};

    // Define smaller rectangles (placed inside the large rectangle)
    Rectangle smallRects[] = {
        {{1, 1}, {3, 2}},  // Small rectangle 1
        {{4, 1}, {6, 3}},  // Small rectangle 2
        {{7, 2}, {9, 4}}   // Small rectangle 3
    };
    int numSmallRects = 3;

    // Find the largest rectangle that fits
    Rectangle largest = findLargestRectangle(largeRect, smallRects, numSmallRects);

    // Output the result
    printf("The largest rectangle that fits is:\n");
    printf("Top Left: (%.2f, %.2f)\n", largest.topLeft.x, largest.topLeft.y);
    printf("Bottom Right: (%.2f, %.2f)\n", largest.bottomRight.x, largest.bottomRight.y);
    printf("Area: %.2f\n", rectangleArea(largest));

    return 0;
}

