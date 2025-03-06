#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define POPULATION_SIZE 100
#define MAX_GENERATIONS 1000
#define MUTATION_RATE 0.1

typedef struct {
    float x, y;     // Center of the rectangle
    float width;    // Width of the rectangle
    float height;   // Height of the rectangle
    float angle;    // Rotation angle (radians)
} RotatedRectangle;

typedef struct {
    float x, y;
} Point;

typedef struct {
    Point corners[4];
} Polygon;

// Random number generator in a range
float random_range(float min, float max) {
    return min + (float)rand() / (float)(RAND_MAX / (max - min));
}

// Generate a random rotated rectangle
RotatedRectangle random_rectangle(float max_width, float max_height) {
    RotatedRectangle rect;
    rect.x = random_range(0, max_width);
    rect.y = random_range(0, max_height);
    rect.width = random_range(0.1, max_width / 2);  // Avoid rectangles too small or too large
    rect.height = random_range(0.1, max_height / 2);
    rect.angle = random_range(0, 2 * M_PI);         // Angle in radians
    return rect;
}

// Calculate the area of a rectangle
float rectangle_area(RotatedRectangle rect) {
    return rect.width * rect.height;
}

// Check if a rectangle is valid (within bounds and outside obstacles)
bool rectangle_is_valid(RotatedRectangle rect, RotatedRectangle large_rect, RotatedRectangle small_rects[], int num_small_rects) {
    // TODO: Implement collision and boundary checks
    return true;
}

// Genetic algorithm function
RotatedRectangle genetic_algorithm(RotatedRectangle large_rect, RotatedRectangle small_rects[], int num_small_rects) {
    RotatedRectangle population[POPULATION_SIZE];
    float fitness[POPULATION_SIZE];
    RotatedRectangle best_solution;
    float best_fitness = -1;

    // Initialize population
    for (int i = 0; i < POPULATION_SIZE; i++) {
        population[i] = random_rectangle(large_rect.width, large_rect.height);
    }

    for (int generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Evaluate fitness
        for (int i = 0; i < POPULATION_SIZE; i++) {
            fitness[i] = rectangle_is_valid(population[i], large_rect, small_rects, num_small_rects)
                             ? rectangle_area(population[i])
                             : 0;

            if (fitness[i] > best_fitness) {
                best_solution = population[i];
                best_fitness = fitness[i];
            }
        }

        // Termination condition
        if (generation > 50 && fabs(best_fitness - fitness[0]) < 1e-3) {
            break;
        }

        // TODO: Selection, Crossover, Mutation
    }

    return best_solution;
}

