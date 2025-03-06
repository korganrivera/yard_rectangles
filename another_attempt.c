#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <float.h>

#define POPULATION_SIZE 100
#define GENERATIONS 1000
#define MUTATION_RATE 0.1

typedef struct {
    float x1, y1, x2, y2;
    float fitness;
} Rectangle;

typedef struct {
    float x1, y1, x2, y2;
} Bounds;

// Check if two rectangles overlap
bool rectanglesOverlap(Rectangle r1, Rectangle r2) {
    return !(r1.x2 <= r2.x1 || r1.x1 >= r2.x2 || r1.y2 <= r2.y1 || r1.y1 >= r2.y2);
}

// Check if a rectangle is valid and does not overlap any obstacles
bool isValid(Rectangle rect, Rectangle *smallRects, int numSmallRects, Bounds bounds) {
    if (rect.x1 < bounds.x1 || rect.y1 < bounds.y1 || rect.x2 > bounds.x2 || rect.y2 > bounds.y2) {
        return false; // Out of bounds
    }
    for (int i = 0; i < numSmallRects; i++) {
        if (rectanglesOverlap(rect, smallRects[i])) {
            return false;
        }
    }
    return true;
}

// Calculate the fitness (area) of a rectangle
float calculateFitness(Rectangle rect, Rectangle *smallRects, int numSmallRects, Bounds bounds) {
    if (!isValid(rect, smallRects, numSmallRects, bounds)) {
        return 0; // Invalid rectangles have zero fitness
    }
    return (rect.x2 - rect.x1) * (rect.y2 - rect.y1);
}

// Generate a random rectangle within bounds
Rectangle generateRandomRectangle(Bounds bounds) {
    float x1 = bounds.x1 + (float)rand() / RAND_MAX * (bounds.x2 - bounds.x1);
    float y1 = bounds.y1 + (float)rand() / RAND_MAX * (bounds.y2 - bounds.y1);
    float x2 = x1 + (float)rand() / RAND_MAX * (bounds.x2 - x1);
    float y2 = y1 + (float)rand() / RAND_MAX * (bounds.y2 - y1);
    return (Rectangle){x1, y1, x2, y2, 0};
}

// Perform crossover between two rectangles
Rectangle crossover(Rectangle parent1, Rectangle parent2) {
    float x1 = (parent1.x1 + parent2.x1) / 2;
    float y1 = (parent1.y1 + parent2.y1) / 2;
    float x2 = (parent1.x2 + parent2.x2) / 2;
    float y2 = (parent1.y2 + parent2.y2) / 2;
    return (Rectangle){x1, y1, x2, y2, 0};
}

// Mutate a rectangle
Rectangle mutate(Rectangle rect, Bounds bounds) {
    if ((float)rand() / RAND_MAX < MUTATION_RATE) {
        rect.x1 += ((float)rand() / RAND_MAX - 0.5) * (bounds.x2 - bounds.x1) * 0.1;
        rect.y1 += ((float)rand() / RAND_MAX - 0.5) * (bounds.y2 - bounds.y1) * 0.1;
        rect.x2 += ((float)rand() / RAND_MAX - 0.5) * (bounds.x2 - bounds.x1) * 0.1;
        rect.y2 += ((float)rand() / RAND_MAX - 0.5) * (bounds.y2 - bounds.y1) * 0.1;
    }
    return rect;
}

// Genetic Algorithm
Rectangle geneticAlgorithm(Rectangle *smallRects, int numSmallRects, Bounds bounds) {
    Rectangle population[POPULATION_SIZE];
    Rectangle best = {0, 0, 0, 0, 0};

    // Initialize population
    for (int i = 0; i < POPULATION_SIZE; i++) {
        population[i] = generateRandomRectangle(bounds);
        population[i].fitness = calculateFitness(population[i], smallRects, numSmallRects, bounds);
    }

    // Evolve over generations
    for (int gen = 0; gen < GENERATIONS; gen++) {
        // Sort by fitness
        for (int i = 0; i < POPULATION_SIZE; i++) {
            for (int j = i + 1; j < POPULATION_SIZE; j++) {
                if (population[j].fitness > population[i].fitness) {
                    Rectangle temp = population[i];
                    population[i] = population[j];
                    population[j] = temp;
                }
            }
        }

        // Keep track of the best solution
        if (population[0].fitness > best.fitness) {
            best = population[0];
        }

        // Create new population
        Rectangle newPopulation[POPULATION_SIZE];
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Rectangle parent1 = population[rand() % (POPULATION_SIZE / 2)];
            Rectangle parent2 = population[rand() % (POPULATION_SIZE / 2)];
            Rectangle child = crossover(parent1, parent2);
            child = mutate(child, bounds);
            child.fitness = calculateFitness(child, smallRects, numSmallRects, bounds);
            newPopulation[i] = child;
        }

        // Replace old population
        for (int i = 0; i < POPULATION_SIZE; i++) {
            population[i] = newPopulation[i];
        }
    }

    return best;
}

int main() {
    srand(time(NULL));

    // Define bounds and obstacles
    Bounds bounds = {0, 0, 10, 8};
    Rectangle smallRects[] = {
        {1, 1, 3, 2},
        {4, 1, 6, 3},
        {7, 2, 9, 4}
    };
    int numSmallRects = 3;

    // Run genetic algorithm
    Rectangle best = geneticAlgorithm(smallRects, numSmallRects, bounds);

    printf("Best rectangle: Bottom-left (%.2f, %.2f), Top-right (%.2f, %.2f), Area: %.2f\n",
           best.x1, best.y1, best.x2, best.y2, best.fitness);

    return 0;
}

