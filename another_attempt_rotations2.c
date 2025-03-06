#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define POPULATION_SIZE 1000
#define MAX_GENERATIONS 100000
#define CONVERGENCE_THRESHOLD 500 // Generations without improvement before termination
#define INITIAL_MUTATION_RATE 0.1
#define LARGE_RECT_WIDTH 10.0
#define LARGE_RECT_HEIGHT 8.0
#define NUM_SMALL_RECTS 3

typedef struct {
    float x, y;
} Point;

typedef struct {
    Point center;
    float width, height;
    float angle; // In radians
} RotatedRectangle;

// Utility Functions
float random_float(float min, float max) {
    return min + ((float)rand() / RAND_MAX) * (max - min);
}

float rectangle_area(RotatedRectangle rect) {
    return rect.width * rect.height;
}

void project_onto_axis(Point vertices[], int num_vertices, Point axis, float *min_proj, float *max_proj) {
    *min_proj = *max_proj = (vertices[0].x * axis.x + vertices[0].y * axis.y);
    for (int i = 1; i < num_vertices; i++) {
        float projection = (vertices[i].x * axis.x + vertices[i].y * axis.y);
        if (projection < *min_proj) *min_proj = projection;
        if (projection > *max_proj) *max_proj = projection;
    }
}



void get_rectangle_corners(RotatedRectangle rect, Point corners[4]) {
    float cos_angle = cos(rect.angle);
    float sin_angle = sin(rect.angle);
    float half_width = rect.width / 2.0;
    float half_height = rect.height / 2.0;

    corners[0] = (Point){
        rect.center.x + (-half_width * cos_angle - -half_height * sin_angle),
        rect.center.y + (-half_width * sin_angle + -half_height * cos_angle)};
    corners[1] = (Point){
        rect.center.x + (half_width * cos_angle - -half_height * sin_angle),
        rect.center.y + (half_width * sin_angle + -half_height * cos_angle)};
    corners[2] = (Point){
        rect.center.x + (half_width * cos_angle - half_height * sin_angle),
        rect.center.y + (half_width * sin_angle + half_height * cos_angle)};
    corners[3] = (Point){
        rect.center.x + (-half_width * cos_angle - half_height * sin_angle),
        rect.center.y + (-half_width * sin_angle + half_height * cos_angle)};
}

bool rectangles_overlap(RotatedRectangle rect1, RotatedRectangle rect2) {
    Point axes[4];
    Point rect1_corners[4], rect2_corners[4];
    get_rectangle_corners(rect1, rect1_corners);
    get_rectangle_corners(rect2, rect2_corners);

    for (int i = 0; i < 2; i++) {
        int next = (i + 1) % 4;
        axes[i].x = -(rect1_corners[next].y - rect1_corners[i].y);
        axes[i].y = rect1_corners[next].x - rect1_corners[i].x;
        axes[i + 2].x = -(rect2_corners[next].y - rect2_corners[i].y);
        axes[i + 2].y = rect2_corners[next].x - rect2_corners[i].x;
    }

    for (int i = 0; i < 4; i++) {
        float min1, max1, min2, max2;
        project_onto_axis(rect1_corners, 4, axes[i], &min1, &max1);
        project_onto_axis(rect2_corners, 4, axes[i], &min2, &max2);

        if (max1 < min2 || max2 < min1) return false;
    }

    return true;
}



bool is_completely_within_large_rectangle(RotatedRectangle rect, RotatedRectangle large_rect) {
    // Simple check: ensure all corners are within the bounding box
    Point corners[4];
    get_rectangle_corners(rect, corners);
    for (int i = 0; i < 4; i++) {
        Point p = corners[i];
        if (p.x < large_rect.center.x - large_rect.width / 2.0 ||
            p.x > large_rect.center.x + large_rect.width / 2.0 ||
            p.y < large_rect.center.y - large_rect.height / 2.0 ||
            p.y > large_rect.center.y + large_rect.height / 2.0) {
            return false;
        }
    }
    return true;
}

float evaluate_fitness(RotatedRectangle rect, RotatedRectangle large_rect, RotatedRectangle small_rects[], int num_small_rects) {
    // Early rejection for obviously invalid rectangles
    if (!is_completely_within_large_rectangle(rect, large_rect)) {
        return 0;
    }

    // Check overlap with small rectangles
    for (int i = 0; i < num_small_rects; i++) {
        if (rectangles_overlap(rect, small_rects[i])) {
            return 0;
        }
    }
    return rectangle_area(rect);
}

RotatedRectangle tournament_selection(RotatedRectangle population[], float fitness[]) {
    int index1 = rand() % POPULATION_SIZE;
    int index2 = rand() % POPULATION_SIZE;

    return fitness[index1] > fitness[index2] ? population[index1] : population[index2];
}

RotatedRectangle mutate(RotatedRectangle rect, float mutation_rate) {
    if (random_float(0.0, 1.0) < mutation_rate) rect.width += random_float(-1.0, 1.0);
    if (random_float(0.0, 1.0) < mutation_rate) rect.height += random_float(-1.0, 1.0);
    if (random_float(0.0, 1.0) < mutation_rate) rect.angle += random_float(-M_PI / 8, M_PI / 8);
    if (random_float(0.0, 1.0) < mutation_rate) {
        rect.center.x += random_float(-1.0, 1.0);
        rect.center.y += random_float(-1.0, 1.0);
    }
    return rect;
}

int main() {
    srand(time(NULL));

    RotatedRectangle large_rect = {{5.0, 4.0}, LARGE_RECT_WIDTH, LARGE_RECT_HEIGHT, 0.0};
    RotatedRectangle small_rects[NUM_SMALL_RECTS] = {
        {{2.0, 2.0}, 2.0, 2.0, 0.785},
        {{7.0, 2.0}, 2.0, 2.0, 0.0},
        {{5.0, 6.0}, 2.0, 2.0, 0.0}};

    RotatedRectangle population[POPULATION_SIZE];
    float fitness[POPULATION_SIZE];

    for (int i = 0; i < POPULATION_SIZE; i++) {
        population[i] = (RotatedRectangle){
            {random_float(0.0, LARGE_RECT_WIDTH), random_float(0.0, LARGE_RECT_HEIGHT)},
            random_float(1.0, 5.0),
            random_float(1.0, 5.0),
            random_float(0.0, 2 * M_PI)};
    }

    RotatedRectangle best_solution;
    float best_fitness = 0;
    int generations_without_improvement = 0;
    float mutation_rate = INITIAL_MUTATION_RATE;

    for (int generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Evaluate fitness
        for (int i = 0; i < POPULATION_SIZE; i++) {
            fitness[i] = evaluate_fitness(population[i], large_rect, small_rects, NUM_SMALL_RECTS);
            if (fitness[i] > best_fitness) {
                best_solution = population[i];
                best_fitness = fitness[i];
                generations_without_improvement = 0; // Reset convergence counter
            }
        }

        // Check for convergence
        generations_without_improvement++;
        if (generations_without_improvement > CONVERGENCE_THRESHOLD) {
            printf("Converged after %d generations.\n", generation);
            break;
        }

        // Adjust mutation rate based on improvement
        mutation_rate = INITIAL_MUTATION_RATE / (1 + generations_without_improvement / 100.0);

        // Generate new population
        RotatedRectangle new_population[POPULATION_SIZE];
        for (int i = 0; i < POPULATION_SIZE; i++) {
            RotatedRectangle parent = tournament_selection(population, fitness);
            new_population[i] = mutate(parent, mutation_rate);
        }
        memcpy(population, new_population, sizeof(new_population));
    }

    printf("Best rectangle: Center (%.2f, %.2f), Width %.2f, Height %.2f, Angle %.2f radians\n",
           best_solution.center.x, best_solution.center.y,
           best_solution.width, best_solution.height,
           best_solution.angle);

    return 0;
}

