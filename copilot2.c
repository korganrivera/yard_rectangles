#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h> // For parallel processing

#define POPULATION_SIZE 1000
#define MAX_GENERATIONS 100000
#define CONVERGENCE_THRESHOLD 500 // Generations without improvement before termination
#define INITIAL_MUTATION_RATE 0.1
#define NUM_SMALL_RECTS 3

typedef struct {
    float x, y;
} Point;

typedef struct {
    Point vertices[4];
} IrregularQuadrilateral;

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

bool is_point_in_quadrilateral(Point p, IrregularQuadrilateral quad) {
    // Function to check if a point is inside an irregular quadrilateral
    // Using the crossing number algorithm or ray-casting method

    int i, j, crossing_number = 0;
    for (i = 0, j = 3; i < 4; j = i++) {
        if (((quad.vertices[i].y > p.y) != (quad.vertices[j].y > p.y)) &&
            (p.x < (quad.vertices[j].x - quad.vertices[i].x) * (p.y - quad.vertices[i].y) / (quad.vertices[j].y - quad.vertices[i].y) + quad.vertices[i].x)) {
            crossing_number++;
        }
    }
    return (crossing_number % 2 == 1);
}

bool is_completely_within_large_quadrilateral(RotatedRectangle rect, IrregularQuadrilateral large_quad) {
    // Check if all corners of the rectangle are within the quadrilateral
    Point corners[4];
    get_rectangle_corners(rect, corners);
    for (int i = 0; i < 4; i++) {
        if (!is_point_in_quadrilateral(corners[i], large_quad)) {
            return false;
        }
    }
    return true;
}

float evaluate_fitness(RotatedRectangle rect, IrregularQuadrilateral large_quad, RotatedRectangle small_rects[], int num_small_rects) {
    // Early rejection for obviously invalid rectangles
    if (!is_completely_within_large_quadrilateral(rect, large_quad)) {
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

void generate_svg(const char *filename, IrregularQuadrilateral large_quad,
                  RotatedRectangle small_rects[], int num_small_rects,
                  RotatedRectangle best_rect) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not create SVG file.\n");
        return;
    }

    // Start SVG
    fprintf(file, "<svg xmlns='http://www.w3.org/2000/svg' width='500' height='500' viewBox='-1 -1 12 10' style='background-color: white;'>\n");
    fprintf(file, "<style> .quad { fill: none; stroke: black; stroke-width: 0.05; } .small { fill: rgba(255, 0, 0, 0.5); stroke: red; stroke-width: 0.05; } .best { fill: rgba(0, 255, 0, 0.5); stroke: green; stroke-width: 0.05; } </style>\n");

    // Draw the large quadrilateral
    fprintf(file, "<polygon class='quad' points='");
    for (int i = 0; i < 4; i++) {
        fprintf(file, "%.2f,%.2f ", large_quad.vertices[i].x, large_quad.vertices[i].y);
    }
    fprintf(file, "' />\n");

    // Draw small rectangles
    for (int i = 0; i < num_small_rects; i++) {
        Point corners[4];
        get_rectangle_corners(small_rects[i], corners);
        fprintf(file, "<polygon class='small' points='");
        for (int j = 0; j < 4; j++) {
            fprintf(file, "%.2f,%.2f ", corners[j].x, corners[j].y);
        }
        fprintf(file, "' />\n");
    }

    // Draw the best rectangle
    Point best_corners[4];
    get_rectangle_corners(best_rect, best_corners);
    fprintf(file, "<polygon class='best' points='");
    for (int j = 0; j < 4; j++) {
        fprintf(file, "%.2f,%.2f ", best_corners[j].x, best_corners[j].y);
    }
    fprintf(file, "' />\n");

    // End SVG
    fprintf(file, "</svg>\n");
    fclose(file);

    printf("SVG file generated: %s\n", filename);
}





int main() {
    srand(time(NULL));

    IrregularQuadrilateral large_quad = {
        .vertices = {{0.0, 0.0}, {10.0, 0.0}, {8.0, 6.0}, {2.0, 8.0}}}; // Define the vertices of the large quadrilateral
    RotatedRectangle small_rects[NUM_SMALL_RECTS] = {
        {{2.0, 2.0}, 2.0, 2.0, 0.785},
        {{7.0, 2.0}, 2.0, 2.0, 0.0},
        {{5.0, 4.0}, 2.0, 2.0, 0.0}}; 

/*
    IrregularQuadrilateral large_quad = {
        .vertices = {{0, 0}, {79.224, 70.317}, {55, 129.33}, {0.55, 123.627}}}; // Define the vertices of the large quadrilateral
    RotatedRectangle small_rects[NUM_SMALL_RECTS] = {
        {{1316.451, 124.449}, 6.702, 6.702, 0.0},
        {{1292.505, 127.971}, 11.354, 11.354, 0.0},
        {{1273.741, 147.889}, 11.354, 11.354, 0.0},
        {{1303.002, 171.790}, 3.354, 3.354, 0.0},
        {{1253.76, 161.859}, 12.5, 14.17, 0.105},
        {{1238.202, 156.861}, 12.5, 7.33, 0.105},
        {{1287.673, 167.546}, 20, 13, 0.105}}; // coop
*/

    RotatedRectangle population[POPULATION_SIZE];
    float fitness[POPULATION_SIZE];

    for (int i = 0; i < POPULATION_SIZE; i++) {
        population[i] = (RotatedRectangle){
            {random_float(0.0, 10.0), random_float(0.0, 8.0)},
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
        #pragma omp parallel for
        for (int i = 0; i < POPULATION_SIZE; i++) {
            fitness[i] = evaluate_fitness(population[i], large_quad, small_rects, NUM_SMALL_RECTS);
            if (fitness[i] > best_fitness) {
                #pragma omp critical
                {
                    if (fitness[i] > best_fitness) {
                        best_solution = population[i];
                        best_fitness = fitness[i];
                        generations_without_improvement = 0; // Reset convergence counter
                    }
                }
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
        #pragma omp parallel for
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

    generate_svg("output.svg", large_quad, small_rects, NUM_SMALL_RECTS, best_solution);

    return 0;
}
